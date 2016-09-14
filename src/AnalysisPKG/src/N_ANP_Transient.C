//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2015 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_Transient.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.237.2.2 $
// Revision Date  : $Date: 2015/11/21 01:21:41 $
// Current Owner  : $Author: tmei $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>
#include <iomanip>

#include <N_ANP_Transient.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Report.h>
#include <N_ERH_Progress.h>
#include <N_IO_CmdParse.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_RestartMgr.h>
#include <N_LAS_System.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_MPDE_Manager.h>
#include <N_NLS_ConductanceExtractor.h>
#include <N_NLS_ReturnCodes.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_OneStep.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Param.h>
#include <N_UTL_Timer.h>

/* KIM*/
#include <iomanip>

namespace Xyce {
namespace Analysis {

void writeConductanceFile(const std::vector<std::string> &device_names, Nonlinear::ConductanceExtractor &conductance_extractor, const std::string &filename);

//-----------------------------------------------------------------------------
// Function      : Transient::Transient(AnalysisManager &manager)
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
Transient::Transient(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::RestartMgr &                      restart_manager,
  OutputAdapter *                       output_adapter,
  HB *                                  hb_analysis,
  N_MPDE_Manager *                      mpde_manager)
  : AnalysisBase(analysis_manager, "Transient"),
    StepEventListener(&analysis_manager),
    comm_(analysis_manager.getPDSManager()->getPDSComm()->comm()),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    restartManager_(restart_manager),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    outputAdapter_(output_adapter),
    tiaParams_(),
    sensFlag_(analysis_manager.getSensFlag()),
    initialIntegrationMethod_(TimeIntg::OneStep::type),
    firstTranOutput_(true),
    isPaused(false),
    dcopFlag_(true),
    startDCOPtime(0.0),
    startTRANtime_(0.0),
    endTRANtime_(0.0),
    gui_(analysisManager_.getCommandLine().argExists("-gui")),
    maxTimeStepExpressionString_(),
    maxTimeStepExpression_(0),
    firstTime(true),
    oldPercentComplete(0.0),
    startSimTime(-1.0),
    nextRestartSaveTime_(0.0),
    dcStats(0),
    tranStats(0),
    exitTime(0.0),
    exitStep(-1),
    integrationMethod(7),
    historyTrackingDepth(25),
    passNLStall(false),
    saveTimeStepsFlag(false),
    condTestFlag(false),
    condTestDeviceNames(),
    hbAnalysis_(hb_analysis),
    mpdeManager_(mpde_manager)
{}

void Transient::notify(const StepEvent &event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    nonlinearManager_.resetAll(Nonlinear::DC_OP);
    analysisManager_.getStepErrorControl().resetAll(tiaParams_);
    dcopFlag_ = false;
    analysisManager_.setNextOutputTime(0.0);

    if (timeQueue_.get_size() != 0)
    {
      // these set_size calls will also "reset" the queue.
      timeQueue_.reset_queue();
      timeStepQueue_.reset_queue();
      stepStatusQueue_.reset_queue();
      estErrorOverTolQueue_.reset_queue();
      nonlinearSolverStatusQueue_.reset_queue();
      nonlinearSolverNumIterationsQueue_.reset_queue();
      nonlinearSolverMaxNormQueue_.reset_queue();
      nonlinearSolverMaxNormIndexQueue_.reset_queue();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .TRAN statement.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/10
//-----------------------------------------------------------------------------
bool Transient::setAnalysisParams(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = *it;

    if (tiaParams_.setAnalysisOption(param))
      ;
    else if (param.uTag() == "NOOP" ||
             param.uTag() == "UIC")
    {
      setNOOP(true);
    }
    else if (param.uTag() == "MAXTIMEEXPRESSION")
    {
      // an expression was given which should specify what
      // max time step to use based on the current simulation
      // time.  The expected format of this is
      // {schedule( t0, max_del_t0, t1, max_del_t1, ...)}
      // we'll just store the expression as a string for now
      // as we need to make sure other classes are up and running
      // before we can create the actual expression.  That will
      // happen in Transient class as it will ultimately use this
      maxTimeStepExpressionString_ = param.stringValue();
    }
    else
      IO::ParamError(option_block, param) << "Unrecognized analysis option";
  }

  if (tiaParams_.finalTime <= tiaParams_.initialOutputTime || tiaParams_.finalTime <= 0.0 || tiaParams_.initialOutputTime < 0.0)
  {
    UserFatal(*this) << "Final time of " << tiaParams_.finalTime
                     << " is earlier or same as start time of " << tiaParams_.initialOutputTime << std::endl
                     << " Check netlist for invalid .TRAN specification";
  }

  // if starting time steps is baloney, set then to a default value.
  // ERK.  This trap is redudant with an identical one in step error
  // control.
  if ( tiaParams_.initialTimeStep <= 0.0 )
  {
    tiaParams_.initialTimeStep = 1.0e-10;
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << std::endl
           << section_divider << std::endl
           << "  Transient simulation parameters" << std::endl
           << "  initial time = " << tiaParams_.initialTime << std::endl
           << "  final time   = " << tiaParams_.finalTime << std::endl
           << "  starting time step = " << tiaParams_.initialTimeStep << std::endl
           << "  initialOutputTime (time of first output) = " << tiaParams_.initialOutputTime << std::endl;

    if (!getNOOP())
    {
      dout() << "  NOOP/UIC is NOT set" << std::endl;
    }
    else
    {
      dout() << "  NOOP/UIC is set" << std::endl;
    }

    dout() << section_divider << std::endl;
  }

  if (!maxTimeStepExpressionString_.empty())
  {
    // set-up expression object for getting a user specified, time dependent max time step value
    delete maxTimeStepExpression_;

    maxTimeStepExpression_ = new Util::ExpressionData(maxTimeStepExpressionString_);
  }

  // historySize_ should can be set from the .options timeint line.  Use default in tiaParams
  // as that will be overwritten if the user set it via .options.  Also, if it's less than
  // or equal to zero assume that the user wants this turned off.
  // We'll store the attempted time step, time step status,
  // estimated error over tolerance, non-linear solver status and non-linear solver norm
  // for each step, up to the last historySize_ steps.

  if (historyTrackingDepth > 0 )
  {
    int size = historyTrackingDepth;
    timeQueue_.set_size( size );
    timeStepQueue_.set_size( size );
    stepStatusQueue_.set_size( size );
    estErrorOverTolQueue_.set_size( size );
    nonlinearSolverStatusQueue_.set_size( size );
    nonlinearSolverNumIterationsQueue_.set_size( size );
    nonlinearSolverMaxNormQueue_.set_size( size );
    nonlinearSolverMaxNormIndexQueue_.set_size( size );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool Transient::setTimeIntegratorOptions(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "METHOD")
    {
      if (param.isInteger())
        integrationMethod = param.getImmutableValue<int>();
      else
      {
        ExtendedString stringVal ( param.stringValue() );
        stringVal.toUpper();

        if (stringVal == "TRAP" || stringVal == "TRAPEZOIDAL")
          integrationMethod = 7;
        else if (stringVal == "BDF")
          integrationMethod = 6;
        else if (stringVal == "GEAR")
          integrationMethod = 8;
        else
        {
          IO::ParamError(option_block, param) << "Unsupported transient method type";
        }
      }

    }
    else if (param.uTag()=="EXITTIME" )
      exitTime = param.getImmutableValue<double>();
    else if (param.uTag()=="EXITSTEP" )
      exitStep = param.getImmutableValue<int>();
    else if (param.uTag() == "HISTORYTRACKINGDEPTH" )
      historyTrackingDepth = param.getImmutableValue<int>();
    else if (param.uTag() == "PASSNLSTALL")
      passNLStall = param.getImmutableValue<bool>();
    else if (param.uTag() == "CONDTEST")
      condTestFlag = static_cast<bool> (param.getImmutableValue<int>());
    else if (param.uTag() == "CONDTESTDEVICENAME")
      condTestDeviceNames.push_back(param.stringValue() );
    else if (param.uTag() == "DAESTATEDERIV" )
      analysisManager_.setDAEStateDerivFlag(static_cast<bool> (param.getImmutableValue<int>()));
    else if (param.uTag() == "DEBUGLEVEL" )
      IO::setTimeIntegratorDebugLevel(analysisManager_.getCommandLine(), param.getImmutableValue<int>());
    else if (nonlinearManager_.setReturnCodeOption(param))
      ;
    else if (tiaParams_.setTimeIntegratorOption(param))
      ;
    else if (setDCOPOption(param))
      ;
    else
      IO::ParamError(option_block, param) << "Not a recognized time integration option";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Transient::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::doRun()
{
  isPaused = false;

  bool bsuccess = false;

  if (analysisManager_.getResumingSimulation())
  {
    bsuccess   = resuming() && doLoopProcess() && doFinish();
  }
  else
  {
    bsuccess   = doInit() && doTranOP () && doLoopProcess() && doFinish();
  }

  if (condTestFlag)
  {
    writeConductanceFile(condTestDeviceNames, nonlinearManager_.getConductanceExtractor(), "conductance.txt");
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::resuming()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::resuming()
{
  bool bsuccess = true;

  if (sensFlag_)
  {
    Stats::StatTop _sensitivityStat("Sensitivity");

    nonlinearManager_.enableSensitivity(*analysisManager_.getDataStore(), *analysisManager_.getPDSManager(), topology_);
  }

  initialIntegrationMethod_ = integrationMethod != TimeIntg::OneStep::type ? integrationMethod : TimeIntg::OneStep::type;

  if (analysisManager_.getTwoLevelMode() == TWO_LEVEL_MODE_TRANSIENT)
  {
    // we're resuming from a previous transient simulation
    // this is an ugly side effect of the analysis manager's deleting
    // the transient object after it returned in a paused state
    // need to fix this in Xyce 6.0 .  RLS 12/22/2009
    dcopFlag_ = false;

    // we default firstTranOutput_ to true to force the first time point to be output
    // however, if we're resuming, then we need this to be false
    // this is really a defect of a paused simulation's Transient object
    // being deleted by the AnalysisManager and then a new one being created
    // when the simulation is resumed.  Need to fully fix this.  RLS 12/22/2009
    firstTranOutput_ = false;
  }

  // in a 2 level problem, we can be asked to resume
  // without first making the integration method.  So
  // catch that.
  if (analysisManager_.getWorkingIntegrationMethod().isTimeIntegrationMethodCreated())
  {
    baseIntegrationMethod_ = integrationMethod;
    analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    dout() << "  transient init called with resume true " << std::endl;
  analysisManager_.setSwitchIntegrator(false);

  // reset min error tracking variables
  // if the step number is less than zero we'll assume there is no valid
  // min. estimated error over tol or an associated time step.  This frees us
  // from putting Machine::Big here and trying to do something with that.
  stepNumberAtMinEstErrorOverTol = -1;
  minEstErrorOverTol = 0.0;
  timeStepAtMinEstErrorOverTol = 0.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::doInit()
{
  bool bsuccess = true;


  if (sensFlag_)
  {
    Stats::StatTop _sensitivityStat("Sensitivity");

    nonlinearManager_.enableSensitivity(*analysisManager_.getDataStore(), *analysisManager_.getPDSManager(), topology_);
  }

  initialIntegrationMethod_ = integrationMethod != TimeIntg::OneStep::type ? integrationMethod : TimeIntg::OneStep::type;

  dcopFlag_ = !getNOOP();

  analysisManager_.setTwoLevelMode(dcopFlag_ ? Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP : Analysis::TWO_LEVEL_MODE_TRANSIENT);

  setDoubleDCOPEnabled(loader_.isPDESystem());

  if (restartManager_.isRestarting())
  {
    loader_.getInitialQnorm(analysisManager_.getDataStore()->innerErrorInfoVec);

    analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
    analysisManager_.getWorkingIntegrationMethod().initialize(tiaParams_);

    dcopFlag_ = false;
    analysisManager_.setTwoLevelMode(Analysis::TWO_LEVEL_MODE_TRANSIENT);

    // Update vectors with off proc values.
    linearSystem_.updateExternValsSolnVector(analysisManager_.getDataStore()->nextSolutionPtr);
    linearSystem_.updateExternValsSolnVector(analysisManager_.getDataStore()->currSolutionPtr);
    linearSystem_.updateExternValsStateVector(analysisManager_.getDataStore()->nextStatePtr);
    linearSystem_.updateExternValsStateVector(analysisManager_.getDataStore()->currStatePtr);
    linearSystem_.updateExternValsStoreVector(analysisManager_.getDataStore()->nextStorePtr);
    linearSystem_.updateExternValsStoreVector(analysisManager_.getDataStore()->currStorePtr);

    // Set the nonlinear solver parameters to those appropriate for the
    // transient solution.  If starting from here, we're not doing a DCOP
    // calculation first - we're jumping right into the transient.
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_TRANSIENT));
  }
  else // not restart
  {
    if (dcopFlag_)
    {
      // Get set to do the operating point.
      baseIntegrationMethod_ = TimeIntg::NoTimeIntegration::type;
      analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
    }
    else // otherwise NOOP/UIC
    {
      baseIntegrationMethod_ = initialIntegrationMethod_;
      analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
    }

    // This setInitialGuess call is to up an initial guess in the
    // devices that have them (usually PDE devices).  This is DIFFERENT
    // than an initial condition.
    loader_.setInitialGuess(analysisManager_.getDataStore()->nextSolutionPtr);
    /* SOLNnot populated here */

    // 12/13/07 tscoffe/tmei:  .ic should not be applied to MPDE currently.
    // A new mechanism must be set up for this, e.g. .mpde_ic
    if (!mpdeManager_ || mpdeManager_->getTransientNeedsToLoadInitialConditionsAndInitializeProblem())
    {
      // If available, set initial solution.  This may be from a DCOP restart,
      // a .IC, or a .NODESET line.  These can be used in both the UIC/NOOP
      // case, as well as the DCOP case, so they need to be outside the if-statement.
      setInputOPFlag(
        initialConditionsManager_.setupInitialConditions(outputManagerAdapter_.getComm(),
                                                         topology_.getSolutionNodeNameMap(),
                                                         *analysisManager_.getDataStore()->nextSolutionPtr,
                                                         *analysisManager_.getDataStore()->flagSolutionPtr));

      if (!dcopFlag_)
      {
        // this "initializeProblem" call is to set the IC's in devices that
        // have them.  This is done IN PLACE of the operating point.
        loader_.initializeProblem
          ((analysisManager_.getDataStore()->nextSolutionPtr),
           (analysisManager_.getDataStore()->currSolutionPtr),
           (analysisManager_.getDataStore()->lastSolutionPtr),
           (analysisManager_.getDataStore()->nextStatePtr),
           (analysisManager_.getDataStore()->currStatePtr),
           (analysisManager_.getDataStore()->lastStatePtr),
           (analysisManager_.getDataStore()->nextStateDerivPtr),
           (analysisManager_.getDataStore()->nextStorePtr),
           (analysisManager_.getDataStore()->currStorePtr),
           (analysisManager_.getDataStore()->lastStorePtr),
           (analysisManager_.getDataStore()->daeQVectorPtr),
           (analysisManager_.getDataStore()->daeFVectorPtr),
           (analysisManager_.getDataStore()->daeBVectorPtr),
           (analysisManager_.getDataStore()->dFdxdVpVectorPtr),
           (analysisManager_.getDataStore()->dQdxdVpVectorPtr) );

        // Do this to populate the q-vector:
        // since we're also skipping the DC OP, this call will
        // also force the loader to propagate the state vector to the
        // device manager which updates state vector seen by the devices.
        analysisManager_.getNonlinearEquationLoader().loadRHS();
      }
    }

    // Set a constant history.
    analysisManager_.getDataStore()->setConstantHistory();
    analysisManager_.getDataStore()->computeDividedDifferences();
    analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();

    // Update vectors with off proc values.
    linearSystem_.updateExternValsSolnVector(analysisManager_.getDataStore()->nextSolutionPtr);
    linearSystem_.updateExternValsSolnVector(analysisManager_.getDataStore()->currSolutionPtr);
    linearSystem_.updateExternValsStateVector(analysisManager_.getDataStore()->nextStatePtr);
    linearSystem_.updateExternValsStateVector(analysisManager_.getDataStore()->currStatePtr);
    linearSystem_.updateExternValsStoreVector(analysisManager_.getDataStore()->nextStorePtr);
    linearSystem_.updateExternValsStoreVector(analysisManager_.getDataStore()->currStorePtr);

    // if we are skipping the DCOP, but .sens was requested, then call this here.
    if (!dcopFlag_ && !mpdeManager_ && sensFlag_)
    {
      nonlinearManager_.icSensitivity (objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

      analysisManager_.getDataStore()->setConstantSensitivityHistory();
    }

    if (!dcopFlag_)
    {
      // Now do a NOOP output.  ERK: 08/30/2007  This isn't the best place to
      // put this, but it will have to do for now.  If this isn't here, then
      // NOOP/UIC simulations don't output at t=0.0 in the *prn file.
      noopOutputs ();
    }

    stepNumber            = 0;
    tranStepNumber        = 0;
    analysisManager_.breakPointRestartStep = 0;
  }
  analysisManager_.setSwitchIntegrator(false);
  // if we're resuming, this was already done prior to pausing, and doing
  // it again screws us up
  double suggestedMaxTime = 0.0;
  if (maxTimeStepExpression_)
  {
    if (maxTimeStepExpression_->setup(comm_,
                                      outputManagerAdapter_.getOutputManager().getOpBuilderManager(),
                                      outputManagerAdapter_.getOutputManager().getMainContextFunctionMap(),
                                      outputManagerAdapter_.getOutputManager().getMainContextParamMap(),
                                      outputManagerAdapter_.getOutputManager().getMainContextGlobalParamMap()) == Util::ExpressionData::READY)
      suggestedMaxTime = maxTimeStepExpression_->evaluate(comm_,
                                                          outputManagerAdapter_.getOutputManager().getCircuitTime(),
                                                          analysisManager_.getDataStore()->currSolutionPtr,
                                                          analysisManager_.getDataStore()->currStatePtr,
                                                          analysisManager_.getDataStore()->currStorePtr);
  }
  analysisManager_.getStepErrorControl().updateMaxTimeStep( comm_, loader_, tiaParams_, suggestedMaxTime );
  analysisManager_.getStepErrorControl().updateMinTimeStep();
  analysisManager_.getStepErrorControl().updateBreakPoints(loader_, tiaParams_.initialTime);

  nextRestartSaveTime_ = analysisManager_.getStepErrorControl().initialTime;

  // reset min error tracking variables
  // if the step number is less than zero we'll assume there is no valid
  // min. estimated error over tol or an associated time step.  This frees us
  // from putting Machine::Big here and trying to do something with that.
  stepNumberAtMinEstErrorOverTol = -1;
  minEstErrorOverTol = 0.0;
  timeStepAtMinEstErrorOverTol = 0.0;

  /* SOLN not populated here */
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::doTranOP ()
// Purpose       : Computes the DCOP calcualtion that precedes transient loop.
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 
//-----------------------------------------------------------------------------
bool Transient::doTranOP ()
{
  bool bsuccess = true;

  Stats::Stat _transientStat(Stats::StatTop::getTop());
  Stats::TimeBlock _transientTimer(_transientStat);

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::TRAN));

  // Transient time stepping loop:
  while (dcopFlag_)
  {
    if (VERBOSE_TIME)
      printStepHeader(Xyce::lout());

    printProgress(Xyce::lout());

    // ------------------------------------------------------------------------
    // If the flag is set to switch integration methods, do that here.
    // For example, switch from operating point to transient backward euler.
    if (analysisManager_.getSwitchIntegrator())
    {
      analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
    }

    // ------------------------------------------------------------------------
    // Set the step size, current time and next time.

    analysisManager_.getStepErrorControl().updateStopTime(
      comm_,
      tiaParams_.bpEnable,
      tiaParams_.initialTime,
      tiaParams_.minTimeStepsBPGiven,
      tiaParams_.minTimeStepsBP);

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << std::endl;
      dout() << "Transient::loopProcess()" << std::endl;
      dout() << "beginningIntegration = " << beginningIntegration << std::endl;
      dout() << "analysisManager_.getStepErrorControl().stepAttemptStatus = " << analysisManager_.getStepErrorControl().stepAttemptStatus << std::endl;
    }

    if (beginningIntegration &&
        analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      // ------------------------------------------------------------------------
      // 07/29/04 TSC:  initial step-size selection is now done in the
      // integration method initialize call under the DAE formulation.  This
      // segregates codes changes better and makes sense to do as an
      // initialization step even if its changed later.
      loader_.getInitialQnorm(analysisManager_.getDataStore()->innerErrorInfoVec);
      double suggestedMaxTime=0.0;
      if (maxTimeStepExpression_)
      {
        suggestedMaxTime = maxTimeStepExpression_->evaluate(comm_, outputManagerAdapter_.getOutputManager().getCircuitTime(),
                                                            analysisManager_.getDataStore()->currSolutionPtr, analysisManager_.getDataStore()->currStatePtr, analysisManager_.getDataStore()->currStorePtr);
      }
      analysisManager_.getStepErrorControl().updateMaxTimeStep( comm_, loader_, tiaParams_, suggestedMaxTime );
      analysisManager_.getWorkingIntegrationMethod().initialize(tiaParams_);
    }

    // ------------------------------------------------------------------------
    // If we've switched the integration method, we need to obtain the
    // corrector derivative only after we've updated the TimeInfo.
    if (analysisManager_.getSwitchIntegrator())
    {
      analysisManager_.setSwitchIntegrator(false);
      analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();
    }

    // Ask the method to update its coefficients
    analysisManager_.getWorkingIntegrationMethod().updateCoeffs();

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::TRAN, analysisManager_.getStepErrorControl().nextTime, analysisManager_.getStepErrorControl().getNumberOfSteps()));
    
    // ------------------------------------------------------------------------
    // Perform the time step:
    takeAnIntegrationStep_();

    // ------------------------------------------------------------------------
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      processSuccessfulDCOP(); /* NOTE KIM X is printed to curX */
    }
    else if (!analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      processFailedDCOP();
      bsuccess = false;
      break;
    }
  }

  endTRANtime_ = analysisManager_.getXyceTranTimer().elapsedTime();
  tranStats = saveLoopInfo();

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::TRAN));

  /* Soln populated here */
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Transient::doLoopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
void writeMatrix(int index, Linear::Matrix * mat, std::string name) {
        std::stringstream sf;
        sf << name;
        sf << "/";
        sf << name;
        sf << "_";
        sf << std::setfill('0');
        sf << std::setw(5);
        sf << index;
        sf << ".csv";
        std::string filename = sf.str();

        std::ofstream os;
        os.open(filename.c_str());
        mat->printPetraObject(os);
        os.close();
}

void writeMatrix(int index, Linear::Vector * mat, std::string name) {
        std::stringstream sf;
        sf << name;
        sf << "/";
        sf << name;
        sf << "_";
        sf << std::setfill('0');
        sf << std::setw(5);
        sf << index;
        sf << ".csv";
        std::string filename = sf.str();

        std::ofstream os;
        os.open(filename.c_str());
        mat->printPetraObject(os);
        os.close();
}

void print(Epetra_MultiVector *vec, char * syntax, int disambiguation) {
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        printf(syntax, i, 0, values[i]);
}


bool Transient::doLoopProcess()
{
  int index = 0;
  bool first = true;
  bool bsuccess = true;

  Stats::Stat _transientStat(Stats::StatTop::getTop());
  Stats::TimeBlock _transientTimer(_transientStat);

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::TRAN));

  // Transient time stepping loop:
  while (!(analysisManager_.getStepErrorControl().isFinished()))
  {
    if (VERBOSE_TIME)
      printStepHeader(Xyce::lout());

    printProgress(Xyce::lout());

    // ------------------------------------------------------------------------
    // If the flag is set to switch integration methods, do that here.
    // For example, switch from operating point to transient backward euler.

    if (analysisManager_.getSwitchIntegrator()) {
      analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
    }

    // ------------------------------------------------------------------------
    // Set the step size, current time and next time.

    analysisManager_.getStepErrorControl().updateStopTime(
      comm_,
      tiaParams_.bpEnable,
      tiaParams_.initialTime,
      tiaParams_.minTimeStepsBPGiven,
      tiaParams_.minTimeStepsBP);

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << std::endl;
      dout() << "Transient::loopProcess()" << std::endl;
      dout() << "beginningIntegration = " << beginningIntegration << std::endl;
      dout() << "analysisManager_.getStepErrorControl().stepAttemptStatus = " << analysisManager_.getStepErrorControl().stepAttemptStatus << std::endl;
    }

    if (beginningIntegration && analysisManager_.getStepErrorControl().stepAttemptStatus) {
      // ------------------------------------------------------------------------
      // 07/29/04 TSC:  initial step-size selection is now done in the
      // integration method initialize call under the DAE formulation.  This
      // segregates codes changes better and makes sense to do as an
      // initialization step even if its changed later.
      loader_.getInitialQnorm(analysisManager_.getDataStore()->innerErrorInfoVec);
      double suggestedMaxTime=0.0;
      if (maxTimeStepExpression_)
      {
        suggestedMaxTime = maxTimeStepExpression_->evaluate(comm_, outputManagerAdapter_.getOutputManager().getCircuitTime(),
                                                            analysisManager_.getDataStore()->currSolutionPtr, analysisManager_.getDataStore()->currStatePtr, analysisManager_.getDataStore()->currStorePtr);
      }
      analysisManager_.getStepErrorControl().updateMaxTimeStep( comm_, loader_, tiaParams_, suggestedMaxTime );
      analysisManager_.getWorkingIntegrationMethod().initialize(tiaParams_);
    }

    // ------------------------------------------------------------------------
    // If we've switched the integration method, we need to obtain the
    // corrector derivative only after we've updated the TimeInfo.
    if (analysisManager_.getSwitchIntegrator())
    {
      analysisManager_.setSwitchIntegrator(false);
      analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();
    }

    if (VERBOSE_TIME)
      analysisManager_.getStepErrorControl().outputTimeInfo(lout());

    // ------------------------------------------------------------------------
    // Set the nonlinear solver parameters to those appropriate for the
    // transient solution, if neccessary.
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_TRANSIENT));

    // Ask the method to update its coefficients
    analysisManager_.getWorkingIntegrationMethod().updateCoeffs();

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::TRAN, analysisManager_.getStepErrorControl().nextTime, analysisManager_.getStepErrorControl().getNumberOfSteps()));

    // ------------------------------------------------------------------------
    // Perform the time step:

    /* EDIT KIM */
    if(first) {
        //first = false;

        TimeIntg::DataStore * ds = analysisManager_.getDataStore();
        Linear::Matrix * dFdx = ds->dFdxMatrixPtr;
        writeMatrix(index, dFdx, "dFdx");
        std::cout << "mat: dFdx\n";
        //dFdx->printPetraObject(std::cout);

        Linear::Matrix * dQdx = ds->dQdxMatrixPtr;
        writeMatrix(index, dQdx, "dQdx");
     //   std::cout << "mat: dQdx\n";
     //   dQdx->printPetraObject(std::cout);
        Linear::Vector * rhs  = ds->RHSVectorPtr;
        writeMatrix(index, rhs, "rhs");
     //   std::cout << "mat: rhs\n";
     //   rhs->printPetraObject(std::cout);

        Linear::Vector * soln = ds->nextSolutionPtr;
        writeMatrix(index, soln, "soln"); /* SOLN NOT CHANGED */
        print(&soln->epetraObj(), "%d,%d,%0.64e\n", 0);
        std::cout << "done\n";
        std::cout << "done\n";
     //   std::cout << "mat: X\n";
      //  soln->printPetraObject(std::cout); /* Represented in prn */

        ++index;
    }
    /* END KIM */

    takeAnIntegrationStep_();
    

    // ------------------------------------------------------------------------
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      processSuccessfulStep(); /* no  correction to soln */
      /* EDIT KIM: print matrices here */
 //     TimeIntg::DataStore * ds = analysisManager_.getDataStore();
 //     Linear::Matrix * dFdx = ds->dFdxMatrixPtr;

 //     writeMatrix(index, dFdx, "dFdx");

 //     Linear::Matrix * dQdx = ds->dQdxMatrixPtr;
 //     writeMatrix(index, dQdx, "dQdx");

 //     Linear::Vector * rhs  = ds->RHSVectorPtr;
 //     writeMatrix(index, rhs, "rhs");
 //     rhs->printPetraObject(std::cout);
      
      TimeIntg::DataStore * ds = analysisManager_.getDataStore();
      Linear::Vector * soln = ds->nextSolutionPtr;
      soln->printPetraObject(std::cout);
   //   writeMatrix(index, soln, "soln");
 //       
 //     std::cout << "mat: X\n";
 //     soln->printPetraObject(std::cout);

      ++index;
      /* DONE */
    }
    else if (passNLStall
             && !analysisManager_.getStepErrorControl().stepAttemptStatus
             && (analysisManager_.getStepErrorControl().currentTimeStep < (4*analysisManager_.getStepErrorControl().minTimeStep)))
    {
      // potentially a VERY dangerous options.
      // if the non-linear solver is stalling, and we're very close to a min
      // time step, then calls this failure a pass
      if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -3)
      {
        {
          UserWarning(*this) << "Nonlinear solver stalled. Calling this a pass";
        }

        doProcessSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is reporting too big of an update, and we're very close to a min
      // time step, then calls this failure a pass
      if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -2)
      {
        {
          UserWarning(*this) << "Update too big. Calling this a pass";
        }

        doProcessSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is not converging in the max number of steps,
      // and we're very close to a min time step, then calls this failure a pass
      //
      // if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -1)
      // {
      //   Report::UserWarning0() << "Too many steps, calling this a pass";
      //   processSuccessfulStep();
      // }
      else
      {
        // process this failed step as we would have by default.
        bool b1 = doProcessFailedStep();
        if (!b1)
        {
          bsuccess = false;
          break;
        }
      }
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      bool b1 = doProcessFailedStep();
      if (!b1)
      {
        bsuccess = false;
        break;
      }
    } // stepAttemptStatus

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << std::endl;
      dout() << "   Here we are, just before checking whether to pause. " << std::endl;
      dout() << "   minTimeStep = " << analysisManager_.getStepErrorControl().minTimeStep << std::endl;
      dout() << "   final time = " << tiaParams_.finalTime << std::endl;
      dout() << "   pause time = " << analysisManager_.getPauseTime() << std::endl;
      dout() << "   initial time = " << analysisManager_.getStepErrorControl().initialTime << std::endl;
      dout() << "   current time = " << analysisManager_.getStepErrorControl().currentTime << std::endl;
      if (analysisManager_.getPauseTime() == analysisManager_.getStepErrorControl().currentTime)
      {
        dout() << "    Pause time and current time equal " << std::endl;
      }
      else
      {
        dout() << "     difference between current and pause times is "
               << analysisManager_.getPauseTime() - analysisManager_.getStepErrorControl().currentTime << std::endl;
      }
      if (analysisManager_.getPauseTime() == analysisManager_.getStepErrorControl().initialTime)
      {
        dout() << "    Pause time and initial time equal " << std::endl;
      }
    }

    if (analysisManager_.getStepErrorControl().isPauseTime())
    {
      // Failure at this point only indicates that the simulation
      // is paused and may be resumed.
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        dout() << "Transient::loopProcess():   pausing simulation " << std::endl;
      }

      analysisManager_.getStepErrorControl().simulationPaused(tiaParams_.initialTime);
      isPaused = true;
      bsuccess = true;
      break;
    }

    // If the exit time has been exceeded, exit.
    if (exitTime != 0.0 && analysisManager_.getStepErrorControl().currentTime > exitTime)
    {
      lout() << "Exit time exceeded.  Exiting transient loop\n" << std::endl;
      bsuccess = true;
      break;
    }

    if (exitStep != -1 && static_cast<int>(stepNumber) == exitStep)
    {
      lout() << "Exit step.  Exiting transient loop\n" << std::endl;
      bsuccess = true;
      break;
    }

  } // end of time loop

  endTRANtime_ = analysisManager_.getXyceTranTimer().elapsedTime();
  tranStats = saveLoopInfo();

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::TRAN));

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Transient::mixedSignalStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool Transient::mixedSignalStep(double maxTimeStepFromHabanero)
{
  preMixedSignalStepDetails(maxTimeStepFromHabanero);

  takeAnIntegrationStep_();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Transient::preMixedSignalStepDetails
// Purpose       :
// Special Notes : Habanero API function
// Scope         : private_
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
void Transient::preMixedSignalStepDetails(
  double        maxTimeStepFromHabanero)
{
  if (VERBOSE_TIME)
    printStepHeader(Xyce::lout());

  printProgress(Xyce::lout());

  // ------------------------------------------------------------------------
  // If the flag is set to switch integration methods, do that here.
  // For example, switch from operating point to transient backward euler.

  if (analysisManager_.getSwitchIntegrator())
  {
    analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);
  }

  // ------------------------------------------------------------------------
  // Set the step size, current time and next time.
  analysisManager_.getStepErrorControl().updateStopTime(
    comm_,
    tiaParams_.bpEnable,
    tiaParams_.initialTime,
    tiaParams_.minTimeStepsBPGiven,
    tiaParams_.minTimeStepsBP);

  // ------------------------------------------------------------------------
  // If a max time is set from Habanero, impose it now that stopTime is updated.
  if (maxTimeStepFromHabanero > 0)
  {
    double currentTimeStep = std::min(maxTimeStepFromHabanero, analysisManager_.getStepErrorControl().currentTimeStep);
    analysisManager_.getStepErrorControl().setTimeStep(currentTimeStep);
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << std::endl;
    dout() << "Transient::loopProcess()" << std::endl;
    dout() << "beginningIntegration = " << beginningIntegration << std::endl;
    dout() << "analysisManager_.getStepErrorControl().stepAttemptStatus = " << analysisManager_.getStepErrorControl().stepAttemptStatus << std::endl;
  }

  if (beginningIntegration &&
      analysisManager_.getStepErrorControl().stepAttemptStatus)
  {
    // ------------------------------------------------------------------------
    // 07/29/04 TSC:  initial step-size selection is now done in the
    // integration method initialize call under the DAE formulation.  This
    // segregates codes changes better and makes sense to do as an
    // initialization step even if its changed later.
    loader_.getInitialQnorm(analysisManager_.getDataStore()->innerErrorInfoVec);
    double suggestedMaxTime=0.0;
    if (maxTimeStepExpression_)
    {
      suggestedMaxTime = maxTimeStepExpression_->evaluate(comm_, outputManagerAdapter_.getOutputManager().getCircuitTime(),
                                                          analysisManager_.getDataStore()->currSolutionPtr, analysisManager_.getDataStore()->currStatePtr, analysisManager_.getDataStore()->currStorePtr);
    }
    analysisManager_.getStepErrorControl().updateMaxTimeStep( comm_, loader_, tiaParams_, suggestedMaxTime );
    analysisManager_.getWorkingIntegrationMethod().initialize(tiaParams_);
  }

  // ------------------------------------------------------------------------
  // If we've switched the integration method, we need to obtain the
  // corrector derivative only after we've updated the TimeInfo.
  if (analysisManager_.getSwitchIntegrator())
  {
    analysisManager_.setSwitchIntegrator(false);
    analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();
  }

  if (VERBOSE_TIME &&!dcopFlag_)
    analysisManager_.getStepErrorControl().outputTimeInfo(lout());

  // ------------------------------------------------------------------------
  // Set the nonlinear solver parameters to those appropriate for the
  // transient solution, if neccessary.
  if (!dcopFlag_)
  {
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_TRANSIENT));
  }

  // Ask the method to update its coefficients
  analysisManager_.getWorkingIntegrationMethod().updateCoeffs();
}

//-----------------------------------------------------------------------------
// Function      : Transient::finalizeMixedSignalStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool Transient::finalizeMixedSignalStep()
{
  bool recoverableFailureFlag = true;

  if (dcopFlag_ && analysisManager_.getStepErrorControl().stepAttemptStatus)
  {
    processSuccessfulDCOP();
  }
  else if (dcopFlag_ && !analysisManager_.getStepErrorControl().stepAttemptStatus)
  {
    processFailedDCOP();
    recoverableFailureFlag = false;
  }
  // Transient
  else
  {
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      doProcessSuccessfulStep();
    }
    else if (passNLStall
             && !analysisManager_.getStepErrorControl().stepAttemptStatus
             && (analysisManager_.getStepErrorControl().currentTimeStep < (4*analysisManager_.getStepErrorControl().minTimeStep)))
    {
      // potentially a VERY dangerous options.
      // if the non-linear solver is stalling, and we're very close to a min
      // time step, then calls this failure a pass

      if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -3)
      {
        Report::UserWarning0() << "Nonlinear solver stalled, calling this a pass";
        doProcessSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is reporting too big of an update, and we're very close to a min
      // time step, then calls this failure a pass
      if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -2)
      {
        Report::UserWarning0() << "Update too big, calling this a pass";
        doProcessSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is not converging in the max number of steps,
      // and we're very close to a min time step, then calls this failure a pass
      //
      // if( analysisManager_.getStepErrorControl().newtonConvergenceStatus == -1)
      // {
      //   Report::UserWarning0() << "?Too many steps, calling this a pass";
      //   processSuccessfulStep();
      // }
      else
      {
        // process this failed step as we would have by default.
        recoverableFailureFlag = doProcessFailedStep();
      }
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      recoverableFailureFlag = doProcessFailedStep();
    }
  } // transient

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << std::endl;
    dout() << "   Here we are, just before checking whether to pause. " << std::endl;
    dout() << "   minTimeStep = " << analysisManager_.getStepErrorControl().minTimeStep << std::endl;
    dout() << "   final time = " << tiaParams_.finalTime << std::endl;
    dout() << "   pause time = " << analysisManager_.getPauseTime() << std::endl;
    dout() << "   initial time = " << analysisManager_.getStepErrorControl().initialTime << std::endl;
    dout() << "   current time = " << analysisManager_.getStepErrorControl().currentTime << std::endl;
    if (analysisManager_.getPauseTime() == analysisManager_.getStepErrorControl().currentTime)
    {
      dout() << "    Pause time and current time equal " << std::endl;
    }
    else
    {
      dout() << "     difference between current and pause times is " << analysisManager_.getPauseTime() - analysisManager_.getStepErrorControl().currentTime << std::endl;
    }
    if (analysisManager_.getPauseTime() == analysisManager_.getStepErrorControl().initialTime)
    {
      dout() << "    Pause time and initial time equal " << std::endl;
    }
  }

  if (analysisManager_.getStepErrorControl().isPauseTime())
  {
    // Failure at this point only indicates that the simulation
    // is paused and may be resumed.
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << "Transient::loopProcess():   pausing simulation " << std::endl;
    }

    analysisManager_.getStepErrorControl().simulationPaused(tiaParams_.initialTime);
    isPaused = true;
    recoverableFailureFlag = false;
  }

  // If the exit time has been exceeded, exit.
  if (exitTime != 0.0 && analysisManager_.getStepErrorControl().currentTime > exitTime)
  {

    lout() << "Exit time exceeded.  Exiting transient loop\n" << std::endl;
    recoverableFailureFlag = false;
  }

  if (exitStep != -1 && static_cast<int>(stepNumber) == exitStep)
  {
    lout() <<"Exit step.  Exiting transient loop\n" << std::endl;
    recoverableFailureFlag = false;
  }

  return recoverableFailureFlag;
}

//-----------------------------------------------------------------------------
// Function      : Transient::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processSuccessfulDCOP()
{
  Stats::StatTop _processSuccessfulDCOPStat("Successful DCOP Steps");
  Stats::TimeBlock _processSuccessfulDCOPTimer(_processSuccessfulDCOPStat);

  bool bsuccess = true;

  loader_.stepSuccess(analysisManager_.getTwoLevelMode());

  if (sensFlag_ && !firstDoubleDCOPStep() )
  {
    nonlinearManager_.calcSensitivity(objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
  }

  stats_.successfulStepsTaken_ += 1;

  // Communicate down to the device level that the step has been accepted.
  // This must be done before the "times" get rotated, and certainly before
  // the vectors get rotated.  The primary target for this call is the
  // transmission line device, which needs to know when it has a valid
  // solution so it can save history.  (And now the ADC/DAC devices too)
  // needs to happen before dcopFlag_ is set to false.
  loader_.acceptStep();

  // Reset some settings (to switch from DCOP to transient, if not the
  // first step of a "double" DCOP.
  if ( firstDoubleDCOPStep() )  // (pde-only)
  {
    dcopFlag_          = true;
    analysisManager_.setTwoLevelMode(Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP);
    baseIntegrationMethod_ = TimeIntg::NoTimeIntegration::type;
  }
  else
  {
    dcopFlag_ = false;
    analysisManager_.setTwoLevelMode(Analysis::TWO_LEVEL_MODE_TRANSIENT);
    analysisManager_.setSwitchIntegrator(true);
    baseIntegrationMethod_   = initialIntegrationMethod_;
    beginningIntegration = true;
  }

  analysisManager_.getDataStore()->setConstantHistory();
  analysisManager_.getDataStore()->computeDividedDifferences();
  analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();

  analysisManager_.getDataStore()->updateSolDataArrays ();

  tranopOutputs(); /* NOTE KIM prn written here */

  // Now that output has been called, update the doubleDCOP step
  // if neccessary. (Only matters for pde problems)
  // doubleDCOPStep_ = lastDCOPStep_;
  nextDCOPStep();

  //Test and save restart if necessary
  if (testRestartSaveTime(analysisManager_, restartManager_, analysisManager_.getStepErrorControl().currentTime, nextRestartSaveTime_))
  {
    if (DEBUG_RESTART)
      dout() << "\n " << analysisManager_.getCommandLine().getArgumentValue("netlist")
             << "  Calling dumpRestart" << std::endl;

    outputManagerAdapter_.dumpRestart(*analysisManager_.getPDSManager()->getPDSComm(),
                                      topology_,
                                      analysisManager_,
                                      restartManager_.getJobName(),
                                      restartManager_.getPack(),
                                      analysisManager_.getStepErrorControl().currentTime);

    if (DEBUG_RESTART)
      dout() << "  Done Calling dumpRestart" << std::endl;
  }

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.
  loader_.outputPlotFiles();

  nonlinearManager_.allocateTranSolver(analysisManager_, analysisManager_.getNonlinearEquationLoader(), linearSystem_, *analysisManager_.getDataStore(), *analysisManager_.getPDSManager(), outputManagerAdapter_.getOutputManager(), topology_);

  analysisManager_.getStepErrorControl().previousCallStepSuccessful = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::doProcessSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::doProcessSuccessfulStep()
{
  Stats::StatTop _processSuccessfulStepStat("Successful Step");
  Stats::TimeBlock _processSuccessfulStepTimer(_processSuccessfulStepStat);

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::TRAN, analysisManager_.getStepErrorControl().nextTime, analysisManager_.getStepErrorControl().getNumberOfSteps()));

  loader_.stepSuccess(analysisManager_.getTwoLevelMode());

  if (timeQueue_.get_size() != 0)
  {
    // store status of this step
    timeQueue_.push_back( analysisManager_.getStepErrorControl().currentTime );
    timeStepQueue_.push_back( analysisManager_.getStepErrorControl().currentTimeStep );
    stepStatusQueue_.push_back( 1 );
    estErrorOverTolQueue_.push_back( analysisManager_.getStepErrorControl().estOverTol_);
    nonlinearSolverStatusQueue_.push_back( analysisManager_.getStepErrorControl().newtonConvergenceStatus );
    nonlinearSolverNumIterationsQueue_.push_back( analysisManager_.getStepErrorControl().nIterations);
    nonlinearSolverMaxNormQueue_.push_back( nonlinearManager_.getNonlinearSolver().getMaxNormF() );
    nonlinearSolverMaxNormIndexQueue_.push_back( nonlinearManager_.getNonlinearSolver().getMaxNormFindex () );
  }

  if (sensFlag_)
  {
    nonlinearManager_.calcSensitivity(objectiveVec_,
                                      dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << "  Transient::processSuccessfulStep()" << std::endl
           << "Newton step succeeded:" << std::endl
           << "nextSolutionPtr: " << std::endl;

    analysisManager_.getDataStore()->nextSolutionPtr->printPetraObject(dout());
    dout() << std::endl;
  }

  // Set things up for the next time step, based on if this one was
  // successful.

  // This output call is for device-specific output (like from a PDE device,
  // outputting mesh-based tecplot files).  It will only work in parallel if on
  // a machine where all processors have I/O capability.
  // Note: this output needs to happen before the "times" get rotated.
  loader_.outputPlotFiles();

  // Communicate down to the device level that the step has been accepted.
  // This must be done before the "times" get rotated, and certainly before
  // the vectors get rotated.  The primary target for this call is the
  // transmission line device, which needs to know when it has a valid
  // solution so it can save history.
  loader_.acceptStep();

  // current time will get updated in completeStep().  We'll save its value
  // for the moment so it can be saved if needed with the rest of the
  // solution if tiaParams_.saveTimeStepsFlag is set.
  // This fixes an off by one bug in getting the right time value and
  // keeps the real solutions associated with that value too.
  double currentTime = analysisManager_.getStepErrorControl().currentTime;
  double suggestedMaxTime = 0.0;
  if (maxTimeStepExpression_)
  {
    suggestedMaxTime = maxTimeStepExpression_->evaluate(comm_, outputManagerAdapter_.getOutputManager().getCircuitTime(),
                                                        analysisManager_.getDataStore()->currSolutionPtr, analysisManager_.getDataStore()->currStatePtr, analysisManager_.getDataStore()->currStorePtr);
  }
  analysisManager_.getStepErrorControl().updateMaxTimeStep( comm_, loader_, tiaParams_, suggestedMaxTime );
  analysisManager_.getStepErrorControl().updateMinTimeStep();
  analysisManager_.getStepErrorControl().updateBreakPoints(loader_, tiaParams_.initialTime);

  if (VERBOSE_TIME && isActive(Diag::TIME_PARAMETERS))
    dout() << "Transient Analysis:  accepting time step" << std::endl;

  analysisManager_.getWorkingIntegrationMethod().completeStep(tiaParams_);

  if (VERBOSE_TIME && tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    dout().precision(15);
    dout() << "ERROROPTION=1: TimeStepLimitedbyBP = " << analysisManager_.getStepErrorControl().TimeStepLimitedbyBP << "\n" << std::endl;
    dout() << "ERROROPTION=1: NL Its =  " << analysisManager_.getStepErrorControl().nIterations << "\n" << std::endl;
    dout() << "ERROROPTION=1: New DeltaT = " << analysisManager_.getStepErrorControl().currentTimeStep << "\n" << std::endl;
  }

  stepNumber     += 1;
  tranStepNumber += 1;
  stats_.successStepsThisParameter_ += 1;
  stats_.successfulStepsTaken_ += 1;

  analysisManager_.getStepErrorControl().numberSuccessiveFailures -= 1;
  if (analysisManager_.getStepErrorControl().numberSuccessiveFailures < 0)
    analysisManager_.getStepErrorControl().numberSuccessiveFailures = 0;

  // --------------------------------------------------------------------
  // Check to see if we have hit the end or a discontinuity point.  The
  // next discontinuity point is represented by the stopTime variable.

  // Note: make sure that when we check for discontinuity point, we take
  // into account the possibility of roundoff error.  Use the same bpTol
  // as is used in function updateBreakPoints.
  double bpTol = analysisManager_.getStepErrorControl().getBreakPointLess().tolerance_;

  if (tiaParams_.bpEnable)
  {
    double timeDiff1 = analysisManager_.getStepErrorControl().currentTime - analysisManager_.getStepErrorControl().stopTime;
    double timeDiff2 = analysisManager_.getStepErrorControl().currentTime - tiaParams_.finalTime;
    timeDiff1 = fabs(timeDiff1);
    timeDiff2 = fabs(timeDiff2);

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << " Checking whether to set breakpointrestartstep" << std::endl;
      dout() << "   current - stop  = " << timeDiff1 << std::endl;
      dout() << "   current - final = " << timeDiff2 << std::endl;
      dout() << "   bpTol           = " << bpTol << std::endl;
      if (timeDiff1 <= bpTol && timeDiff2 > bpTol)
        dout() << "    setting breakPointRestartStep to " << tranStepNumber;
    }

    if (timeDiff1 <= bpTol && timeDiff2 > bpTol)
    {
      analysisManager_.breakPointRestartStep = tranStepNumber;
    }
  }

  if (analysisManager_.breakPointRestartStep == tranStepNumber)
  {
    beginningIntegration = true;
  }
  else
  {
    beginningIntegration = false;
  }

  if (tiaParams_.errorAnalysisOptionResetCount > 1)
  {
    tiaParams_.errorAnalysisOptionResetCount--;
  }
  else if(tiaParams_.errorAnalysisOptionResetCount == 1)
  {
    tiaParams_.errorAnalysisOption = TimeIntg::LOCAL_TRUNCATED_ESTIMATES;
    tiaParams_.errorAnalysisOptionResetCount = 0;
  }

  if (saveTimeStepsFlag)
  {
    analysisManager_.getDataStore()->timeSteps.push_back(currentTime);
    analysisManager_.getDataStore()->timeStepsBreakpointFlag.push_back(beginningIntegration);
    Linear::Vector * aVecPtr = new Linear::Vector( *(analysisManager_.getDataStore()->currSolutionPtr) );
    analysisManager_.getDataStore()->fastTimeSolutionVec.push_back( aVecPtr );
    aVecPtr = new Linear::Vector( *(analysisManager_.getDataStore()->currStatePtr) );
    analysisManager_.getDataStore()->fastTimeStateVec.push_back( aVecPtr );
    aVecPtr = new Linear::Vector( *(analysisManager_.getDataStore()->daeQVectorPtr) );
    analysisManager_.getDataStore()->fastTimeQVec.push_back( aVecPtr );
    aVecPtr = new Linear::Vector( *(analysisManager_.getDataStore()->currStorePtr) );
    analysisManager_.getDataStore()->fastTimeStoreVec.push_back( aVecPtr );
  }

  // 03/16/04 tscoffe:  This is where the solution pointers are rotated.
  analysisManager_.getDataStore()->updateSolDataArrays();

  if (DEBUG_ANALYSIS & DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
    analysisManager_.getDataStore()->outputSolDataArrays(Xyce::dout());

  {
    if  (DEBUG_ANALYSIS)
    {
      Stats::StatTop _outputStat("Output");
      Stats::TimeBlock _outputTimer(_outputStat);
    }

    tranStepOutputs ();
  }

  // Test and save restart if necessary
  if (testRestartSaveTime(analysisManager_, restartManager_, analysisManager_.getStepErrorControl().currentTime, nextRestartSaveTime_))
    outputManagerAdapter_.dumpRestart(*analysisManager_.getPDSManager()->getPDSComm(),
                                      topology_,
                                      analysisManager_,
                                      restartManager_.getJobName(),
                                      restartManager_.getPack(),
                                      analysisManager_.getStepErrorControl().currentTime);

  analysisManager_.getStepErrorControl().previousCallStepSuccessful = true;

  // reset min error tracking variables
  // if the step number is less than zero we'll assume there is no valid
  // min. estimated error over tol or an associated time step.  This frees us
  // from putting Machine::Big here and trying to do something with that.
  stepNumberAtMinEstErrorOverTol = -1;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Transient::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::doProcessFailedStep()
{
  Stats::StatTop _processFailedStat("Failed Steps");
  Stats::TimeBlock _processFailedTimer(_processFailedStat);

  if (analysisManager_.getStepErrorControl().newtonConvergenceStatus <= 0) 
  {
    Stats::StatTop _nonlinearConvergenceFailureStat("Nonlinear Failure");
    Stats::TimeBlock _nonlinearConvergenceFailureTimer(_nonlinearConvergenceFailureStat);
  }

  bool bsuccess = true;

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::TRAN, analysisManager_.getStepErrorControl().nextTime, analysisManager_.getStepErrorControl().getNumberOfSteps()));

  if (timeQueue_.get_size() != 0)
  {
    // store status of this step
    timeQueue_.push_back( analysisManager_.getStepErrorControl().currentTime );
    timeStepQueue_.push_back( analysisManager_.getStepErrorControl().currentTimeStep );
    stepStatusQueue_.push_back( 0 );
    estErrorOverTolQueue_.push_back( analysisManager_.getStepErrorControl().estOverTol_);
    nonlinearSolverStatusQueue_.push_back( analysisManager_.getStepErrorControl().newtonConvergenceStatus );
    nonlinearSolverNumIterationsQueue_.push_back( analysisManager_.getStepErrorControl().nIterations);
    nonlinearSolverMaxNormQueue_.push_back( nonlinearManager_.getNonlinearSolver().getMaxNormF() );
    nonlinearSolverMaxNormIndexQueue_.push_back( nonlinearManager_.getNonlinearSolver().getMaxNormFindex() );
  }

  // save some info about this step
  double estOverTol = analysisManager_.getStepErrorControl().getEstOverTol();
  if( (stepNumberAtMinEstErrorOverTol < 0) || ( estOverTol < minEstErrorOverTol ) )
  {
    // our first failed step, so automatically save info
    stepNumberAtMinEstErrorOverTol = stepNumber;
    minEstErrorOverTol = analysisManager_.getStepErrorControl().getEstOverTol();
    timeStepAtMinEstErrorOverTol = analysisManager_.getStepErrorControl().currentTimeStep;
  }

  loader_.stepFailure(analysisManager_.getTwoLevelMode());

  // DO NOT REMOVE THIS OUTPUT LINE.  It is depended upon by the TIA/ERROROPTION
  // test case, which will fail if this output doesn't happen.
  if (VERBOSE_TIME)
  {
    dout() << "Transient Analysis:  rejecting time step" << std::endl;
  }

  analysisManager_.getWorkingIntegrationMethod().rejectStep(tiaParams_);

  if (VERBOSE_TIME && tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    dout().precision(15);
    dout() << "ERROROPTION=1: TimeStepLimitedbyBP = " 
      << analysisManager_.getStepErrorControl().TimeStepLimitedbyBP << "\n" << std::endl;
    dout() << "ERROROPTION=1: NL Its =  " 
      << analysisManager_.getStepErrorControl().nIterations << "\n" << std::endl;
    dout() << "ERROROPTION=1: New DeltaT = " 
      << analysisManager_.getStepErrorControl().currentTimeStep << "\n" << std::endl;
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << "  Transient::processFailedStep" << std::endl
           << "Newton step failed:" << std::endl
           << "nextSolutionPtr: " << std::endl;
    analysisManager_.getDataStore()->nextSolutionPtr->printPetraObject(dout());
    dout() << std::endl;
  }

  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  if (analysisManager_.getStepErrorControl().currentTimeStep <= 
      analysisManager_.getStepErrorControl().minTimeStep)
  {
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      dout() << "currTimeStep: " << analysisManager_.getStepErrorControl().currentTimeStep << std::endl;
      dout() << "minTimeStep:  " << analysisManager_.getStepErrorControl().minTimeStep << std::endl;
    }

    // before we exit with a time step too small error, check if minTimeStepRecoveryCounter is greater than zero
    // if it is, then the user wants us to try to accept the step that had the minimum
    // estimated error over tol.
    if( tiaParams_.minTimeStepRecoveryCounter > 0 )
    {
      lout() << "Attempting to retake and accept step where estimated error over tolerance was: " 
        << minEstErrorOverTol
        << " and time step was: " << timeStepAtMinEstErrorOverTol << std::endl;

      tiaParams_.minTimeStepRecoveryCounter--;
      bsuccess = retakeAndAcceptTimeStep( timeStepAtMinEstErrorOverTol );
    }
    else
    {
      logQueuedData();
      lout() << "Time step too small near step number: " <<  stepNumber 
        << "  Exiting transient loop.\n" << std::endl;

      bsuccess = false;
    }
  }

  if (tiaParams_.constantTimeStepFlag)
  {
    logQueuedData();
    lout() << "Newton solver failed in constant time step mode.  Exiting transient loop.\n" << std::endl;

    bsuccess = false;
  }

  if (exitStep != -1 && static_cast<int>(stats_.successStepsThisParameter_) == exitStep - 1)
  {
    logQueuedData();
    lout() << "Exit Step.  Exiting transient loop\n" << std::endl;
    bsuccess = false;
  }

#if 1
  if (VERBOSE_TIME && bsuccess)
  {
    outputFailedStepData();
  }
#endif

  if (VERBOSE_TIME && !bsuccess)
  {
    endTRANtime_ = analysisManager_.getXyceTranTimer().elapsedTime();
    tranStats = saveLoopInfo();
    finalVerboseOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processFailedDCOP()
{
  Stats::StatTop _processFailedDCOPStat("Failed DCOP Steps");
  Stats::TimeBlock _processFailedDCOPTimer(_processFailedDCOPStat);

  bool bsuccess = true;

  loader_.stepFailure(analysisManager_.getTwoLevelMode());

  // DC Operating Point failed.
  stats_.failedStepsAttempted_++;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures++;

  lout() << "DC Operating Point Failed.  Exiting transient loop" << std::endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::doFinish()
{
  bool bsuccess = true;

  if (saveTimeStepsFlag && hbAnalysis_)
  {
    TimeIntg::DataStore *dsPtr_ = analysisManager_.getDataStore();
    dsPtr_->timeSteps.push_back(analysisManager_.getStepErrorControl().currentTime);
    dsPtr_->timeStepsBreakpointFlag.push_back(beginningIntegration);
    dsPtr_->fastTimeSolutionVec.push_back(new Linear::Vector( *dsPtr_->currSolutionPtr));
    dsPtr_->fastTimeStateVec.push_back(new Linear::Vector( *dsPtr_->currStatePtr));
    dsPtr_->fastTimeQVec.push_back(new Linear::Vector( *dsPtr_->daeQVectorPtr));
    dsPtr_->fastTimeStoreVec.push_back(new Linear::Vector( *dsPtr_->currStorePtr));
  }

  if (!isPaused)
  {
    if (DEBUG_ANALYSIS)
      dout() << "Calling finishOutput" << std::endl;

    outputManagerAdapter_.finishOutput();

    // This output call is for device-specific output (like from a PDE device,
    // outputting mesh-based tecplot files).  It will only work in parallel if on
    // a machine where all processors have I/O capability, as devices are
    // local to a processor.
    loader_.finishOutput();

    finalVerboseOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::doHandlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool Transient::doHandlePredictor()
{
  /* EDIT KIM */
  TimeIntg::DataStore * ds = analysisManager_.getDataStore();
  Linear::Vector * soln = ds->nextSolutionPtr;
  soln->printPetraObject(std::cout); /* not changed */

  analysisManager_.getDataStore()->setErrorWtVector(tiaParams_);
  /* EDIT KIM */
  soln->printPetraObject(std::cout); /* not changed */
  analysisManager_.getWorkingIntegrationMethod().obtainPredictor();
  /* EDIT KIM */
  soln->printPetraObject(std::cout); /* changed */
  analysisManager_.getWorkingIntegrationMethod().obtainPredictorDeriv();

  /* EDIT KIM */
  soln->printPetraObject(std::cout); /* changed */

  if (DEBUG_ANALYSIS & DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
  {
    dout() << "  Transient::doHandlePredictor" << std::endl;
    analysisManager_.getDataStore()->outputPredictedSolution(Xyce::dout());
    analysisManager_.getDataStore()->outputPredictedDerivative(Xyce::dout());
  }

  // Now, in case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  bool        beginIntegrationFlag = analysisManager_.getBeginningIntegrationFlag();           // system_state.beginIntegrationFlag;
  double      nextTimeStep = analysisManager_.getStepErrorControl().currentTimeStep;           // system_state.nextTimeStep;
  double      nextTime = analysisManager_.getStepErrorControl().nextTime;                      // system_state.nextTime;
  int         currentOrder = analysisManager_.getWorkingIntegrationMethod().getOrder();        // system_state.currentOrder;

  loader_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Transient::resetForHB()
// Purpose       : When doing initial transient run, some analyses require
//                 a reset function
//                 they can fill this in if needed.
// Special Notes :
// Scope         : public
// Creator       : T. Mei, SNL, Parallel Computational Sciences
// Creation Date : 2/26/09
//-----------------------------------------------------------------------------
bool Transient::resetForHB()
{
  dcopFlag_ = false;
  if (timeQueue_.get_size() != 0)
  {
    // these set_size calls will also "reset" the queue.
    timeQueue_.reset_queue();
    timeStepQueue_.reset_queue();
    stepStatusQueue_.reset_queue();
    estErrorOverTolQueue_.reset_queue();
    nonlinearSolverStatusQueue_.reset_queue();
    nonlinearSolverNumIterationsQueue_.reset_queue();
    nonlinearSolverMaxNormQueue_.reset_queue();
    nonlinearSolverMaxNormIndexQueue_.reset_queue();
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Transient::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::finalVerboseOutput()
{
  bool bsuccess = true;

  lout() << "***** Problem read in and set up time: " << analysisManager_.getSolverStartTime() << " seconds" << std::endl;

  if (analysisManager_.getAnalysisMode() == ANP_MODE_TRANSIENT)
  {
    double time;
    if (startTRANtime_ > startDCOPtime)
    {
      time = startTRANtime_ - startDCOPtime;
    }
    else
    {
      time = endTRANtime_ - startDCOPtime;
    }
    lout() << " ***** DCOP time: " << time << " seconds.  Breakdown follows:" << std::endl;

    printLoopInfo(0, dcStats);
  }

  if (analysisManager_.getAnalysisMode() == ANP_MODE_TRANSIENT && endTRANtime_ >= startTRANtime_)
  {
    lout() << " ***** Transient Stepping time: " << endTRANtime_ - startTRANtime_ << " seconds.  Breakdown follows:" << std::endl;

    printLoopInfo(dcStats, tranStats);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::logQueuedData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
void Transient::logQueuedData()
{
  if (timeQueue_.get_size() != 0)
  {
    // get the current non-linear solver return codes
    Nonlinear::ReturnCodes nlReturnCodes = nonlinearManager_.getReturnCodes();

    lout() << " *** Transient failure history: " << std::endl;
    if (tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
    {
      // truncation error is not used here so est err over tol is not useful in the output
      lout() << "Time        Time      Step      Non-Linear Solver      node    node" << std::endl;
      lout() << "(sec)       Step     Status   Status   Iters   ||F||   index   name" << std::endl;
    }
    else
    {
      lout() << "Time        Time      Step   EstErr      Non-Linear Solver      node     node" << std::endl;
      lout() << "(sec)       Step     Status  OverTol   Status  Iters   ||F||    index    name" << std::endl;
    }

    for (int i = 0; i < timeQueue_.get_size(); i++ )
    {
      int fieldWidth=10;
      lout() << std::scientific 
             << std::setprecision(fieldWidth-7) 
             << std::setfill(' ') 
             << std::right 
             << std::setw( fieldWidth )
             << timeQueue_.at_from_tail(i) << "  "
             << timeStepQueue_.at_from_tail(i) << "  ";

      if( stepStatusQueue_.at_from_tail(i) == 1 )
      {
        lout() << "pass  ";
      }
      else
      {
        lout() << "fail  ";
      }
      if (tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
      {
      }
      else
      {
        lout() << estErrorOverTolQueue_.at_from_tail(i) << "  ";
      }
      int nlStatus = nonlinearSolverStatusQueue_.at_from_tail(i);
      lout() << std::setw(7) << std::right;
      if( nlStatus == nlReturnCodes.normTooSmall )
      {
        lout() << "P:s nrm";
      }
      else if( nlStatus == nlReturnCodes.normalConvergence)
      {
        lout() << "pass   ";
      }
      else if( nlStatus == nlReturnCodes.nearConvergence )
      {
        lout() << "P:near ";
      }
      else if( nlStatus == nlReturnCodes.smallUpdate )
      {
        lout() << "P:s up ";
      }
      else if( nlStatus == nlReturnCodes.nanFail )
      {
        lout() << "F:NaN  ";
      }
      else if( nlStatus == nlReturnCodes.tooManySteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.tooManyTranSteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.updateTooBig )
      {
        lout() << "F:big u";
      }
      else if( nlStatus == nlReturnCodes.stalled )
      {
        lout() << "F:stall";
      }
      else if( nlStatus == nlReturnCodes.wrmsExactZero )
      {
        lout() << "F:n zro";
      }
      else if( nlStatus == nlReturnCodes.innerSolveFailed )
      {
        lout() << "F:in Fl";
      }
      else
      {
        lout() << "code=" <<
          nonlinearSolverStatusQueue_.at_from_tail(i) << "  ";
      }

      lout() << std::right << std::setw( 4 )
             << nonlinearSolverNumIterationsQueue_.at_from_tail(i) << "  "
             << nonlinearSolverMaxNormQueue_.at_from_tail(i) ;

      int outIndex = nonlinearSolverMaxNormIndexQueue_.at_from_tail(i) ;
      lout() << std::right << std::fixed << std::setw( 7 ) << outIndex;

      const std::vector<const std::string *> &name_vec = topology_.getSolutionNodeNames();

      lout() << "    " << std::left << ((outIndex < name_vec.size() && outIndex >= 0) ? *name_vec[outIndex] : "N/A") << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Transient::outputFailedStepData()
// Purpose       : outputs same information as logQueuedData, but only for
//                 a single time step.  Also different in that it is called 
//                 for any failed step, not just the final exit.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/9/2015
//-----------------------------------------------------------------------------
void Transient::outputFailedStepData()
{
  if (timeQueue_.get_size() != 0)
  {
    // get the current non-linear solver return codes
    Nonlinear::ReturnCodes nlReturnCodes = nonlinearManager_.getReturnCodes();

    lout() << "************\n*** Transient failure information: " << std::endl;
    if (tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
    {
      // truncation error is not used here so est err over tol is not useful in the output
      lout() << "Time        Time      Step      Non-Linear Solver      node    node" << std::endl;
      lout() << "(sec)       Step     Status   Status   Iters   ||F||   index   name" << std::endl;
    }
    else
    {
      lout() << "Time        Time      Step   EstErr      Non-Linear Solver      node     node" << std::endl;
      lout() << "(sec)       Step     Status  OverTol   Status  Iters   ||F||    index    name" << std::endl;
    }

    {
      int i=timeQueue_.get_size()-1;
      int fieldWidth=10;
      lout() << std::scientific 
             << std::setprecision(fieldWidth-7) 
             << std::setfill(' ') 
             << std::right 
             << std::setw( fieldWidth )
             << timeQueue_.at_from_tail(i) << "  "
             << timeStepQueue_.at_from_tail(i) << "  ";

      if( stepStatusQueue_.at_from_tail(i) == 1 )
      {
        lout() << "pass  ";
      }
      else
      {
        lout() << "fail  ";
      }
      if (tiaParams_.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
      {
      }
      else
      {
        lout() << estErrorOverTolQueue_.at_from_tail(i) << "  ";
      }
      int nlStatus = nonlinearSolverStatusQueue_.at_from_tail(i);
      lout() << std::setw(7) << std::right;
      if( nlStatus == nlReturnCodes.normTooSmall )
      {
        lout() << "P:s nrm";
      }
      else if( nlStatus == nlReturnCodes.normalConvergence)
      {
        lout() << "pass   ";
      }
      else if( nlStatus == nlReturnCodes.nearConvergence )
      {
        lout() << "P:near ";
      }
      else if( nlStatus == nlReturnCodes.smallUpdate )
      {
        lout() << "P:s up ";
      }
      else if( nlStatus == nlReturnCodes.nanFail )
      {
        lout() << "F:NaN  ";
      }
      else if( nlStatus == nlReturnCodes.tooManySteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.tooManyTranSteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.updateTooBig )
      {
        lout() << "F:big u";
      }
      else if( nlStatus == nlReturnCodes.stalled )
      {
        lout() << "F:stall";
      }
      else if( nlStatus == nlReturnCodes.wrmsExactZero )
      {
        lout() << "F:n zro";
      }
      else if( nlStatus == nlReturnCodes.innerSolveFailed )
      {
        lout() << "F:in Fl";
      }
      else
      {
        lout() << "code=" <<
          nonlinearSolverStatusQueue_.at_from_tail(i) << "  ";
      }

      lout() << std::right << std::setw( 4 )
             << nonlinearSolverNumIterationsQueue_.at_from_tail(i) << "  "
             << nonlinearSolverMaxNormQueue_.at_from_tail(i) ;

      int outIndex = nonlinearSolverMaxNormIndexQueue_.at_from_tail(i) ;
      lout() << std::right << std::fixed << std::setw( 7 ) << outIndex;

      const std::vector<const std::string *> &name_vec = topology_.getSolutionNodeNames();

      lout() << "    " << std::left << ((outIndex < name_vec.size() && outIndex >= 0) ? *name_vec[outIndex] : "N/A") << std::endl;
    }
    lout() << "************" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::takeAnIntegrationStep_
// Purpose       : Take a transient integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void print(Epetra_MultiVector *vec, std::string output, int disambiguation) {
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        printf(output.c_str(), i, 0, values[i]);
}
void
Transient::takeAnIntegrationStep_()
{

    /* EDIT KIM */
    TimeIntg::DataStore * ds = analysisManager_.getDataStore();
    Linear::Vector * soln = ds->nextSolutionPtr;
   // soln->printPetraObject(std::cout); /* not changed */
    print(&soln->epetraObj(),"%d,%d,%0.64e\n", 0); 
    std::cout << "done!!!\n";
    

  { // Integration step predictor
    if  (DEBUG_ANALYSIS)
    {
      Stats::StatTop _predictorStat("Predictor");
      Stats::TimeBlock _predictorTimer(_predictorStat);
    }

    doHandlePredictor();
  }
    /* EDIT KIM */
    soln->printPetraObject(std::cout); /* Changed */

  { // Load B/V source devices with time data
 
    if  (DEBUG_ANALYSIS)
    {
      Stats::StatTop _updateDeviceSourceStat("Update Device Sources");
      Stats::TimeBlock _updateDeviceSourceTimer(_updateDeviceSourceStat);
    }

    loader_.updateSources();

    /* EDIT KIM */
    soln->printPetraObject(std::cout);
  }

  { // Nonlinear solve
    Stats::StatTop _nonlinearSolveStat("Nonlinear Solve");
    Stats::TimeBlock _nonlinearSolveTimer(_nonlinearSolveStat);

//    analysisManager_.getDataStore()->dFdxMatrixPtr->printPetraObject(std::cout);
    analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve(); /* soln has changed */
 //   analysisManager_.getDataStore()->dFdxMatrixPtr->printPetraObject(std::cout);
 //   NOTE dFdx does not correspond to the strings here
  }

  {
//    Stats::StatTop _errorStat("Error Estimation");
//    Stats::TimeBlock _errorTimer(_errorStat);

    // Add change to solution
    analysisManager_.getDataStore()->stepLinearCombo();

    gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
    analysisManager_.getStepErrorControl().nIterations = nonlinearManager_.getNonlinearSolver().getNumIterations();

    // Evaluate result


    analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::twoLevelStep
//
// Purpose       : Take a transient integration step, for inner 2-level solve.
//
// Special Notes : Same as takeAnIntegrationStep, but does not do the
//                 prediction. (that is in handlePredictor).
//
//                 The prediction is handled separately, as for a 2-level
//                 solve, you only want to do the prediction on the first
//                 solve of the attempted time step.
//
// Scope         : public
// Creator       :
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::twoLevelStep()
{
  loader_.updateSources();
  analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  analysisManager_.getDataStore()->stepLinearCombo ();

  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);

  endTRANtime_ = analysisManager_.getXyceTranTimer().elapsedTime();
  tranStats = saveLoopInfo();

  return analysisManager_.getStepErrorControl().stepAttemptStatus;
}

//-----------------------------------------------------------------------------
// Function      : Transient::retakeAndAcceptTimeStep
// Purpose       : Do a requested time step and accept it
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
bool Transient::retakeAndAcceptTimeStep( double aTimeStep )
{
  bool bsuccess=true;
  // This function was put in place to handle the following situation.
  // If a series of time steps are rejected because the estimated error
  // over tol. is too high, Xyce may exit if the next time step falls under
  // the minTimeStep.  A user could adjust the tolerances of a simulation
  // to try and get around this, but in UQ studies where there are 1,000s of
  // simulations, it can be difficult to get all of them to run.  So,
  // if the user has set the netlist option timeint MINTIMESTEPRECOVERY=<int>
  // and we're about to exit because the time step has fallen under the minTimeStep
  // we will try and retake the time step that had the min. est. error over tol.
  // and accept that.
  //
  // At this point, Transient::processFailedStep() has already determined
  // that the preconditions outlined above are in place (i.e. the user requested
  // this and we're about to exit with a time step too small error), so lets
  // try the step that had the min. est error over tol.

  // set time step
  analysisManager_.getStepErrorControl().currentTimeStep = timeStepAtMinEstErrorOverTol;

  // take step
  takeAnIntegrationStep_();

  // can't accept step if the non-linear solver failed
  if(analysisManager_.getStepErrorControl().nIterations==0)
  {
    lout() << "Time step too small near step number: " <<  stepNumber << "  Exiting transient loop.\n" << std::endl;
    bsuccess = false;
  }
  else
  {
    doProcessSuccessfulStep();
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void Transient::printStepHeader(std::ostream &os)
{
  if (VERBOSE_TIME) 
  {
    os << "***** " << (DEBUG_ANALYSIS ? analysisManager_.getCommandLine().getArgumentValue("netlist") : "") << "  ";

    if (dcopFlag_)
    {
      os << "Start of DCOP STEP                        # ";
    }
    else
    {
      if (beginningIntegration)
      {
        if (analysisManager_.getStepErrorControl().currentTime == tiaParams_.initialTime)
        {
          os << "Start of Time Step (INITIAL STEP)         # ";
        }
        else
        {
          os << "Start of Time Step (DISCONTINUITY STEP)   # ";
        }
      }
      else
      {
        if (!analysisManager_.getStepErrorControl().stepAttemptStatus)
        {
          os << "Start of Time Step (AFTER FAILED STEP)    # ";
        }
        else
        {
          os <<  "Start of Time Step (AFTER SUCCESS STEP)   # ";
        }
      }
    }

    os << stats_.successStepsThisParameter_ + 1 << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::printProgress()
// Purpose       : Outputs run completion percentage and estimated
//                 time-to-completion.
//
// Special Notes : This will need some fixing to work with .STEP.
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/07/2002
//-----------------------------------------------------------------------------
void Transient::printProgress(std::ostream &os)
{
  if (analysisManager_.getProgressFlag())
  {
    // Percent of the overall requested simulation time which is completed.
    double percentComplete;
    // Average CPU time per time step.
    double aveCPUTimePerStep;
    // Average clock time per time step.
    double aveSimTimePerStep;
    // Estimated CPU time to complete the simulation.
    double estCompletionTime = 0.0;

    // Report the beginning of the DC OP calculation.  First call in OP.
    if (dcopFlag_)
    {
      startDCOPtime = analysisManager_.getXyceTranTimer().elapsedTime();
      if (gui_)
      {
        Report::signalProgress("Xyce in DC Operating Point Calculation");
      }
      os << "***** Beginning DC Operating Point Calculation...\n" << std::endl;
    }
    else if (firstTime && stats_.successStepsThisParameter_ == 1)
    {
      startTRANtime_ = analysisManager_.getXyceTranTimer().elapsedTime();
      dcStats = saveLoopInfo();
      firstTime = false;
      if (gui_)
      {
        Report::signalProgress("Xyce in Transient Calculation");
      }
      os << "***** Beginning Transient Calculation...\n" << std::endl;
    }
    if (analysisManager_.getAnalysisMode() == ANP_MODE_TRANSIENT && stats_.successStepsThisParameter_  > 0)
    {
      if (startSimTime == -1.0)
        startSimTime = analysisManager_.getStepErrorControl().currentTime;

      double diff1 = fabs(analysisManager_.getStepErrorControl().currentTime - tiaParams_.initialTime);
      double diff2 = fabs(tiaParams_.finalTime - tiaParams_.initialTime);

      // 02/05/08 tscoffe:  optimization differences after the TIA_REFACTOR
      // caused a floating point difference here and this resolves that
      // difference.  After the refactor is merged back onto the main trunk, this
      // stuff should be removed.
      percentComplete = 100.0 * diff1/diff2;

      if (fabs(percentComplete - oldPercentComplete) > 1.0)
      {
        oldPercentComplete = percentComplete;

        aveCPUTimePerStep = analysisManager_.getXyceTranTimer().elapsedTime() /
                            (stats_.successfulStepsTaken_ + 1);
        aveSimTimePerStep = (analysisManager_.getStepErrorControl().currentTime - startSimTime) /
                            (stats_.successfulStepsTaken_ + 1);

        if (aveSimTimePerStep > Util::MachineDependentParams::MachineEpsilon())
          estCompletionTime = aveCPUTimePerStep *
                              fabs(tiaParams_.finalTime - analysisManager_.getStepErrorControl().currentTime) /
                              aveSimTimePerStep;

        if (!gui_)
        {
          os << "***** Percent complete: " <<  percentComplete <<  " %" << std::endl;
        }

        if (estCompletionTime > Util::MachineDependentParams::MachineEpsilon())
        {
          unsigned int days, hours, minutes, seconds;
          days    = static_cast<int> (estCompletionTime / 86400);
          hours   = static_cast<int> ((estCompletionTime - days * 86400) / 3600);
          minutes = static_cast<int> ((estCompletionTime - days * 86400 - hours * 3600) / 60);
          seconds = static_cast<int> (estCompletionTime - days * 86400 - hours * 3600 - minutes * 60);

          char timeStr[256];
          for (char *c = timeStr; c != timeStr + sizeof(timeStr); ++c)
            *c = 0;

          if (Parallel::rank(comm_) == 0) 
          {
            // get current local system time
            time_t t = time( NULL );
            struct tm * now = localtime( &t );

            // format and display output
            if ( ( t != (time_t)-1 ) && ( strftime( timeStr, 255, "%c", now ) != 0 ) )
            {
              os << "***** Current system time: " << timeStr << std::endl;
            }

            else
            {
              os << "***** Current system time could not be determined." << std::endl;
            }
          }

          if (days > 0)
            sprintf(timeStr, "%3d days, %2d hrs., %2d min., %2d sec.", days, hours, minutes, seconds);
          else if (hours > 0)
            sprintf(timeStr, "%2d hrs., %2d min., %2d sec.", hours, minutes, seconds);
          else if (minutes > 0)
            sprintf(timeStr, "%2d min., %2d sec.", minutes, seconds);
          else
            sprintf(timeStr, "%2d sec.", seconds);

          if (gui_)
          {
            std::ostringstream ost;
            ost << "Xyce transient ";
            if (percentComplete < 10)
              ost.precision(2);
            else
              ost.precision(3);
            ost << percentComplete
                << "%% complete ... Estimated time to completion: "
                << timeStr << std::endl;
            Report::signalProgress(ost.str());
          }
          else
          {
            os << "***** Estimated time to completion: " << timeStr << std::endl << std::endl;
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::noopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::noopOutputs()
{
  if (testOutputTime(analysisManager_.getStepErrorControl().currentTime, analysisManager_.getNextOutputTime(), tiaParams_.initialOutputTime))
  {
    if ( !firstDoubleDCOPStep() )
    {
      outputManagerAdapter_.tranOutput(
        analysisManager_.getStepErrorControl().currentTime,
        *analysisManager_.getDataStore()->currSolutionPtr,
        *analysisManager_.getDataStore()->currStatePtr,
        *analysisManager_.getDataStore()->currStorePtr,
        *analysisManager_.getDataStore()->currLeadCurrentPtr,
        *analysisManager_.getDataStore()->currLeadDeltaVPtr,
        *analysisManager_.getDataStore()->currLeadCurrentQDerivPtr,
        objectiveVec_,
        dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

      if (outputAdapter_)
        outputAdapter_->outputMPDE(
          analysisManager_.getStepErrorControl().currentTime,
          analysisManager_.getDataStore()->nextSolutionPtr);
    }
    double next_output_time = updateOutputTime(analysisManager_.getStepErrorControl().currentTime,
                                               analysisManager_.getNextOutputTime(),
                                               tiaParams_.finalTime,
                                               outputManagerAdapter_.getInitialOutputInterval(),
                                               outputManagerAdapter_.getOutputIntervals());
    analysisManager_.setNextOutputTime(next_output_time);
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::tranopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::tranopOutputs ()
{
  //Test and output if necessary
  if (testOutputTime(analysisManager_.getStepErrorControl().currentTime, analysisManager_.getNextOutputTime(), tiaParams_.initialOutputTime))
  {
    // Make sure this isn't the NLP step of a PDE DCOP.
    if ( !firstDoubleDCOPStep() )
    {
      if (DEBUG_ANALYSIS)
      {
        dout() << "Transient::tranopOutputs:" << std::endl;
        dout() << "current lead current vector:" << std::endl;

        analysisManager_.getDataStore()->currLeadCurrentPtr->printPetraObject(dout());

        dout() << "current lead current deltaV vector:" << std::endl;
        analysisManager_.getDataStore()->currLeadDeltaVPtr->printPetraObject(dout());
      }

      /* NOTE KIM the dFdxMatrix is not similar to the strings here. This might be because of the dcop estimation. */
      //analysisManager_.getDataStore()->dFdxMatrixPtr->printPetraObject(std::cout); /* Print EDIT KIM */ 
      
      outputManagerAdapter_.tranOutput(analysisManager_.getStepErrorControl().currentTime, /* NOTE KIM output written here */
                                       *analysisManager_.getDataStore()->currSolutionPtr,
                                       *analysisManager_.getDataStore()->currStatePtr,
                                       *analysisManager_.getDataStore()->currStorePtr,
                                       *analysisManager_.getDataStore()->currLeadCurrentPtr,
                                       *analysisManager_.getDataStore()->currLeadDeltaVPtr,
                                       *analysisManager_.getDataStore()->currLeadCurrentQDerivPtr,
                                       objectiveVec_,
                                       dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

      if (outputAdapter_)
      {
        outputAdapter_->outputMPDE(
          analysisManager_.getStepErrorControl().currentTime,
          analysisManager_.getDataStore()->nextSolutionPtr);
      }
    }

    double next_output_time = updateOutputTime(analysisManager_.getStepErrorControl().currentTime,
                                               analysisManager_.getNextOutputTime(),
                                               tiaParams_.finalTime,
                                               outputManagerAdapter_.getInitialOutputInterval(),
                                               outputManagerAdapter_.getOutputIntervals());
    analysisManager_.setNextOutputTime(next_output_time);
  }

  // SAVE and DCOP restart:
  if (initialConditionsManager_.getDCOPRestartFlag()
    || testSaveOutputTime(analysisManager_, initialConditionsManager_))
  {
    initialConditionsManager_.outputDCOP(outputManagerAdapter_.getComm(), topology_.getSolutionNodeNameMap(), *analysisManager_.getDataStore()->currSolutionPtr);
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::tranStepOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::tranStepOutputs()
{
  // The printOutputSolution function will sometimes perform
  // interpolations between the current step and the previous step,
  // if the integrator is using a high order of integration during
  // a new-DAE simulation.

  // If the user has specified tstart, and this happens to be
  // the first output of the simulation (at t=tstart, instead of
  // t=0) the force the interpolator to *not* interpolate, even
  // if it is doing so normally.
  bool doNotInterpolate = beginningIntegration && tiaParams_.initialOutputTimeGiven && firstTranOutput_;

  if (testOutputTime(analysisManager_.getStepErrorControl().currentTime, analysisManager_.getNextOutputTime(), tiaParams_.initialOutputTime))
  {
    if (outputAdapter_)
    {
      outputAdapter_->outputMPDE(
        analysisManager_.getStepErrorControl().currentTime,
        analysisManager_.getDataStore()->nextSolutionPtr);

      if (!outputAdapter_->outputFunkyMPDE())
      {
        // try normal new dae output
        std::vector<double> output_interpolation_times =
          computeOutputInterpolationTimes(analysisManager_.getStepErrorControl().currentTime,
                                          analysisManager_.getNextOutputTime(),
                                          tiaParams_.finalTime,
                                          outputManagerAdapter_.getInitialOutputInterval(),
                                          outputManagerAdapter_.getOutputIntervals());

        analysisManager_.getWorkingIntegrationMethod().printOutputSolution(
          outputManagerAdapter_,
          tiaParams_,
          analysisManager_.getStepErrorControl().currentTime,
          analysisManager_.getDataStore()->currSolutionPtr,
          doNotInterpolate,
          output_interpolation_times,
          false) ;
      }
    }
    else
    {
      std::vector<double> output_interpolation_times =
        computeOutputInterpolationTimes(analysisManager_.getStepErrorControl().currentTime,
                                        analysisManager_.getNextOutputTime(),
                                        tiaParams_.finalTime,
                                        outputManagerAdapter_.getInitialOutputInterval(),
                                        outputManagerAdapter_.getOutputIntervals());

      analysisManager_.getWorkingIntegrationMethod().printOutputSolution(
        outputManagerAdapter_,
        tiaParams_,
        analysisManager_.getStepErrorControl().currentTime,
        analysisManager_.getDataStore()->currSolutionPtr,
        doNotInterpolate,
        output_interpolation_times,
        false) ;
    }

    firstTranOutput_ = false;

    double next_output_time = updateOutputTime(analysisManager_.getStepErrorControl().currentTime,
                                               analysisManager_.getNextOutputTime(),
                                               tiaParams_.finalTime,
                                               outputManagerAdapter_.getInitialOutputInterval(),
                                               outputManagerAdapter_.getOutputIntervals());
    analysisManager_.setNextOutputTime(next_output_time);
  }
  else
  {
    std::vector<double> output_interpolation_times =
      computeOutputInterpolationTimes(analysisManager_.getStepErrorControl().currentTime,
                                      analysisManager_.getNextOutputTime(),
                                      tiaParams_.finalTime,
                                      outputManagerAdapter_.getInitialOutputInterval(),
                                      outputManagerAdapter_.getOutputIntervals());

    analysisManager_.getWorkingIntegrationMethod().printOutputSolution(
      outputManagerAdapter_, tiaParams_, analysisManager_.getStepErrorControl().currentTime,
      analysisManager_.getDataStore()->currSolutionPtr,
      doNotInterpolate,
      output_interpolation_times,
      true) ;
  }

  // .SAVE output
  if (testSaveOutputTime(analysisManager_, initialConditionsManager_)) 
  {
    analysisManager_.getWorkingIntegrationMethod().saveOutputSolution(
      outputManagerAdapter_.getComm(),
      initialConditionsManager_,
      topology_.getSolutionNodeNameMap(),
      tiaParams_,
      analysisManager_.getDataStore()->currSolutionPtr,
      initialConditionsManager_.getSaveTime(),
      doNotInterpolate);
  }
}

//-----------------------------------------------------------------------------
// Function      : testSaveOutputTime
// Purpose       : Similar to testOutputTime, except that this is for
//                 .SAVE files.
// Special Notes : Only outputs 1x.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool
testSaveOutputTime(
  Analysis::AnalysisManager &           analysis_manager,
  IO::InitialConditionsManager &        initial_conditions_manager)
{
  bool flag = true;

  if( !initial_conditions_manager.getSaveFlag())
  {
    flag = false;
  }
  else if (analysis_manager.getStepErrorControl().currentTime < initial_conditions_manager.getSaveTime())
  {
    flag = false;
  }
  else if (analysis_manager.getSavedAlready())
  {
    flag = false;
  }

  if (flag)
  {
    analysis_manager.setSavedAlready(true);
    Xyce::dout() << "Calling SAVE outputs!" <<std::endl;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : testOutputTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool
testOutputTime(
  double                        current_time,
  double                        next_output_time,
  double                        start_time)
{
  bool flag = true;

  if (current_time < start_time)
  {
    flag = false;
  }
  else if ((current_time < next_output_time)&& (fabs(current_time - next_output_time) >= 2 * Util::MachineDependentParams::MachinePrecision()))
  {
    flag = false;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : updateOutputTime
// Purpose       : Advance output time so it is the next one after currTime
// Special Notes : formerly part of testOutputTime
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
double
updateOutputTime(
  double                        current_time,
  double                        next_output_time,
  double                        final_output_time,
  double                        initial_output_interval,
  const IO::IntervalVector &    output_intervals)
{
  // We only need to bother with this if the user has specified output control
  // options
  if (initial_output_interval > 0.0)
  {
    if (output_intervals.empty())
    {
      // Only an initial interval was specified
      while ((next_output_time < current_time)|| (fabs(current_time -   next_output_time) < 2 * Util::MachineDependentParams::MachinePrecision()))
        next_output_time = next_output_time + initial_output_interval;
    }
    else if (current_time < output_intervals[0].first)
    {
      // Multiple intervals, but we're still in the first one
      while (next_output_time <= current_time)
        next_output_time = next_output_time + initial_output_interval;
      if (next_output_time > output_intervals[0].first)
        next_output_time = output_intervals[0].first;
    }
    else
    {
      // Multiple intervals specified, we're past the first one
      std::pair<double, double> currInterval, nextInterval;
      int size = output_intervals.size();
      for (int i = 0; i < size; ++i)
        if (output_intervals[i].first <= current_time)
        {
          currInterval = output_intervals[i];
          if ((i+1) < static_cast<int>(output_intervals.size()))
            nextInterval = output_intervals[i+1];
        }
      int step = static_cast<int> ((current_time-currInterval.first) /
                                   currInterval.second);
      next_output_time = currInterval.first + (step+1)*currInterval.second;
      if (nextInterval.first && (nextInterval.first!=currInterval.first) && (next_output_time>=nextInterval.first))
        next_output_time = nextInterval.first;
    }

    if (next_output_time >= final_output_time)
      next_output_time =final_output_time;
  }

  return next_output_time;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testRestartSaveTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool
testRestartSaveTime(
  AnalysisManager &     analysis_manager,
  IO::RestartMgr &      restart_manager,
  double                current_time,
  double &              next_restart_save_time)
{
  double initial_restart_interval = restart_manager.getInitialRestartInterval();
  const IO::IntervalVector &restart_intervals = restart_manager.getRestartIntervals();

  if (DEBUG_RESTART)
  {
    Xyce::dout() << "TESTING FOR RESTART SAVE" << std::endl
                 << Xyce::subsection_divider << std::endl
                 << "current_time: " << current_time << std::endl
                 << "nextSaveTime: " << next_restart_save_time << std::endl
                 << "initial_restart_interval: " << initial_restart_interval << std::endl;
    if (!restart_intervals.empty())
    {
      Xyce::dout() << "First restart interval: " << restart_intervals[0].first << std::endl;
    }
    else
    {
      Xyce::dout() << "restart_intervals is empty" << std::endl;
    }
  }

  bool flag = false;
  if (initial_restart_interval == 0.0)
  {
    flag = false;
  }
  else if (current_time < next_restart_save_time)
  {
    flag = false;
  }
  else if (restart_intervals.empty())
  {
    while (next_restart_save_time <= current_time)
    {
      next_restart_save_time += initial_restart_interval;
    }
    flag = true;
  }
  else if (current_time < restart_intervals[0].first)
  {
    while (next_restart_save_time <= current_time)
    {
      next_restart_save_time += initial_restart_interval;
    }
    if (next_restart_save_time > restart_intervals[0].first)
    {
      next_restart_save_time = restart_intervals[0].first;
    }
    flag = true;
  }
  else
  {
    std::pair<double, double> currInterval, nextInterval;
    int size = restart_intervals.size();
    for (int i = 0; i < size; ++i)
    {
      if (restart_intervals[i].first <= current_time)
      {
        currInterval = restart_intervals[i];
        if ((i+1) < (int)restart_intervals.size())
        {
          nextInterval = restart_intervals[i+1];
        }
      }
    }
    int step = static_cast <int> ((current_time - currInterval.first) /
                                  currInterval.second);
    next_restart_save_time = currInterval.first + (step+1)*currInterval.second;

    if (nextInterval.first && (nextInterval.first!=currInterval.first)
        && (next_restart_save_time >= nextInterval.first))
    {
      next_restart_save_time = nextInterval.first;
    }
    flag = true;
  }

  if (DEBUG_RESTART)
    Xyce::dout() << "new nextSaveTime: " << next_restart_save_time << std::endl
                 << "restart flag: " << flag << std::endl
                 << Xyce::subsection_divider << std::endl;

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : computeOutputInterpolationTimes
// Purpose       : When we pass "nextOutputTime_", we might have skipped
//                 over points where output was requested.  Make a list of
//                 those times so we can interpolate to them.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling.
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
std::vector<double>
computeOutputInterpolationTimes(
  double                        current_time,
  double                        next_output_time,
  double                        final_output_time,
  double                        initial_output_interval,
  const IO::IntervalVector &    output_intervals)
{
  std::vector<double>          output_interpolation_times;

  // if there are no output control options specified, we do nothing.
  if (initial_output_interval > 0.0)
  {
    if (output_intervals.empty() || current_time <= output_intervals[0].first)
    {
      // we're either in the first interval, or there only *is* one interval
      double t = next_output_time;
      while ((t < current_time) && ( current_time -t ) > 2 * Util::MachineDependentParams::MachinePrecision())
      {
        output_interpolation_times.push_back(t);
        t += initial_output_interval; // outputManagerAdapter_.getInitialOutputInterval();
      }

      if ((t - current_time) <= 2 * Util::MachineDependentParams::MachinePrecision() )
      {
        output_interpolation_times.push_back(current_time);
      }

      if ((t - final_output_time) > 2 * Util::MachineDependentParams::MachinePrecision())
      {
        // tack on finalTime or we'll miss it
        output_interpolation_times.push_back(final_output_time);
      }
    }
    else // we're out of the first interval, and there is more than one
    {
      int outInt,lastInt;

      lastInt=output_intervals.size()-1;

      // find which interval nextOutputTime_ is in
      for (outInt=0;
           outInt<lastInt&&output_intervals[outInt+1].first<=next_output_time;
           ++outInt) ;

      double t = next_output_time;
      while (t <= current_time)
      {
        output_interpolation_times.push_back(t);
        t += output_intervals[outInt].second;
        if (outInt != lastInt && t >= output_intervals[outInt+1].first)
        {
          ++outInt;
          t = output_intervals[outInt].first;
        }
      }
    }
  }

  return output_interpolation_times;
}

//-----------------------------------------------------------------------------
// Function      : writeConductanceFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
void writeConductanceFile(const std::vector<std::string> &device_names, Nonlinear::ConductanceExtractor &conductance_extractor, const std::string &filename)
{
  std::map<std::string,double> inputMap;

  // load inputMap from tiaParam.device_names option
  for (std::vector<std::string>::const_iterator it = device_names.begin(), end = device_names.end(); it != end; ++it) 
    inputMap[*it] = 0.0;

  // if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  // {
  //   Xyce::dout() << "AnalysisManager::conductanceTest()" << std::endl;
  //   for (std::list<std::string>::const_iterator it = tia_params.device_names.begin(), end = tia_params.device_names.end(); it != end; ++it)
  //     Xyce::dout() << "current device name = \"" << *it
  //                  << "\" added to inputMap[ " << *it << " ] = " << inputMap[ *it ] << std::endl;
  // }

  int isize = inputMap.size();
  std::vector<double> outputVector(isize, 0.0);
  std::vector< std::vector<double> > jacobian(isize);
  for (int i = 0; i < isize; ++i)
  {
    jacobian[i].resize(isize, 0.0);
  }

  bool b1 = conductance_extractor.extract(inputMap, outputVector, jacobian);

  int iE1, iE2;
  int numElectrodes = isize;

  FILE *fp1;
  fp1 = fopen(filename.c_str(), "w");

  fprintf(fp1, "%s", "Conductance array: \n");
  fprintf(fp1,"%s", "              ");
  if (b1)
  {
    std::map<std::string,double>::const_iterator iterM = inputMap.begin();
    std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
    for (iE2 = 0; iE2 < numElectrodes; ++iE2, ++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1, "\t%14s", srcname.c_str());
    }
    fprintf(fp1, "%s", "\n");

    iterM = inputMap.begin();
    for (iE1 = 0; iE1 < numElectrodes; ++iE1, ++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1,"%14s",srcname.c_str());
      for (iE2 = 0; iE2 < numElectrodes; ++iE2)
      {
        fprintf(fp1,"\t%14.4e",jacobian[iE1][iE2]);
      }
      fprintf(fp1,"%s", "\n");
    }
    fprintf(fp1,"%s", "\n");
  }
  else
  {
    fprintf(fp1,"%s", "\nConductance calculation failed!\n");
  }

  fclose(fp1);
}

namespace {

typedef Util::Factory<AnalysisBase, Transient>  TransientFactoryBase;

//-----------------------------------------------------------------------------
// Class         : TransientFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing Transient parameters from the netlist and creating Transient analysis.
///
class TransientFactory: public TransientFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : TransientFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the Transient analysis factory
  ///
  /// @invariant Stores the results of parsing, so if more than one of the analysis and
  /// associated lines are parsed, the second options simply overwrite the previously parsed
  /// values.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager
  /// @param linear_system
  /// @param nonlinear_manager
  /// @param topology
  ///
  TransientFactory(
    Analysis::AnalysisManager &         analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager)
    : TransientFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      restartManager_(restart_manager)
  {}

  virtual ~TransientFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:59:00 2015
  //-----------------------------------------------------------------------------
  ///
  /// Create a new Transient analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new Transient analysis object
  ///
  Transient *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_TRANSIENT);
    Transient *transient = new Transient(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_);
    transient->setAnalysisParams(transientAnalysisOptionBlock_);
    transient->setTimeIntegratorOptions(timeIntegratorOptionBlock_);

    return transient;
  }

  //-----------------------------------------------------------------------------
  // Function      : setTransientAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setTransientAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    transientAnalysisOptionBlock_ = option_block;
  }

  //-----------------------------------------------------------------------------
  // Function      : setTimeIntegratorOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:01:27 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the time integrator parsed option block.
  ///
  /// @invariant Overwrites any previously specified time integrator option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setTimeIntegratorOptionBlock(const Util::OptionBlock &option_block)
  {
    timeIntegratorOptionBlock_ = option_block;

    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;

private:
  Util::OptionBlock     transientAnalysisOptionBlock_;
  Util::OptionBlock     timeIntegratorOptionBlock_;
};

// .TRAN
struct TransientAnalysisReg : public IO::PkgOptionsReg
{
  TransientAnalysisReg(
    TransientFactory &        factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setTransientAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  TransientFactory &          factory_;
};


//-----------------------------------------------------------------------------
// Function      : extractTRANData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine_.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool
extractTRANData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("TRAN");

  int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".TRAN line has an unexpected number of fields";
    return false;
  }
  else
  {
    int linePosition = 1;   // Start of parameters on .param line.
    int endPosition = numFields;
    
    Util::Param parameter("", "");
    
    // TSTEP and TSTOP are required, get them now.
    parameter.setTag( "TSTEP" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.
    
    parameter.setTag( "TSTOP" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.
    
    // at this point we can have
    // [ <TSTART> [<DTMAX>] ] [NOOP|UIC] [{expression for time step schedule}]
    //
    // because these are optional, we'll need to test the type of each
    // parameter before deciding how to handle it
    //
    // need these flags to figure out what any untaged doubles are on the .tran line
    bool tstartFound=false;
    bool dtmaxFound=false;
    while( linePosition < endPosition )
    {
      // convert the next item to a parameter and set the
      // parameter's tag and value.  We set both because
      // in some cases we need to change the tag and in
      // others we need to change the value
      parameter.setTag( parsed_line[linePosition].string_ );
      parameter.setVal( parsed_line[linePosition].string_ );
      
      if(parameter.hasExpressionValue())
      {
        // found a max time step schedule expression
        parameter.setTag("MAXTIMEEXPRESSION");
        option_block.addParam( parameter );
      }
      else
      {
        // could be TSTART DTMAX or NOOP|UIC
        if ( parameter.uTag() == "NOOP" || parameter.uTag() == "UIC" )
        {
          parameter.setVal( "1" );
          option_block.addParam( parameter );
        }
        else
        {
          if( !tstartFound )
          {
            parameter.setTag( "TSTART" );
            option_block.addParam( parameter );
            tstartFound=true;
          }
          else if( !dtmaxFound )
          {
            parameter.setTag( "DTMAX" );
            option_block.addParam( parameter );
            dtmaxFound=true;
          }
          else
          {
            Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_)
              << "expected NOOP/UIC field on .TRAN line but found" << parameter.usVal();
          }
        }
      }
      linePosition++;
    }
    
    circuit_block.addOptions(option_block);
    
    return true; // Only get here on success.
  }
}

} // namespace <unnamed>


bool
registerTransientFactory(
  FactoryBlock &        factory_block)
{
  TransientFactory *factory = new TransientFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_, factory_block.restartManager_);

  addAnalysisFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandParser(".TRAN", extractTRANData);
  factory_block.optionsManager_.addCommandParser(".TR", extractTRANData);

  factory_block.optionsManager_.addCommandProcessor("TRAN", new TransientAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions(*factory, &TransientFactory::setTimeIntegratorOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce
