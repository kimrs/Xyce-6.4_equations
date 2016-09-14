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
// Filename       : $RCSfile: N_ANP_OutputMgrAdapter.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.59 $
//
// Revision Date  : $Date: 2015/11/03 17:27:14 $
//
// Current Owner  : $Author: peshola $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_OutputMgrAdapter.h>

#include <N_ANP_NoiseData.h>

#include <N_DEV_DeviceMgr.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_Objective.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputROM.h>
#include <N_IO_RestartMgr.h>
#include <N_DEV_Op.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::OutputMgrAdapter( )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
OutputMgrAdapter::OutputMgrAdapter(
  Parallel::Machine                     comm,
  IO::OutputMgr &                       output_manager,
  IO::Measure::Manager &                measure_manager,
  IO::FourierMgr &                      fourier_manager,
  IO::ObjectiveManager &                objective_manager,
  Device::DeviceMgr &                   device_manager)
  : comm_(comm),
    outputManager_(output_manager),
    measureManager_(measure_manager),
    fourierManager_(fourier_manager),
    objectiveManager_(objective_manager),
    deviceManager_(device_manager),
    tempOp_(new Device::ArtificialParameterOp("TEMP", deviceManager_, *(*deviceManager_.getArtificialParameterMap().find("TEMP")).second, "TEMP")),
    stepSweepVector_(),
    dcSweepVector_(),
    stepAnalysisStepNumber_(0),
    stepAnalysisMaxSteps_(0),
    dcAnalysisStepNumber_(0),
    dcAnalysisMaxSteps_(0)
{}

OutputMgrAdapter::~OutputMgrAdapter()
{
  delete tempOp_;
}

void
OutputMgrAdapter::notify(
  const StepEvent &     event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
    stepAnalysisStepNumber_ = event.count_;
  else if (event.state_ == StepEvent::INITIALIZE)
    stepAnalysisMaxSteps_ = event.count_;
}


void
OutputMgrAdapter::setStepSweepVector(
  const Analysis::SweepVector &       sweep_vector)
{
  stepSweepVector_ = sweep_vector;
}

void
OutputMgrAdapter::setDCSweepVector(
  const Analysis::SweepVector &       sweep_vector)
{
  dcSweepVector_ = sweep_vector;
}

void
OutputMgrAdapter::dumpRestart(
  Parallel::Communicator &      parallel_communicator,
  Topo::Topology &              topology,
  Analysis::AnalysisManager &   analysis_manager,
  const std::string &           job_name,
  bool                          pack,
  double                        current_time) const
{
  IO::dumpRestartData(
    parallel_communicator,
    topology,
    analysis_manager,
    deviceManager_,
    job_name,
    pack,
    current_time);
}


void
OutputMgrAdapter::tranOutput(
  double                time,
  Linear::Vector &      currSolutionPtr,
  Linear::Vector &      stateVecPtr,
  Linear::Vector &      storeVecPtr,
  Linear::Vector &      lead_current_vector,
  Linear::Vector &      junction_voltage_vector,
  Linear::Vector &      lead_current_dqdt_vector,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_,
  bool                  skipPrintLineOutput)
{
//  /* PRINT KIM */ currSolutionPtr.printPetraObject(std::cout);

  fourierManager_.updateFourierData(comm_, time, &currSolutionPtr, &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, &lead_current_dqdt_vector);

  measureManager_.updateTranMeasures(comm_, time, &currSolutionPtr, &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, &lead_current_dqdt_vector);

  // Device::ArtificialParameterOp *op = new Device::ArtificialParameterOp("TEMP", deviceManager_, *(*deviceManager_.getArtificialParameterMap().find("TEMP")).second, "TEMP");
  Util::Op::OpData op_data;
  double temp = (*tempOp_)(comm_, op_data).real();

  outputManager_.output(
    comm_, time, temp, // getParamAndReduce(comm_, deviceManager_, "TEMP"),
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    dcAnalysisStepNumber_, dcAnalysisMaxSteps_, dcSweepVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_,
    skipPrintLineOutput);

  // transient values for the objective function call.
  double arg1 = time;
  double arg2 = 0.0;

  objectiveManager_.outputObjective(comm_, outputManager_.getOpBuilderManager(), outputManager_, outputManager_.getCircuitTime(), arg1, arg2, currSolutionPtr, stateVecPtr, storeVecPtr);
}

void
OutputMgrAdapter::dcOutput(
  int                   dcStepNumber,
  Linear::Vector &      currSolutionPtr,
  Linear::Vector &      stateVecPtr,
  Linear::Vector &      storeVecPtr,
  Linear::Vector &      lead_current_vector,
  Linear::Vector &      junction_voltage_vector,
  Linear::Vector &      lead_current_dqdt_vector,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_)
{
  measureManager_.updateDCMeasures(comm_, dcSweepVector_, &currSolutionPtr, &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, &lead_current_dqdt_vector);

  Util::Op::OpData op_data;
  double temp = (*tempOp_)(comm_, op_data).real();

  outputManager_.output(
    comm_, 0.0, temp, // getParamAndReduce(comm_, deviceManager_, "TEMP"),
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    dcStepNumber, dcAnalysisMaxSteps_, dcSweepVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

  // dc values for the objective function call.
  double arg1 = 0.0;
  double arg2 = 0.0;

  if (dcSweepVector_.size() > 0)
    arg1 = dcSweepVector_[0].currentVal;
  if (dcSweepVector_.size() > 1)
    arg2 = dcSweepVector_[1].currentVal;

  objectiveManager_.outputObjective(comm_, outputManager_.getOpBuilderManager(), outputManager_, outputManager_.getCircuitTime(), arg1, arg2, currSolutionPtr, stateVecPtr, storeVecPtr);
}

void
OutputMgrAdapter::finishOutput()
{
  outputManager_.finishOutput();
}

void
OutputMgrAdapter::outputMPDE(
  double                        time,
  const std::vector<double> &   fast_time_points,
  const Linear::BlockVector &   solution_vector)
{
  outputManager_.outputMPDE(comm_, time, fast_time_points, solution_vector);
}

void
OutputMgrAdapter::outputHB (
  const std::vector<double> &   timePoints,
  const std::vector<double> &   freqPoints,
  const Linear::BlockVector &   timeDomainSolutionVec,
  const Linear::BlockVector &   freqDomainSolutionVecReal,
  const Linear::BlockVector &   freqDomainSolutionVecImaginary,
  const Linear::BlockVector &   timeDomainStoreVec,
  const Linear::BlockVector &   freqDomainStoreVecReal,
  const Linear::BlockVector &   freqDomainStoreVecImaginary,
  const Linear::BlockVector &   timeDomainLeadCurrentVec,
  const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
  const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
  const Linear::BlockVector &   timeDomainJunctionVoltageVec,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary )

{
  outputManager_.outputHB(
    comm_,
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    timePoints, freqPoints,
    timeDomainSolutionVec, freqDomainSolutionVecReal, freqDomainSolutionVecImaginary, 
    timeDomainStoreVec, freqDomainStoreVecReal, freqDomainStoreVecImaginary,
    timeDomainLeadCurrentVec, freqDomainLeadCurrentVecReal, freqDomainLeadCurrentVecImaginary,
    timeDomainJunctionVoltageVec, freqDomainJunctionVoltageVecReal, freqDomainJunctionVoltageVecImaginary );
}

void
OutputMgrAdapter::outputAC (
  double                frequency,
  const Linear::Vector &  solnVecRealPtr,
  const Linear::Vector &  solnVecImaginaryPtr)
{
  measureManager_.updateACMeasures(comm_, frequency, &solnVecRealPtr, &solnVecImaginaryPtr);
  
  outputManager_.outputAC(comm_, frequency, solnVecRealPtr, solnVecImaginaryPtr);

}
void
OutputMgrAdapter::outputNoise (
  double                frequency, 
  double totalOutputNoiseDens_, double totalInputNoiseDens_, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_)
{
  
  outputManager_.outputNoise(comm_, frequency, 
      totalOutputNoiseDens_, totalInputNoiseDens_, noiseDataVec_);
}

// void
// OutputMgrAdapter::outputROM(
//   const Teuchos::SerialDenseMatrix<int, double> &       Ghat,
//   const Teuchos::SerialDenseMatrix<int, double> &       Chat,
//   const Teuchos::SerialDenseMatrix<int, double> &       Bhat,
//   const Teuchos::SerialDenseMatrix<int, double> &       Lhat )
// {
//   IO::outputROM(comm_, outputManager_.getNetlistFilename(), Ghat, Chat, Bhat, Lhat );
// }

// void
// OutputMgrAdapter::outputROM(
//   const Linear::Matrix &                                  Ghat,
//   const Linear::Matrix &                                  Chat,
//   const Teuchos::SerialDenseMatrix<int, double> &       Bhat,
//   const Teuchos::SerialDenseMatrix<int, double> &       Lhat )
// {
//   IO::outputROM(comm_, outputManager_.getNetlistFilename(),  Ghat, Chat, Bhat, Lhat );
// }

double
OutputMgrAdapter::getInitialOutputInterval() const
{
  return outputManager_.getInitialOutputInterval();
}

const IO::IntervalVector &
OutputMgrAdapter::getOutputIntervals() const
{
  return outputManager_.getOutputIntervals();
}

void
OutputMgrAdapter::outputHomotopy(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals,
  Linear::Vector &                        solnVecPtr )
{
  outputManager_.outputHomotopy (comm_, paramNames, paramVals, solnVecPtr);
}

} // namespace Analysis
} // namespace Xyce
