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
// Filename       : $RCSfile: N_LOA_NonlinearEquationLoader.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.38.2.1 $
//
// Revision Date  : $Date: 2015/11/20 00:37:50 $
//
// Current Owner  : $Author: tmei $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Algorithm.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternDevice.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::NonlinearEquationLoader
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
NonlinearEquationLoader::NonlinearEquationLoader(
  TimeIntg::DataStore &                 ds,
  Loader &                              loader,
  Device::DeviceMgr &                   device_manager,
  TimeIntg::WorkingIntegrationMethod &  wim,
  bool                                  daeStateDerivFlag)
  : ds_(ds),
    loader_(loader),
    deviceManager_(device_manager),
    wim_(wim),
    residualTimerPtr_(0),
    jacobianTimerPtr_(0),
    daeStateDerivFlag_(daeStateDerivFlag)
{
  residualTimerPtr_ = new Util::Timer();
  jacobianTimerPtr_ = new Util::Timer();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::~NonlinearEquationLoader
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
NonlinearEquationLoader::~NonlinearEquationLoader()
{
  delete residualTimerPtr_;
  delete jacobianTimerPtr_;
}

bool
NonlinearEquationLoader::analyticSensitivitiesAvailable (std::string & name)
{
  return deviceManager_.analyticSensitivitiesAvailable (name);
}

void
NonlinearEquationLoader::getAnalyticSensitivities(
  std::string &         name,
  std::vector<double> & dfdpVec,
  std::vector<double> & dqdpVec,
  std::vector<double> & dbdpVec,
  std::vector<int> &    FindicesVec,
  std::vector<int> &    QindicesVec,
  std::vector<int> &    BindicesVec) const
{
  return deviceManager_.getAnalyticSensitivities(name, dfdpVec, dqdpVec, dbdpVec, FindicesVec, QindicesVec, BindicesVec);
}

bool NonlinearEquationLoader::setParam(std::string & name, double val, bool overrideOriginal)
{
  return deviceManager_.setParam(name, val, overrideOriginal);
}

bool NonlinearEquationLoader::getParamAndReduce(Parallel::Machine comm, const std::string & name, double & val) const
{
  return Device::getParamAndReduce(comm, deviceManager_, name, val);
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::loadRHS
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the residual (RHS).
//
// Special Notes : All the contributions to the RHS come (in the new DAE
//                 form) from the device package, as Q, F, and B.  The
//                 RHS needs dQdt + F - B.  As dQdt is determined by the
//                 time integration package, the final summation should be
//                 managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::loadRHS ()
{

  /* EDIT KIM */
  std::ofstream ofs;
  ofs.open("load_F.txt");
  ds_.daeFVectorPtr->printPetraObject(ofs);
  ofs.close();

  Linear::Vector fvec(*ds_.daeFVectorPtr);
  

  /* END KIM */

  // Stats::StatTop _residualStat("Residual Load");
  // Stats::TimeBlock _residualTimer(_residualStat);

  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  residualTimerPtr_->resetStartTime();

  ds_.daeQVectorPtr->putScalar(0.0);
  ds_.daeFVectorPtr->putScalar(0.0);
  ds_.daeBVectorPtr->putScalar(0.0);

  ds_.dFdxdVpVectorPtr->putScalar(0.0);
  ds_.dQdxdVpVectorPtr->putScalar(0.0);

  // Update the state. Note - underneath this call, most of the calculations
  // pertaining to the currents, conductances, etc. will happen.
  {
//    Stats::StatTop _deviceResidualStat("Device Load");
//    Stats::TimeBlock _deviceResidualTimer(_deviceResidualStat);

    tmpBool = loader_.updateState(
      ds_.nextSolutionPtr,
      ds_.currSolutionPtr,
      ds_.lastSolutionPtr,
      ds_.nextStatePtr,
      ds_.currStatePtr,
      ds_.lastStatePtr,
      ds_.nextStorePtr,
      ds_.currStorePtr,
      ds_.lastStorePtr);

    bsuccess = bsuccess && tmpBool;

    if (daeStateDerivFlag_)
    {
      wim_.updateStateDeriv ();
    }

    // first load the 2 components: Q, F and B
    tmpBool = loader_.loadDAEVectors(
      ds_.nextSolutionPtr,
      ds_.currSolutionPtr,
      ds_.lastSolutionPtr,
      ds_.nextStatePtr,
      ds_.currStatePtr,
      ds_.lastStatePtr,
      ds_.nextStateDerivPtr,
      ds_.nextStorePtr,
      ds_.currStorePtr,
      ds_.lastStorePtr,
      ds_.nextStoreLeadCurrQPtr,
      ds_.nextLeadCurrentPtr,
      ds_.currLeadCurrentPtr,
      ds_.lastLeadCurrentPtr,
      ds_.nextLeadCurrentQPtr,
      ds_.nextLeadDeltaVPtr,
      ds_.currLeadDeltaVPtr,
      ds_.lastLeadDeltaVPtr,
      ds_.daeQVectorPtr,
      ds_.daeFVectorPtr,
      ds_.daeBVectorPtr,
      ds_.dFdxdVpVectorPtr,
      ds_.dQdxdVpVectorPtr);

    bsuccess = bsuccess && tmpBool;
  }


  wim_.updateLeadCurrent();

  // Now determine dQdt:
  // now sum them all together, to create the total.
  // f(x) is given by:
  //
  //    f(x) = dQ/dt + F(x) - B(t)= 0
  //
  // Note, the nonlinear solver is expecting the RHS vector to
  // contain -f(x).  Or, possibly -f(x) + J*dx, if voltage
  // limiting is on.

  {
//    Stats::StatTop _timeResidualStat("Time Integrator Assembler");
//    Stats::TimeBlock _timeResidualTimer(_timeResidualStat);

    wim_.obtainResidual();
  }

  /* EDIT KIM */
  ofs.open("load_F.txt");
  ds_.daeFVectorPtr->printPetraObject(ofs);
  ofs.close();

  fvec.daxpy(fvec, -1, *ds_.daeFVectorPtr);

  ofs.open("load_diff.txt");
  fvec.printPetraObject(ofs);
  ofs.close();
  /* END KIM */

  // Update the total load time
  residualTime_ = residualTimerPtr_->elapsedTime();


  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::loadJacobian
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the Jacobian.
//
// Special Notes : All the contributions to the Jacobian
//                 come (in the new DAE form)
//                 from the device package, as dQdx and dFdx.
//
//                 The Jacobian is:
//
//                 J = df/dx = d(dQdt)/dx + dF/dx
//
//                 As dQdt is determined by the time integration package,
//                 the final summation should be managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::loadJacobian ()
{
  // Stats::StatTop _jacobianStat("Jacobian Load");
  // Stats::TimeBlock _jacobianTimer(_jacobianStat.getTop());

  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();


  ds_.dQdxMatrixPtr->put(0.0);
  ds_.dFdxMatrixPtr->put(0.0);

  // first load the 2 components: dQdx and dFdx
  {
//    Stats::StatTop _deviceJacobianStat("Device Load");
//    Stats::TimeBlock _deviceJacobianTimer(_deviceJacobianStat);

    tmpBool = loader_.loadDAEMatrices( /* NOTE KIM: the new matrices appear to come from this call */
      ds_.nextSolutionPtr,
      ds_.nextStatePtr,
      ds_.nextStateDerivPtr,
      ds_.nextStorePtr,
      ds_.dQdxMatrixPtr,
      ds_.dFdxMatrixPtr);
    bsuccess = bsuccess && tmpBool;
  }

  // Now determine the d(dQdt)/dx stuff:
  {
//    Stats::StatTop _timeJacobianStat("Time Integrator Assembler");
//    Stats::TimeBlock _timeJacobianTimer(_timeJacobianStat);

    // now sum them all together, to create the total:
    // J = alpha/dt * dQ/dx + dF/dx
    wim_.obtainJacobian();
  }

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();
  
  /* EDIT KIM */
  std::ofstream ofs("currJ.txt");
  ds_.JMatrixPtr->printPetraObject(ofs);
  ofs.close();
  std::cout << "new jacobian\n";
  /* END KIM */

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::applyJacobian
//
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::applyJacobian(
  const Linear::Vector& input,
  Linear::Vector& result)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();

  // first load the 2 components: dQdx and dFdx
  tmpBool = loader_.applyDAEMatrices(
    ds_.nextSolutionPtr,
    ds_.nextStatePtr,
    ds_.nextStateDerivPtr,
    ds_.nextStorePtr,
    input,
    ds_.dQdxVecVectorPtr,
    ds_.dFdxVecVectorPtr);
  bsuccess = bsuccess && tmpBool;

  // Now determine the d(dQdt)/dx stuff:

  // now sum them all together, to create the total:
  // J = alpha/dt * dQ/dx + dF/dx
  wim_.applyJacobian(input, result);

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::loadSensitivityResiduals
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the sensitivity residual (RHS).
//
// Special Notes : All the contributions to the vector come (in the new DAE
//                 form) from the device package, as dQp, dFp, and dBdp.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::loadSensitivityResiduals ()
{
  bool bsuccess = true;
  wim_.obtainSensitivityResiduals ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::loadSensitivityResiduals
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the sensitivity residual (RHS).
//
// Special Notes : All the contributions to the vector come (in the new DAE
//                 form) from the device package, as dQp, dFp, and dBdp.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 09/02/2014
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::loadFinalSensitivityDerivatives ()
{
  bool bsuccess = true;
  wim_.loadFinalSensitivityDerivatives ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::isLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::isLinearSystem() const
{
  return deviceManager_.isLinearSystem();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
int NonlinearEquationLoader::enablePDEContinuation ()
{
  return deviceManager_.enablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::disablePDEContinuation ()
{
  return deviceManager_.disablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::getNumInterfaceNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
void NonlinearEquationLoader::getNumInterfaceNodes
  (std::vector<int> & numINodes)
{
  deviceManager_.getNumInterfaceNodes (numINodes);
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::loadCouplingRHS
  (int iSubProblem, int iCouple, Linear::Vector * dfdvPtr)
{
  return deviceManager_.loadCouplingRHS (iSubProblem, iCouple, dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::calcCouplingTerms
  (int iSubProblem, int iCouple, const Linear::Vector * dxdvPtr)
{
  return deviceManager_.calcCouplingTerms (iSubProblem, iCouple, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::getHomotopyBlockSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
int NonlinearEquationLoader::getHomotopyBlockSize() const
{
  return deviceManager_.getHomotopyBlockSize();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void NonlinearEquationLoader::homotopyStepSuccess(
  const std::vector<std::string> & paramNames,
  const std::vector<double> & paramVals)
{
  const Device::InstanceVector &extern_devices = deviceManager_.getDevices(Device::ExternDevice::Traits::modelType());

  for (Device::InstanceVector::const_iterator it = extern_devices.begin(); it != extern_devices.end(); ++it)
  {
    Device::ExternDevice::Instance &extern_device = static_cast<Device::ExternDevice::Instance &>(*(*it));

    extern_device.homotopyStepSuccess (paramNames, paramVals);
  }
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void NonlinearEquationLoader::homotopyStepFailure ()
{
  const Device::InstanceVector &extern_devices = deviceManager_.getDevices(Device::ExternDevice::Traits::modelType());

  for (Device::InstanceVector::const_iterator it = extern_devices.begin(); it != extern_devices.end(); ++it)
  {
    Device::ExternDevice::Instance &extern_device = static_cast<Device::ExternDevice::Instance &>(*(*it));

    extern_device.homotopyStepFailure ();
  }
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::allDevicesConverged(Parallel::Machine comm)  
{
  return deviceManager_.allDevicesConverged(comm);
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::innerDevicesConverged(Parallel::Machine comm)
{
  return Device::devicesConverged(comm, deviceManager_.getDevices(Device::ExternDevice::Traits::modelType()));
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::getVoltageLimiterStatus
// Purpose       : voltage limiter toggle functions
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/22/2015
//-----------------------------------------------------------------------------
bool NonlinearEquationLoader::getVoltageLimiterStatus()
{
  return deviceManager_.getVoltageLimiterStatus();
}

//-----------------------------------------------------------------------------
// Function      : NonlinearEquationLoader::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/22/2015
//-----------------------------------------------------------------------------
void NonlinearEquationLoader::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  return deviceManager_.setVoltageLimiterStatus(voltageLimterStatus);
}

} // namespace Loader
} // namespace Xyce
