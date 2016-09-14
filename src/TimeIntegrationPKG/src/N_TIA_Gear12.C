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
// Filename      : $RCSfile: N_TIA_Gear12.C,v $
//
// Purpose       : This file contains the functions which define the
//		             backward differentiation, order 1-2, class.
//
// Special Notes :
//
// Creator       : Ting Mei
//
// Creation Date : 2/16/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.64 $
//
// Revision Date  : $Date: 2015/09/03 22:11:40 $
//
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>
/* EDIT KIM */
#include <sstream> 
#include <iostream>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_Gear12.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

using std::abs;

namespace Xyce {
namespace TimeIntg {

const char *
Gear12::name = "Gear 12";

//-----------------------------------------------------------------------------
// Function      : Gear12::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :  10/31/07 
//-----------------------------------------------------------------------------
TimeIntegrationMethod *
Gear12::factory(
    const TIAParams &   tia_params,
    StepErrorControl &  step_error_control,
    DataStore &         data_store)
{
  return new Gear12(tia_params, step_error_control, data_store);
}

//-----------------------------------------------------------------------------
// Function      : Gear::Gear
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------
Gear12::Gear12(
  const TIAParams & tia_params,
  StepErrorControl & secTmp,
  DataStore & dsTmp)
  : TimeIntegrationMethod(),
    leadingCoeff(1.0),
    sec(secTmp),
    ds(dsTmp)
{
  leadingCoeff = 1;
  sec.maxOrder_=(std::min(2,tia_params.maxOrder));
  sec.minOrder_=(std::max(1,tia_params.minOrder));

  if (sec.minOrder_ > sec.maxOrder_)
  {
    sec.minOrder_ = sec.maxOrder_;
  }
  //  sec.maxOrder_ = 2;
  timept_ = -1.0;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainPredictor
// Purpose       : Calculate predictor 
// Special Notes : stored in ds.xn0Ptr,qn0Ptr,qpn0Ptr
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
        elems.push_back(item);
}

template <typename T>
void readBeta(std::vector< T > &out) {
    const std::string FILENAME = "beta.txt";
    std::ifstream ist(FILENAME.c_str());
    if(ist.is_open()) {
        std::string line;
        while ( std::getline(ist, line)) {
            std::cout << line << "\n";
            std::vector<std::string> elems;
            split(line, ',', elems);
            T beta;
            for(int i  = 0; i < elems.size(); ++i)
                beta.push_back(atof(elems[i].c_str()));
            out.push_back(beta);
        }
    }
}

template <typename T>
void writeBeta(T vec) {
    const std::string FILENAME = "beta.txt";
    std::vector< T > betas;
    readBeta(betas);
    betas.push_back(vec);
    
    std::ofstream ofs("beta.txt");
    for(int i = 0; i < betas.size(); ++i) {
        T beta = betas[i];
        std::stringstream ss;
        for(int j = 0; j < beta.size(); j++) {
           ss << std::setprecision(64) << beta[j]; 

           std::cout << "wrote beta: " << std::setprecision(64) << beta[j] << "\n"; 
           if((j+1) < beta.size())
               ss << ",";
        }
        ss << "\n";
        ofs << ss.str();
    }
}

void print(Epetra_MultiVector *vec, char * filename) {
    FILE * file = fopen(filename, "w");
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        fprintf(file, "%d,%d,%0.64e\n", i, 0, values[i]);
    
    fclose(file);
}

void print(Epetra_MultiVector *vec, char * syntax, int disambiguation) {
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        printf(syntax, i, 0, values[i]);
}

void Gear12::obtainPredictor()
{ 
  print(&ds.nextSolutionPtr->epetraObj(), "%d,%d,%0.64e\n", 0);
  // evaluate predictor
  *ds.sn0Ptr = *(ds.sHistory[0]);
  ds.xn0Ptr->putScalar(0.0);
  ds.qn0Ptr->putScalar(0.0);
  *ds.stoQn0Ptr = *(ds.stoLeadCurrQHistory[0]);
  *ds.leadCurrentQn0Ptr = *(ds.leadCurrentQHistory[0]);
  
  ds.spn0Ptr->putScalar(0.0);

  /* EDIT KIM */
  std::cout << "*********** OBTAIN PREDICTOR ***************\n";
  writeBeta(sec.beta_);
  std::cout << "beta: ";
  for(int i = 0; i < sec.beta_.size(); ++i)
      std::cout << sec.beta_[i] << "\t";
  std::cout << "\n";
  std::cout << "order: " << sec.currentOrder_ << "\n";

  std::cout << "******************* HISTORY ******************\n";
  print(&ds.xHistory[0]->epetraObj(), "last_x.txt");
  for(int i = 0; i < 3; ++i) {
      std::cout << "nr: " << i << "\n";
      print(&ds.xHistory[i]->epetraObj(), "%d,%d,%0.64e\n", 0);
  }
  std::cout << "*************** END OF HISTORY ***************\n";


  for (int i=0;i<=sec.currentOrder_;++i)
  {
    std::cout << "******* " << sec.beta_[i] << "\n";
    std::cout << "history\n";
    print(&ds.xHistory[i]->epetraObj(), "%d,%d,%0.64e\n", 0);
    std::cout << "soln\n";
    print(&ds.xn0Ptr->epetraObj(), "%d,%d,%0.64e\n", 0);
    std::cout << "******* BETA: " << sec.beta_[i] << "\t" << "order: " << sec.currentOrder_ << "\n";
    printf("BETA: %0.64e\n", sec.beta_[i]);
    ds.xn0Ptr->linearCombo(sec.beta_[i],*(ds.xHistory[i]),1.0,*ds.xn0Ptr);
    std::cout << "new soln\n";
    print(&ds.xn0Ptr->epetraObj(), "%d,%d,%0.64e\n", 0);

    ds.qn0Ptr->linearCombo(sec.beta_[i],*(ds.qHistory[i]),1.0,*ds.qn0Ptr);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_PREDICTOR))
  {
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainPredictor" << std::endl
      << "\n currentOrder = " << sec.currentOrder_ << std::endl
      << "\n sec.nscsco_: " << sec.nscsco_ << std::endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
      Xyce::dout() << "\n sec.beta_[" << i << "] = " << sec.beta_[i] << "\n" << std::endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << "\n xn0: \n" << std::endl;
    ds.xn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qn0: \n" << std::endl;
    ds.qn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qpn0: \n" << std::endl;
    ds.qpn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n sn0: \n" << std::endl;
    ds.sn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n spn0: \n" << std::endl;
    ds.spn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  // copy the prediction into the next solution:
  /* EDIT KIM */
  print(&ds.xn0Ptr->epetraObj(), "%d,%d,%0.64e\n", 0);

  *(ds.nextSolutionPtr) = *(ds.xn0Ptr);

  obtainSensitivityPredictors();

  return;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainSensitivityPredictors
// Purpose       : Calculate predictor 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::obtainSensitivityPredictors()
{}

  // Computes the step addjustment.
  // 2/16/04 tscoffe:  I'm not exactly sure what this routine is for...
double
Gear12::computeExpoStepAdjust(
  double        stepadjust)
{
  return pow(stepadjust, 1.0 / 3.0);
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainResidual
// Purpose       : Calculate Residual
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
int index = 0; 
void Gear12::obtainResidual()
{
  // output: ds.RHSVectorPtr
  // Note:  ds.nextSolutionPtr is used to get Q,F,B in Analysis::AnalysisManager::loadRHS.
  ds.RHSVectorPtr->linearCombo(sec.alpha_[0],*ds.daeQVectorPtr, sec.alpha_[1],*(ds.qHistory[0]) );
  /* EDIT KIM */
  std::cout << "aplha 1: " << sec.alpha_[0] << "\talpha 2: " << sec.alpha_[1] << "\n";
  std::ofstream ofs;
  ofs.open("gear12_qHistory.txt");
  ds.qHistory[0]->printPetraObject(ofs);
  ofs.close();
  ofs.open("gear12_Q_vec.txt");
  ds.daeQVectorPtr->printPetraObject(ofs);
  ofs.close();
  ofs.open("gear12_1_RHS.txt");
  ds.RHSVectorPtr->printPetraObject(ofs);
  ofs.close();
  /* END KIM */


  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainResidual" << std::endl
      << "\n t = " << sec.nextTime << "\n" << std::endl
      << "\n solution: \n" << std::endl;
    ds.nextSolutionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n daeQVector: \n" << std::endl;
    ds.daeQVectorPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n qn0: \n" << std::endl;
    ds.qn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n qpn0: \n" << std::endl;
    ds.qpn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n sec.alphas_/hn: " << sec.alphas_/sec.currentTimeStep << "\n" << std::endl
      << "\n daeFVector: \n" << std::endl;
    ds.daeFVectorPtr->printPetraObject(Xyce::dout());

    Xyce::dout() << "\n dQdt-vector: \n" << std::endl;
    ds.RHSVectorPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }

  if (sec.currentOrder_  == 2)
  {
    ds.RHSVectorPtr->linearCombo(1.0, *ds.RHSVectorPtr, sec.alpha_[2],*(ds.qHistory[1]));
  }

  if(index != 1) {
    ds.RHSVectorPtr->linearCombo(1.0/sec.currentTimeStep,*ds.RHSVectorPtr,+1.0,*ds.daeFVectorPtr);
    ++index;
  } else {

    ds.RHSVectorPtr->scale(1.0/sec.currentTimeStep);
    ofs.open("gear12_2_RHS.txt");
    ds.RHSVectorPtr->printPetraObject(ofs);
    ofs.close();
    ofs.open("gear12_F.txt");
    ds.daeFVectorPtr->printPetraObject(ofs);
    ofs.close();
  /* END KIM*/
    ds.RHSVectorPtr->linearCombo(+1.0,*ds.RHSVectorPtr,+1.0,*ds.daeFVectorPtr);
  }

  /* EDIT KIM*/
  std::cout << "1.0 / current timestep: " << 1.0/sec.currentTimeStep << "\n";
  ofs.open("gear12_F.txt");
  ds.daeFVectorPtr->printPetraObject(ofs);
  ofs.close();

  ofs.open("gear12_3_RHS.txt");
  ds.RHSVectorPtr->printPetraObject(ofs);
  ofs.close();
  /* END KIM*/

  ds.RHSVectorPtr->linearCombo(1.0,*ds.RHSVectorPtr,-1.0,*ds.daeBVectorPtr);
  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  /* EDIT KIM*/
  ofs.open("gear12_4_RHS.txt");
  ds.RHSVectorPtr->printPetraObject(ofs);
  ofs.close();
  /* END KIM*/

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.dQdxdVpVectorPtr)->scale( sec.alpha_[0]/sec.currentTimeStep );
    //        double qscalar(sec.alpha_[0]/sec.currentTimeStep);

    (ds.RHSVectorPtr)->daxpy(
        *(ds.RHSVectorPtr), +1.0, *(ds.dQdxdVpVectorPtr));

    (ds.RHSVectorPtr)->daxpy(
        *(ds.RHSVectorPtr), +1.0, *(ds.dFdxdVpVectorPtr));
    /*EDIT KIM*/
    ofs.open("dQdxdVp.txt");
    ds.dQdxdVpVectorPtr->printPetraObject(ofs);
    ofs.close();

    ofs.open("dFdxdVp.txt");
    ds.dFdxdVpVectorPtr->printPetraObject(ofs);
    ofs.close();
    /*END KIM*/
    
  }

  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << "\n Residual-vector: \n" << std::endl
      << "-(qpn0-(sec.alpha_s/h)*(Q-qn0)+F-B) \n" << std::endl;
    ds.RHSVectorPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << Xyce::section_divider << std::endl
      << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainSensitivityResiduals 
// Purpose       : Calculate sensitivity residual
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::obtainSensitivityResiduals()
{
  int numParams = ds.sensRHSPtrVector.size();
  for (int ip=0; ip<numParams;++ip)
  {
    Linear::Vector & RHSVec  = *(ds.sensRHSPtrVector[ip]);
    Linear::Vector & dfdpVec = *(ds.nextDfdpPtrVector[ip]);
    Linear::Vector & dqdpVec = *(ds.nextDqdpPtrVector[ip]);
    Linear::Vector & dbdpVec = *(ds.nextDbdpPtrVector[ip]);

    Linear::Vector & currDXdpVec = *(ds.currDXdpPtrVector[ip]);
    Linear::Vector & lastDXdpVec = *(ds.lastDXdpPtrVector[ip]);

    std::vector<Linear::Vector*> & dqdpHistoryVec = ds.dqdpHistory[ip];

    Linear::Vector & currDQdxDXdpVec = *(ds.currDQdxDXdpPtrVector[ip]);
    Linear::Vector & lastDQdxDXdpVec = *(ds.lastDQdxDXdpPtrVector[ip]);

    RHSVec.linearCombo(sec.alpha_[0],dqdpVec, sec.alpha_[1],*(dqdpHistoryVec[0]) );

    if (sec.currentOrder_  == 2)
    {
      RHSVec.linearCombo(1.0, RHSVec, sec.alpha_[2],*(dqdpHistoryVec[1]));
    }

    RHSVec.linearCombo(1.0/sec.currentTimeStep,RHSVec,+1.0,dfdpVec);
    RHSVec.linearCombo(1.0,RHSVec,-1.0,dbdpVec);

    // since the nonlinear solver is expecting a -f, scale by -1.0:
    RHSVec.scale(-1.0);

    // deal with the matvec corrections
    // This assumes that the matvecs have been properly stored.
    double qscalar1(sec.alpha_[1]/sec.currentTimeStep);
    RHSVec.linearCombo(1.0,RHSVec, -qscalar1, currDQdxDXdpVec);

    if (sec.currentOrder_  == 2)
    {
      double qscalar2(sec.alpha_[2]/sec.currentTimeStep);
      RHSVec.linearCombo(1.0,RHSVec, -qscalar2, lastDQdxDXdpVec);
    }

#ifdef DEBUG_SENS
    Xyce::dout() << "obtainSensitivityResiduals: RHS Vector, ip = " << ip << ":\n";
    RHSVec.printPetraObject(Xyce::dout());
#endif
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainJacobian
// Purpose       : Calculate Jacobian
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void Gear12::obtainJacobian()
{

  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainJacobian" << std::endl;
  }

  // output: ds.JMatrixPtr

  // This function returns the following matrix:
  // $-(sec.alphas_/hn)dQdx(x)+dFdx$

  // Note:  ds.nextSolutionPtr is used to get dQdx,dFdx in Analysis::AnalysisManager::loadJacobian.

  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);
  Linear::Matrix & dFdx = *(ds.dFdxMatrixPtr);
  Linear::Matrix & Jac = *(ds.JMatrixPtr);

  double qscalar(sec.alpha_[0]/sec.currentTimeStep);
  double fscalar(1.0);

  Jac.linearCombo( qscalar, dQdx, fscalar, dFdx );
  /* EDIT KIM */
  std::cout << "qscalar: " << qscalar << "\tfscalar: " << fscalar << "\n";
  std::ofstream ofs;

  ofs.open("gear12_dQdx.txt");
  dQdx.printPetraObject(ofs);
  ofs.close();

  ofs.open("gear12_dFdx.txt");
  dFdx.printPetraObject(ofs);
  ofs.close();

  ofs.open("Jacobian.txt");
  Jac.printPetraObject(ofs);
  ofs.close();
  /* END KIM */

  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << "\n dFdx:" <<std::endl;
    dFdx.printPetraObject(Xyce::dout());
    Xyce::dout() << "\n Total Jacobian:" <<std::endl;
    Jac.printPetraObject(Xyce::dout());
    //    for (int i=0;i<3;++i)
    //    {
    //      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
    //    }

    Xyce::dout() << Xyce::section_divider << std::endl << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::interpolateSolution
// Purpose       : Interpolate solution approximation at prescribed time point.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
bool Gear12::interpolateSolution(double timepoint, 
    Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec)

{
  // this is a very course approximation to determine if we are too 
  // close to the actual time step to do an interpolation.
  // it could use more work. 
  double dtr = timepoint - sec.currentTime;  // the delta time requested.
  if( -dtr < 100 * Util::MachineDependentParams::MachinePrecision() )
  {
    *tmpSolVectorPtr = *(historyVec[0]);
    return false;
  }

  tmpSolVectorPtr->linearCombo(1.0, *(historyVec[0]), -1.0, *(historyVec[1]));

  if( sec.usedOrder_ <= 2)
  {
    // do first order interpolation
    // X_interp = X + delta_t_requested * delta_X/delta_t[last step]
    dtr = dtr / sec.lastTimeStep;
    tmpSolVectorPtr->linearCombo(1.0, *(historyVec[0]), dtr, *tmpSolVectorPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::interpolateMPDESolution
// Purpose       : Interpolate solution approximation at prescribed time points.
// Special Notes : This routine computes the solution at the output 
//               : timepoints by intepolation of the history using the order
//               : used for the most recent completed step, orderUsed.
//               : The output is put into provided Linear::Vector pointer.
//               : The interpolation is as follows:
//               : tmpSolVectorPtr->block(i) is interpolated at timepoint(i)
//               : Therefore, if you want them all interpolated at the same time, 
//               : then use timepoint(i) = timepoint(0) forall i
//               : or use interpolateSolution. 
// Scope         : public
// Creator       : Ting Mei, Eric Keiter, SNL 
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool Gear12::interpolateMPDESolution(std::vector<double>& timepoint, 
    Linear::Vector * tmpSolVectorPtr)
{
  Linear::BlockVector * blockTempSolVectorPtr = 
    dynamic_cast<Linear::BlockVector*>(tmpSolVectorPtr);
  if (blockTempSolVectorPtr == NULL)
  {
    std::string msg = "Gear12::interpolateMPDESolution: ";
    msg += "Linear::Vector tmpSolVectorPtr is not of type Linear::BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }

  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t_{n-1}
  int numblocks = timepoint.size();
  int blockCount = blockTempSolVectorPtr->blockCount();
  if (numblocks > blockCount)
  {
    std::string msg = "Gear12::interpolateMPDESolution: ";
    msg += "Number of time points requested is greater than number of fast time points in MPDE block vector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  double delt;  
  double c = 1.0;
  double gam;  
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = sec.usedStep_;
  int kused = sec.usedOrder_;
  double uround = 0.0;  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  for (int i=0; i<numblocks ; ++i)
  {
    if ( (timepoint[i] - tp)*hh < 0.0 ) 
      return false;
  }

  *tmpSolVectorPtr = *(ds.xHistory[0]);

  Linear::Vector * solVectorPtr;
  Linear::Vector * xHistoryVectorPtr;
  // Loop over blocks first so that maximal order can be maintained
  for (int i=0; i < numblocks ; ++i)
  {
    if ((kused == 0) || (timepoint[i] == tn)) { kord = 1; }
    else { kord = kused; }
    solVectorPtr = &(blockTempSolVectorPtr->block(i));
    c = 1.0;
    delt = timepoint[i] - tn;
    gam = delt/sec.psi_[0];
    for (int j=1 ; j <= kord ; ++j)
    {
      c = c*gam;
      gam = (delt + sec.psi_[j-1])/sec.psi_[j];
      Linear::BlockVector * blockXHistoryVectorPtr = 
        dynamic_cast<Linear::BlockVector*>(ds.xHistory[j]);
      if (blockXHistoryVectorPtr == NULL)
      {
        Xyce::Report::DevelFatal0().in("Gear12::interpolateMPDESolution") << "Linear::Vector ds.xHistory[j] is not of type Linear::BlockVector\n j = " << j;
        return(false);
      }
      xHistoryVectorPtr = &(blockXHistoryVectorPtr->block(i));
      solVectorPtr->linearCombo(1.0,*solVectorPtr,c,*xHistoryVectorPtr);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::printMPDEOutputSolution()
// Purpose       : Print transient output from MPDE simulation
// Special Notes : This routine uses interpolateMPDESolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool Gear12::printMPDEOutputSolution(
    Analysis::OutputMgrAdapter & outputManagerAdapter,
    const double time,
    Linear::Vector * solnVecPtr,
    const std::vector<double> & fastTimes )
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Gear12::printWaMPDEOutputSolution()
// Purpose       : Print transient output from WaMPDE simulation
// Special Notes : This routine uses interpolateSolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool Gear12::printWaMPDEOutputSolution(
    Analysis::OutputMgrAdapter & outputManagerAdapter,
    const double time,
    Linear::Vector * solnVecPtr,
    const std::vector<double> & fastTimes,
    const int phiGID )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::printOutputSolution()
// Purpose       : Print output that is dumbed down in order.
// Special Notes : This routine picks smaller time steps to approximate first
//               : order integration from the perspective of the output.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/16/07 
//-----------------------------------------------------------------------------

bool Gear12::printOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter,
  const TIAParams &             tia_params, 
  const double                  time,
  Linear::Vector *                solnVecPtr,
  const bool                    doNotInterpolate,
  const std::vector<double> &   outputInterpolationTimes,
  bool                          skipPrintLineOutput)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Gear12::printOutputSolution" << std::endl
                 << "usedOrder_ = " << sec.usedOrder_ << std::endl;
  }

  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  bool dointerp = true;
  double hh = timestep/(sec.usedOrder_);

  if (hh <= 10*sec.minTimeStep)
  {
    dointerp = false;
  }

  if (!tia_params.interpOutputFlag)
  {
    dointerp = false;
  }

  if (doNotInterpolate)
  {
    dointerp = false;
  }

  if (dointerp && !outputInterpolationTimes.empty())
  {
    for (unsigned int i=0;i<outputInterpolationTimes.size();++i)
    {
      interpolateSolution(outputInterpolationTimes[i], ds.tmpSolVectorPtr, ds.xHistory);    // interpolate solution vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStaVectorPtr, ds.sHistory);    // interpolate state vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStoVectorPtr, ds.stoHistory);  // interpolate store vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentVectorPtr, ds.leadCurrentHistory);  // interpolate store vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadDeltaVPtr, ds.leadDeltaVHistory);  // interpolate store vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentQDerivVectorPtr, ds.leadCurrentQDerivHistory);  // interpolate store vector
      outputManagerAdapter.tranOutput(outputInterpolationTimes[i], *ds.tmpSolVectorPtr, 
                                      *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
                                      ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_, 
                                      ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_,
                                      skipPrintLineOutput);

    }
  }

  // Either way, do an output on the actual computed time step, but only
  // if we weren't given a list of specific times *or* we were told not to
  // interpoloate.
  if (outputInterpolationTimes.empty() || doNotInterpolate)
  {
    outputManagerAdapter.tranOutput(time, *ds.currSolutionPtr, 
                                    *ds.currStatePtr, *ds.currStorePtr, *ds.currLeadCurrentPtr, *ds.currLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
                                    ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_, 
                                    ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_,
                                    skipPrintLineOutput);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    Xyce::dout() << Xyce::section_divider << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::saveOutputSolution
// Purpose       : This is similar to printOutputSolution, but is in support of
//                 the .SAVE capability, rather than .PRINT.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
bool
Gear12::saveOutputSolution( 
  Parallel::Machine                     comm,
  IO::InitialConditionsManager &        initial_conditions_manager,
  const NodeNameMap &                   node_name_map,
  const TIAParams &                     tia_params,
  Linear::Vector *                      solnVecPtr,
  const double                          saveTime,
  const bool                            doNotInterpolate)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::saveOutputSolution" << std::endl;
  }

  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  bool dointerp = true;
  double hh = timestep/(sec.usedOrder_);

  // outputManagerAdapter.outputDCOP( *(ds.currSolutionPtr) );
  initial_conditions_manager.outputDCOP(comm, node_name_map, *ds.currSolutionPtr);

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    Xyce::dout() << Xyce::section_divider << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateHistory
// Purpose       : Update history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void Gear12::updateHistory()
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Gear12::updateHistory" << std::endl
                 << "\n Before updates \n" << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  if (sec.currentOrder_ == 2)
  {
    *(ds.xHistory[2]) = *(ds.xHistory[1]);
  }
  *(ds.qHistory[1]) = *(ds.qHistory[0]);
  *(ds.xHistory[1]) = *(ds.xHistory[0]);

  *(ds.stoLeadCurrQHistory[1]) = *(ds.stoLeadCurrQHistory[0]);
  *(ds.sHistory[1]) = *(ds.sHistory[0]);
  *(ds.stoHistory[1]) = *(ds.stoHistory[0]);

  *(ds.leadCurrentHistory[1]) = *(ds.leadCurrentHistory[0]);
  *(ds.leadCurrentQHistory[1]) = *(ds.leadCurrentQHistory[0]);
  *(ds.leadDeltaVHistory[1]) = *(ds.leadDeltaVHistory[0]);

  *(ds.xHistory[0]) = *ds.nextSolutionPtr;
  /* EDIT KIM */ ds.nextSolutionPtr->printPetraObject(std::cout);
  *(ds.qHistory[0]) =  *ds.daeQVectorPtr;
  *(ds.sHistory[0]) =  *ds.nextStatePtr;    
  *(ds.stoLeadCurrQHistory[0]) = *ds.nextStoreLeadCurrQPtr;
  *(ds.stoHistory[0]) = *ds.nextStorePtr;

  *(ds.leadCurrentHistory[0]) = *ds.nextLeadCurrentPtr;
  *(ds.leadCurrentQHistory[0]) = *ds.nextLeadCurrentQPtr;
  *(ds.leadDeltaVHistory[0]) = *ds.nextLeadDeltaVPtr;

  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << "\n After updates \n" << std::endl;
    Xyce::dout() << "\n newtonCorrectionPtr: " << std::endl;
    ds.newtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n qnewtonCorrectionPtr: " << std::endl;
    ds.qNewtonCorrectionPtr->printPetraObject(Xyce::dout());
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << "\n sNewtonCorrectionPtr: " << std::endl;
    ds.sNewtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n nextStatePtr: " << std::endl;
    ds.nextStatePtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  updateSensitivityHistory();
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateSensitivityHistory
// Purpose       : Update sensitivity history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::updateSensitivityHistory()
{
  int numParams = ds.sensRHSPtrVector.size();
  for (int ip=0; ip<numParams;++ip)
  {
    std::vector<Linear::Vector*> & dqdpHistory = ds.dqdpHistory[ip];
    std::vector<Linear::Vector*> & dXdpHistory = ds.dXdpHistory[ip];
    std::vector<Linear::Vector*> & dQdxdXdpHistory = ds.dQdxdXdpHistory[ip];

    Linear::Vector & nextDQdxDXdpVec = *(ds.nextDQdxDXdpPtrVector[ip]);
    Linear::Vector & nextDXdpVec = *(ds.nextDXdpPtrVector[ip]);
    Linear::Vector & nextDqdpVec = *(ds.nextDqdpPtrVector[ip]);

    if (sec.currentOrder_ == 2)
    {
      *(dQdxdXdpHistory[2]) = *(dQdxdXdpHistory[1]);
      *(dXdpHistory[2]) = *(dXdpHistory[1]);
      *(dqdpHistory[1]) = *(dqdpHistory[0]);
    }

    *(dQdxdXdpHistory[1]) = *(dQdxdXdpHistory[0]);
    *(dXdpHistory[1]) = *(dXdpHistory[0]);

    *(dQdxdXdpHistory[0]) = nextDQdxDXdpVec;
    *(dXdpHistory[0]) = nextDXdpVec;
    *(dqdpHistory[0]) = nextDqdpVec;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::restoreHistory()
{
  for (int i=1;i<=sec.currentOrder_;++i)
  {
    sec.psi_[i-1] = sec.psi_[i];
  }
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::restoreHistory" << std::endl;
    for (int i=1;i<=sec.currentOrder_;++i)
      Xyce::dout() << "\n sec.psi_[i] = " << sec.psi_[i] << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->printPetraObject(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
} 

//-----------------------------------------------------------------------------
// Function      : Gear12::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void Gear12::updateCoeffs()
{
  // synchronize with Step Error Control
  //  sec.psi_[0] = sec.currentTimeStep;
  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::updateCoeffs" << std::endl
      << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
      << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl
      << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl
      << "  nscsco_ = " <<  sec.nscsco_ << std::endl
      << "  psi_[0] = " <<  sec.psi_[0] << std::endl;
  }

  double temp1 = sec.currentTimeStep;

  if (sec.currentOrder_ == 2)
  {
    sec.psi_[2] = sec.psi_[1];
  }
  sec.psi_[1] = sec.psi_[0];
  sec.psi_[0] = temp1;

  //    sec.beta_[0] = 1.0;
  //    sec.alpha_[0] = 1.0; 
  sec.ck_ = 1.0;
  sec.alphas_ = -1.0;

  if (sec.currentOrder_ == 2)
  {
    // the coeffs of predictor
    sec.beta_[2] = temp1/sec.psi_[2] * (temp1 + sec.psi_[1])/(sec.psi_[1] + sec.psi_[2]);
    sec.beta_[1] = -temp1/sec.psi_[1] - sec.beta_[2] * (sec.psi_[1] + sec.psi_[2])/sec.psi_[1];
    sec.beta_[0] = 1.0 - sec.beta_[2] - sec.beta_[1];

    sec.alpha_[2] = -temp1/sec.psi_[1] * temp1/(2 * temp1 + sec.psi_[1]); 
    sec.alpha_[1] = 1 - sec.alpha_[2];
    sec.alpha_[0] = -sec.alpha_[1] - sec.alpha_[2] * (1 + sec.psi_[1]/temp1);

    sec.alpha_[2] = sec.alpha_[2]/sec.alpha_[0];
    sec.alpha_[1] = sec.alpha_[1]/sec.alpha_[0];
    sec.alpha_[0] = -1/sec.alpha_[0];

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1] + sec.psi_[2]);
  }
  else
  {
    sec.beta_[0] = 1.0 + temp1/sec.psi_[1];
    sec.beta_[1] = -temp1/sec.psi_[1];
    sec.alpha_[0] = 1.0;
    sec.alpha_[1] = -1.0;

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1]);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << "  nscsco_ = " <<  sec.nscsco_ << std::endl
      << "  beta_[0] = " <<  sec.beta_[0] << std::endl
      << "  beta_[1] = " <<  sec.beta_[1] << std::endl
      << "  beta_[2] = " <<  sec.beta_[2] << std::endl
      << "  beta_[3] = " <<  sec.beta_[3] << std::endl
      << "  beta_[4] = " <<  sec.beta_[4] << std::endl
      << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl
      << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl
      << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl
      << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl
      << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl
      << "  alphas_ = " <<  sec.alphas_ << std::endl
      << "  alpha0_ = " <<  sec.alpha0_ << std::endl
      << "  gamma_[0] = " <<  sec.gamma_[0] << std::endl
      << "  gamma_[1] = " <<  sec.gamma_[1] << std::endl
      << "  gamma_[2] = " <<  sec.gamma_[2] << std::endl
      << "  gamma_[3] = " <<  sec.gamma_[3] << std::endl
      << "  gamma_[4] = " <<  sec.gamma_[4] << std::endl
      << "  psi_[0] = " <<  sec.psi_[0] << std::endl
      << "  psi_[1] = " <<  sec.psi_[1] << std::endl
      << "  psi_[2] = " <<  sec.psi_[2] << std::endl
      << "  psi_[3] = " <<  sec.psi_[3] << std::endl
      << "  psi_[4] = " <<  sec.psi_[4] << std::endl
      << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl
      << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl
      << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl
      << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl
      << "  sigma_[4] = " <<  sec.sigma_[4] << std::endl
      << "  ck_ = " <<  sec.ck_ << std::endl
      << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void Gear12::initialize(const TIAParams &tia_params)
{
  // we assume the solution vector is available here 
  // Note that I'm using currSolutionPtr instead of
  // nextSolutionPtr because this is the first step.

  // Update next stop time from StepErrorControl:
  // ERK.  Commenting this out, as it is already called from Analysis::AnalysisManager,
  // right before this initialize call.  It should not be called 2x, as 
  // it is history dependent (unfortunately), so calling it 2x in a row changes 
  // the stop time to a different number.
  // sec.updateStopTime();

  // Choose initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;
  double currentTimeStep;

  sec.TimeStepLimitedbyBP =  false; 

  if (tia_params.constantTimeStepFlag)
  {
    currentTimeStep = 0.1 * time_to_stop;
    currentTimeStep = std::min(sec.startingTimeStep, currentTimeStep);
    sec.currentTimeStep = currentTimeStep;
  }
  else
  {
    // compute an initial step-size based on rate of change in the 
    // solution initially
    double dnorm_q = ds.delta_x_errorNorm_q1();
    if (dnorm_q > 0.0)  // time-dependent DAE
    {
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = std::min(sec.h0_max_factor_*abs(time_to_stop),sqrt(2.0)/(sec.h0_safety_*dnorm_q));
      else
        currentTimeStep = 0.1* std::min(sec.savedTimeStep, abs(time_to_stop));
    } 
    else  // non-time-dependent DAE
    {
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = sec.h0_max_factor_*abs(time_to_stop);
      else
        currentTimeStep = 0.1* std::min(sec.savedTimeStep, abs(time_to_stop));
    }

    // choose min of user specified value and our value:
    if (sec.startingTimeStep > 0.0 && (sec.currentTime == sec.initialTime))
      currentTimeStep = std::min(sec.startingTimeStep, currentTimeStep);

    // check for maximum step-size:
    double rh = abs(currentTimeStep)*sec.h_max_inv_; 
    if (rh>1.0)
      currentTimeStep = currentTimeStep/rh;

    sec.currentTimeStep = currentTimeStep;
  }

  sec.currentTimeStepRatio = 1.0;
  sec.currentTimeStepSum   = 2.0*sec.currentTimeStep;

  sec.lastTimeStep      = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio;
  sec.lastTimeStepSum   = sec.currentTimeStepSum;

  sec.numberSuccessiveFailures = 0;
  sec.stepAttemptStatus        = true;

  if (VERBOSE_TIME && tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    Xyce::dout() << "ERROROPTION=1:  DeltaT Grow = 2" << "\n" << std::endl
      << "ERROROPTION=1:  DeltaT Cut = 0.125" << "\n" << std::endl
      << "ERROROPTION=1:  NL MIN = " << tia_params.NLmin << "\n" << std::endl
      << "ERROROPTION=1:  NL MAX = " << tia_params.NLmax << "\n" << std::endl
      << "ERROROPTION=1:  DELMAX = " << sec.maxTimeStep << "\n" << std::endl;
  }

  //  sec.tolAimFac_ = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

  //  if (sec.currentTime == sec.initialTime)
  {
    // x history
    *(ds.xHistory[0]) = *(ds.currSolutionPtr);
    *(ds.xHistory[1]) = *(ds.currSolutionPtr);

    // q history
    *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
    //  *(ds.qHistory[1]) = *(ds.daeFVectorPtr);
    //  (ds.qHistory[1])->scale(-sec.currentTimeStep);
    //  (ds.qHistory[1])->putScalar(0.0);
    *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

    // state history
    *(ds.sHistory[0]) = *(ds.currStatePtr);
    (ds.sHistory[1])->putScalar(0.0);

    // lead current Q compontent history
    *(ds.stoLeadCurrQHistory[0]) = *(ds.currStoreLeadCurrQPtr);
    *(ds.stoLeadCurrQHistory[1]) = *(ds.currStoreLeadCurrQPtr);

    // store history
    *(ds.stoHistory[0]) = *(ds.currStorePtr);
    (ds.stoHistory[1])->putScalar(0.0);

    // lead current history
    *(ds.leadCurrentHistory[0]) = *(ds.currLeadCurrentPtr);
    (ds.leadCurrentHistory[1])->putScalar(0.0);
    *(ds.leadCurrentQHistory[0]) = *(ds.currLeadCurrentQPtr);
    *(ds.leadCurrentQHistory[1]) = *(ds.currLeadCurrentQPtr);
    *(ds.leadDeltaVHistory[0]) = *(ds.currLeadDeltaVPtr);
    (ds.leadDeltaVHistory[1])->putScalar(0.0);
  }

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.currentOrder_ = 1;
  //  sec.usedOrder_ = 1;
  //  if (sec.currentTime == sec.initialTime)
  //  {
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  //  }
  sec.nscsco_ = 0;
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::initialize" << std::endl
      << "\n xHistory: \n" << std::endl;
    (ds.xHistory[0])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.xHistory[1])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n qHistory: \n" << std::endl;
    (ds.qHistory[0])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.qHistory[1])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n sHistory: \n" << std::endl;
    (ds.sHistory[0])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.sHistory[1])->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << std::endl
      << "\n" << "time_to_stop = " << time_to_stop << "\n" << std::endl
      << Xyce::section_divider << std::endl;
  }

  initializeSensitivities();
}

//-----------------------------------------------------------------------------
// Function      : Gear12::initializeSensitivities
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::initializeSensitivities()
{
  int numParams = ds.sensRHSPtrVector.size();
  for (int ip=0; ip<numParams;++ip)
  {
    std::vector<Linear::Vector*> & dqdpHistory = ds.dqdpHistory[ip];
    std::vector<Linear::Vector*> & dXdpHistory = ds.dXdpHistory[ip];
    std::vector<Linear::Vector*> & dQdxdXdpHistory = ds.dQdxdXdpHistory[ip];

    Linear::Vector * currDQdxDxdpPtr = ds.currDQdxDXdpPtrVector[ip];
    Linear::Vector * currDxdpPtr = ds.currDXdpPtrVector[ip];
    Linear::Vector * currDqdpPtr = ds.currDqdpPtrVector[ip];

    // dQdx*dXdp matvec history
    *(dQdxdXdpHistory[0]) = *(currDQdxDxdpPtr);
    *(dQdxdXdpHistory[1]) = *(currDQdxDxdpPtr);

    // dXdp history
    *(dXdpHistory[0]) = *(currDxdpPtr);
    *(dXdpHistory[1]) = *(currDxdpPtr);

    // dqdp history
    *(dqdpHistory[0]) = *(currDqdpPtr);
    *(dqdpHistory[1]) = *(currDqdpPtr);
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::setTwoLevelTimeInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
void
Gear12::setTwoLevelTimeInfo()
{
  // set initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;


  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  *(ds.xHistory[1]) = *(ds.currSolutionPtr); 

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0); 

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::checkReduceOrder()
// Purpose       : check whether to reduce order independent of local error test
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::checkReduceOrder()
{

}

//-----------------------------------------------------------------------------
// Function      : Gear12::rejectStep()
// Purpose       : code to restore history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void
Gear12::rejectStep(
  const TIAParams &     tia_params)
{
  // This routine puts its output in newTimeStep_ and sec.newOrder_

  // This routine changes the following variables:
  //    lastAttemptedTimeStep, sec.initialPhase_, sec.nef_, sec.psi_, newTimeStep_,
  //    sec.newOrder_, sec.currentOrder_, currentTimeStep_, ds.xHistory,
  //    ds.qHistory, nextTimePt, nextTime, currentTimeStepRatio,
  //    currentTimeStepSum, nextTimePt

  // This routine reades but does not change the following variables:
  //    stepAttemptStatus, sec.r_factor_, sec.r_safety_, sec.Est_, sec.r_fudge_, sec.r_min_, sec.r_max_,
  //    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep

  sec.TimeStepLimitedbyBP = false;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::rejectStep" << std::endl;
  }

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = !tia_params.constantTimeStepFlag;

  sec.lastAttemptedTimeStep = sec.currentTimeStep;


  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
    {
      newTimeStep_ = sec.currentTimeStep/8;
    }
    else
    {
      sec.initialPhase_ = false;
      sec.nef_++;
      restoreHistory();
      // restore sec.psi_
      //    for (int i=1;i<=sec.currentOrder_;++i)
      //      sec.psi_[i-1] = sec.psi_[i] - sec.currentTimeStep;

      if (sec.nef_ >= sec.max_LET_fail_)  
      {
        std::string msg = "Gear12::rejectStep: ";
        msg += "  Maximum number of local error test failures.  ";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }

      if ((sec.newtonConvergenceStatus <= 0))
      {
        /// 11/11/05 erkeite:  If the Newton solver fails, don't 
        // rely on the error estimate - it may be full of Nan's.
        //        rr = sec.r_min_;

        newTimeStep_ = sec.currentTimeStep/8;
        //        sec.currentOrder_ = 1;
        sec.currentOrder_ = sec.minOrder_;
        //      if (sec.nef_ > 2) sec.newOrder_ = 1;//consistent with block below.
      }
      else
      {
        // 03/11/04 tscoffe:  Here is the block for choosing order & 
        // step-size when the local error test FAILS (but Newton 
        // succeeded). 
        if (sec.nef_ == 1) // first local error test failure
        {
          //	sec.estOverTol_ 
          sec.Est_ = sec.estOverTol_;

          rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
          rr = pow(rr, 1.0/(sec.currentOrder_+1.0));
          rr = std::max(sec.r_min_,std::min(sec.r_max_,rr));

          newTimeStep_ = rr * sec.currentTimeStep;

          //          sec.currentOrder_ = 1; 
          //          sec.currentOrder_ = sec.minOrder_;
        }
        else // if (sec.nef_ == 2) // second Dae failure
        {
          rr = sec.r_min_;
          newTimeStep_ = rr * sec.currentTimeStep;

          //          sec.currentOrder_ = 1;
          sec.currentOrder_ = sec.minOrder_;
        }
      }
      if (DEBUG_TIME && isActive(Diag::TIME_STEP) && isActive(Diag::TIME_COEFFICIENTS))
      {
        Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
          << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl
          << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl
          << "  nscsco_ = " <<  sec.nscsco_ << std::endl
          << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl
          << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl
          << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl
          << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl
          << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl
          << "  psi_[0] = " <<  sec.psi_[0] << std::endl
          << "  psi_[1] = " <<  sec.psi_[1] << std::endl
          << "  psi_[2] = " <<  sec.psi_[2] << std::endl
          << "  psi_[3] = " <<  sec.psi_[3] << std::endl
          << "  psi_[4] = " <<  sec.psi_[4] << std::endl
          << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl
          << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl
          << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl
          << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl
          << "  sigma_[4] = " <<  sec.sigma_[4] << std::endl
          << "  rr = " <<  rr << std::endl
          << "  r_factor_ = " <<  sec.r_factor_ << std::endl
          << "  r_safety_ = " <<  sec.r_safety_ << std::endl
          << "  Est_ = " <<  sec.Est_ << std::endl
          << "  r_fudge_ = " <<  sec.r_fudge_ << std::endl
          << "  newOrder_ = " <<  sec.newOrder_ << std::endl
          << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
          << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      }
    }
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    std::string tmp = "  Gear12:rejectStep: Warning: Local error test failed with constant step-size.\n";
    Xyce::dout() << tmp << std::endl;
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep_ = std::max(newTimeStep_, sec.minTimeStep);
    newTimeStep_ = std::min(newTimeStep_, sec.maxTimeStep);

    double nextTimePt = sec.currentTime + newTimeStep_;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt  = sec.stopTime;
      newTimeStep_ = sec.stopTime - sec.currentTime;
      sec.TimeStepLimitedbyBP = true;
    }

    sec.nextTime = nextTimePt;

    sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
    sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl
        << "  nextTime = " <<  sec.nextTime << std::endl;
    }

    sec.currentTimeStep = newTimeStep_;
  }
  else // if time step is constant for this step:
  {
    double nextTimePt = sec.currentTime + sec.currentTimeStep;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt      = sec.stopTime;
      sec.currentTimeStep = sec.stopTime - sec.currentTime;
    }

    sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
    sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

    sec.nextTime = nextTimePt;
  }
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::rejectStepForHabanero
// Purpose       : step rejection, but from an outside program (Habanero API)
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/11/09
//-----------------------------------------------------------------------------
void Gear12::rejectStepForHabanero()
{
  restoreHistory();
  sec.setTimeStep(sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : Gear12::completeStep()
// Purpose       : code to update history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::completeStep(const TIAParams &tia_params)
{
  sec.TimeStepLimitedbyBP = false;

  sec.numberOfSteps_ ++;
  sec.nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::completeStep" << std::endl;
  }

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = !tia_params.constantTimeStepFlag;

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  int orderDiff = sec.currentOrder_ - sec.usedOrder_;
  sec.usedOrder_ = sec.currentOrder_;
  sec.usedStep_ = sec.currentTimeStep;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "  initialPhase_ = " <<  sec.initialPhase_ << std::endl
      << "  rr = " <<  rr << std::endl
      << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
      << "  currentTime = " <<  sec.currentTime << std::endl
      << "  nextTime = " <<  sec.nextTime << std::endl
      << "  newTimeStep_ = " <<  newTimeStep_ << std::endl
      << "  minTimeStep = " <<  sec.minTimeStep << std::endl
      << "  maxTimeStep = " <<  sec.maxTimeStep << std::endl;
  }

  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.


  /*   if (sec.numberOfSteps_ >= 2)
       {
       sec.currentOrder_ = 2;
  //     (ds.relErrTolPtr)->putScalar(1e-2);	
  } 
  */ 
  if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {

    if (sec.numberOfSteps_ >= 2 &&  sec.maxOrder_ == 2)
    {
      sec.currentOrder_ = 2;   
    } 

    rr = 1;

    if (sec.nIterations <= tia_params.NLmin)
      rr = 2;  

    if (sec.nIterations > tia_params.NLmax)
      rr = 1.0/8;

    newTimeStep_ = rr*sec.currentTimeStep;
  }
  else
  {
    rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
    rr = pow(rr, 1.0/(sec.currentOrder_+1.0));

    if (sec.numberOfSteps_ >= 2 && sec.maxOrder_ == 2)
    {
      if (sec.currentOrder_ == 1)
      { 
        sec.currentOrder_ = 2;   
        rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
        rr = pow(rr, 1.0/(sec.currentOrder_+1.0)); 

        if (rr <= 1.05)
        {
          //            sec.currentOrder_ = 1; 
          sec.currentOrder_ = sec.minOrder_;
        }
      }
    } 
    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl;
      Xyce::dout() << "  r_safety = " <<  sec.r_safety_ << std::endl;
      Xyce::dout() << "  r_fudge_ = " <<  sec.r_fudge_ << std::endl;
      Xyce::dout() << "  r_hincr_ = " <<  sec.r_hincr_ << std::endl;
      Xyce::dout() << "  r_hincr_test_ = " <<  sec.r_hincr_test_ << std::endl;
      Xyce::dout() << "  Est = " <<  sec.Est_ << std::endl;
      Xyce::dout() << "  raw rr = " <<  rr << std::endl;
    }

    if (rr >= sec.r_hincr_test_)
    {
      rr = sec.r_hincr_;
      newTimeStep_ = rr*sec.currentTimeStep;
    }
    else if (rr <= 1)
    {
      rr = std::max(sec.r_min_,std::min(sec.r_max_,rr));
      newTimeStep_ = rr*sec.currentTimeStep;
    }
  }

  updateHistory();

  newTimeStep_ = std::max(newTimeStep_, sec.minTimeStep);
  newTimeStep_ = std::min(newTimeStep_, sec.maxTimeStep);

  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  //  if (sec.currentTime < sec.stopTime)
  if ((sec.stopTime - sec.currentTime) >= sec.minTimeStep) 
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      //      newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
      //      newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

      double nextTimePt = sec.currentTime + newTimeStep_;

      if (nextTimePt > sec.stopTime)
      {

        sec.savedTimeStep = newTimeStep_;

        nextTimePt  = sec.stopTime;
        newTimeStep_ = sec.stopTime - sec.currentTime;
        sec.TimeStepLimitedbyBP = true;
      }

      sec.nextTime = nextTimePt;

      sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
      sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

      if (DEBUG_TIME && isActive(Diag::TIME_STEP))
      {
        Xyce::dout() << "  nextTime = " <<  sec.nextTime << std::endl;
        Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      }


    sec.currentTimeStep = newTimeStep_;
    }
    else // if time step is constant for this step:
    {
      double nextTimePt = sec.currentTime + sec.currentTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        sec.savedTimeStep = sec.currentTimeStep;

        nextTimePt      = sec.stopTime;
        sec.currentTimeStep = sec.stopTime - sec.currentTime;
      }

      sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
      sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

      sec.nextTime = nextTimePt;
    }
  }

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

//  sec.currentTimeStep = newTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateStateDeriv
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
void Gear12::updateStateDeriv ()
{
  ds.nextStateDerivPtr->
    linearCombo(sec.alpha_[0],*ds.nextStatePtr, sec.alpha_[1],*ds.sn0Ptr);

  if (sec.currentOrder_ == 2)
  {
    ds.nextStateDerivPtr->
      linearCombo(1.0,*ds.nextStateDerivPtr, sec.alpha_[2], *(ds.sHistory[1]));
  }

  ds.nextStateDerivPtr->scale(1.0/sec.currentTimeStep);

  if (DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
  {
    Xyce::dout() << "\n next state Ptr: \n" << std::endl;
    ds.nextStatePtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;

    Xyce::dout() << "\n sn0: \n" << std::endl;
    ds.sn0Ptr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;

    Xyce::dout() << "\n next State Deriv: \n" << std::endl;
    ds.nextStateDerivPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateLeadCurrent
// Purpose       : calculates lead currents in store vector with 
//                 the storeLeadCurrQVec. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 03/22/2013
//-----------------------------------------------------------------------------
void Gear12::updateLeadCurrent ()
{
  ds.nextStoreLeadCurrQDerivPtr->
    linearCombo(sec.alpha_[0],*ds.nextStoreLeadCurrQPtr, sec.alpha_[1],*ds.stoQn0Ptr);

  if (sec.currentOrder_ == 2)
  {
    ds.nextStoreLeadCurrQDerivPtr->
      linearCombo(1.0,*ds.nextStoreLeadCurrQDerivPtr, sec.alpha_[2], *(ds.stoLeadCurrQHistory[1]));
  }

  ds.nextStoreLeadCurrQDerivPtr->scale(1.0/sec.currentTimeStep);

  ds.nextStorePtr->linearCombo(1.0,*ds.nextStorePtr,1.0,*ds.nextStoreLeadCurrQDerivPtr);

  if (DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
  {
    Xyce::dout() << "\n next store Ptr: \n" << std::endl;
    ds.nextStorePtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateLeadCurrentVec
// Purpose       : calculates lead currents in lead current vector with 
//                 the leadCurrQVec. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 03/22/2013
//-----------------------------------------------------------------------------
void Gear12::updateLeadCurrentVec ()
{
  ds.nextLeadCurrentQDerivPtr->
    linearCombo(sec.alpha_[0],*ds.nextLeadCurrentQPtr, sec.alpha_[1],*ds.leadCurrentQn0Ptr);

  if (sec.currentOrder_ == 2)
  {
    ds.nextLeadCurrentQDerivPtr->
      linearCombo(1.0,*ds.nextLeadCurrentQDerivPtr, sec.alpha_[2], *(ds.leadCurrentQHistory[1]));
  }

  ds.nextLeadCurrentQDerivPtr->scale(1.0/sec.currentTimeStep);

  ds.nextLeadCurrentPtr->linearCombo(1.0,*ds.nextLeadCurrentPtr,1.0,*ds.nextLeadCurrentQDerivPtr);

  if (DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
  {
    Xyce::dout() << "\n next leadcurrent q Ptr: \n" << std::endl;
    ds.nextLeadCurrentPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : Gear12::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void Gear12::getInitialQnorm(TwoLevelError & tle) const
{
  tle.q1HistorySum = ds.partialSum_q1();
}

//-----------------------------------------------------------------------------
// Function      : Gear12::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void Gear12::getTwoLevelError(TwoLevelError & tle) const
{
  tle.xErrorSum    = ds.partialErrorNormSum ();
  tle.qErrorSum    = ds.partialQErrorNormSum ();
  tle.xErrorSum_m1 = ds.partialSum_m1 (sec.currentOrder_);
  tle.xErrorSum_m2 = ds.partialSum_m2 (sec.currentOrder_);
  tle.xErrorSum_p1 = ds.partialSum_p1 (sec.currentOrder_, sec.maxOrder_);
  tle.innerSize    = ds.globalLength ();

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
    Xyce::dout() << tle;
}

double Gear12::partialTimeDeriv() const
{
  if (sec.currentTimeStep < 1e-30) 
  {
    Xyce::Report::UserWarning() << "Excessively small current time step, incorrectly returning with large value";

    return leadingCoeff * 1.e+30;
  }
  
  return leadingCoeff / sec.currentTimeStep;
}

} // namespace TimeIntg
} // namespace Xyce
