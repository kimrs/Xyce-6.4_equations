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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_SharedSystem.C,v $
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.44 $
//
// Revision Date  : $Date: 2015/07/11 23:40:30 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_SharedSystem.h"
#include "N_NLS_NOX_Interface.h"
#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "N_LAS_System.h"
#include "N_LAS_Builder.h"
#include "N_ERH_ErrorMgr.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include <N_NLS_NOX_Group.h>
#include <NOX_Epetra_Group.H>
#include <N_NLS_LOCA_Group.h>

// ----------   NOX Includes   ----------

// ---------- Namespaces ---------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

/* EDIT KIM */
void print(const Epetra_MultiVector *vec, const char * filename) {
    FILE * file = fopen(filename, "w");
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        fprintf(file, "%d,%d,%0.16e\n", i, 0, values[i]);
    
    fclose(file);
}

void print(const Epetra_CrsMatrix * A, const char * filename) {
   int n_rows = A->NumMyRows();
   int n_maxCols = A->NumGlobalCols();
   FILE * file = fopen(filename, "w");
   for(int i = 0; i < n_rows; ++i) {
       int n_entries;
       double * values = new double[n_maxCols];
       int * indices = new int[n_maxCols];
       A->ExtractMyRowCopy(i, n_rows, n_entries, values, indices); 
       for(int j = 0; j < n_entries; ++j) 
         fprintf(file, "%d,\t%d,\t%0.16f\n", i, indices[j], values[j]); 

       delete values;
       delete indices;
   }
   fclose(file);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::SharedSystem
// Purpose       : Constructor. Creates a shared system containing
//                 the soln vector, the previous solution vector,
//                 the RHS vector, the Newton vector, the Jacobian
//                 matrix, and a reference to the interface used to
//                 call the evaluation functions.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::SharedSystem(Linear::Vector& soln,
			   Linear::Vector& f,
			   Linear::Matrix& jacobian,
			   Linear::Vector& newton,
			   Linear::Vector& gradient,
			   Linear::System& lasSys,
			   Interface& interface) :

  xyceSolnPtr_(0),
  xyceFPtr_(0),
  xyceJacobianPtr_(0),
  xyceNewtonPtr_(0),
  xyceGradientPtr_(0),
  xyceLasSysPtr_(0),
  xyceInterfacePtr_(0),
  matrixFreeFlag_(false),
  ownerOfJacobian_(0),
  ownerOfStateVectors_(0),
  ifpackGraphPtr_(0),
  ifpackPreconditionerPtr_(0)
{
  reset(soln, f, jacobian, newton, gradient, lasSys, interface);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::~SharedSystem
// Purpose       : Destructor.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::~SharedSystem()
{
  deletePreconditioner();
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;
}

/* EDIT KIM */
void print(Epetra_MultiVector *vec, char * syntax, int disambiguation) {
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        printf(syntax, i, 0, values[i]);
}

void print(Epetra_CrsMatrix * A, char * syntax, int disambiguation) {
   int n_rows = A->NumMyRows();
   int n_maxCols = A->NumGlobalCols();
   for(int i = 0; i < n_rows; ++i) {
       int n_entries;
       double * values = new double[n_maxCols];
       int * indices = new int[n_maxCols];
       A->ExtractMyRowCopy(i, n_rows, n_entries, values, indices); 
       for(int j = 0; j < n_entries; ++j) 
         printf(syntax, i, indices[j], values[j]); 

       delete values;
       delete indices;
   }
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::reset
// Purpose       : reset the Xyce fill objects - pointers may have changed!
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::reset(Linear::Vector& x,
			 Linear::Vector& f,
			 Linear::Matrix& jacobian,
			 Linear::Vector& newton,
			 Linear::Vector& gradient,
			 Linear::System& lasSys,
			 Interface& interface)
{
  // Clear out old views
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;

  xyceJacobianPtr_ = &jacobian;
  xyceLasSysPtr_ = &lasSys;
  xyceInterfacePtr_ = &interface;
  matrixFreeFlag_ = xyceInterfacePtr_->getMatrixFreeFlag();

  // Create views of the data used for fills in xyce
  xyceSolnPtr_ = new Vector(x, lasSys);
  xyceFPtr_ = new Vector(f, lasSys);
  xyceNewtonPtr_ = new Vector(newton, lasSys);
  xyceGradientPtr_ = new Vector(gradient, lasSys);

  /* EDIT KIM */
  Xyce::Linear::Vector * vec = xyceSolnPtr_->getNativeVectorPtr();
  print(&vec->epetraObj(),"%d,%d,%0.64e\n", 0);
  std::cout << "done\n"; /* soln has changed */


  // Wipe the preconditioner clean
  deletePreconditioner();
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeF
// Purpose       : Compute the F corresponding to the current
//                 primary solution vector. Makes the primary
//                 solution vector owner in to the owner of the F.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
/* EDIT KIM */
std::string genFileName(std::string filename) {
  const std::string PATH = "csv_files/";
  const std::string SUFFIX = ".csv";
  std::stringstream ss; 
  ss << PATH;
  ss << filename;
  ss << "_";

  std::ifstream fin("curr_timestep");
  std::copy(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>(), std::ostreambuf_iterator<char>(ss));
  ss << SUFFIX;
  std::string str = ss.str();
  str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
  return str;
};

bool SharedSystem::computeF(const Vector& solution, Vector& F,
			    const Group* grp)
{
  ownerOfStateVectors_ = grp;
  

 // solution.getNativeVectorRef().printPetraObject(std::cout); /* EDIT KIM */

#ifdef Xyce_NOX_USE_VECTOR_COPY
  *xyceSolnPtr_ = solution;
  bool status = xyceInterfacePtr_->computeF();
#else
  bool status = xyceInterfacePtr_->computeF(F, solution);
#endif

//  solution.getNativeVectorRef().printPetraObject(std::cout); /* EDIT KIM */
  //solution.printPetraObject(std::cout); 
  ///* EDIT KIM */
 // print(&solution.getNativeVectorRef().epetraObj(), genFileName("sol").c_str());
 // print(&F.getNativeVectorRef().epetraObj(), genFileName("F").c_str());

  if (status == false) {
    const std::string message = "Error: SharedSystem::computeF() - compute F failed!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

#ifdef Xyce_NOX_USE_VECTOR_COPY
  F = *xyceFPtr_;
#endif
  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeJacobian
// Purpose       : Compute the Jacobian corresponding to the current
//                 primary solution vector. 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeJacobian(Group* grp)
{
  ownerOfJacobian_ = grp;

#ifdef Xyce_NOX_USE_VECTOR_COPY
  *xyceSolnPtr_ = grp->getX();
#endif

  if (!areStateVectors(grp)) {
#ifdef Xyce_VERBOSE_NOX
    if (1) { //RPP: Need to add priting utilities to group ctor
      dout() << "Warning: SharedSystem::computeJacobian() - State "
	   << "Vectors are not valid wrt solution!" << std::endl;
      dout() << "Calling computeResidual to fix this!" << std::endl;
    }
#endif
    // RPP: This line is not needed since we now call the group
    //ownerOfStateVectors_ = grp;
    
    NOX::Abstract::Group::ReturnType status = grp->computeF();

    if (status != NOX::Abstract::Group::Ok) {
      const std::string message = "SharedSystem::computeJacobian() - compute F failed!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
    }

  }

#ifdef Xyce_NOX_USE_VECTOR_COPY
  bool status = xyceInterfacePtr_->computeJacobian();
#else
  bool status = xyceInterfacePtr_->computeJacobian
    (dynamic_cast<const Vector &> (grp->getX()));
#endif

  if (status == false) {
    const std::string message = "SharedSystem::computeJacobian() - compute Jac failed!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }


  /* EDIT KIM */
//  print(&xyceJacobianPtr_->epetraObj(), genFileName("Jac").c_str());

  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeNewton
// Purpose       : Compute the Newton corresponding to the current
//                 primary solution vector. 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeNewton(const Vector& F, Vector& Newton,
				 Teuchos::ParameterList& params)
{
  /* EDIT KIM */
  std::ofstream ofs;
  ofs.open("computeNewton_F.txt"); /* F has some stuff removed here */
  Xyce::Nonlinear::N_NLS_LOCA::Group soln(*(xyceInterfacePtr_->getSolutionGroup()));
  soln.getF().print(ofs);
  ofs.close();
  /* END KIM */

  *xyceFPtr_ = F;
  // Zero out the Newton vector
  xyceNewtonPtr_->scale(0.0);
  bool status = xyceInterfacePtr_->computeNewton(params);
  Newton = *xyceNewtonPtr_;

  /* EDIT KIM */
  ofs.open("curr_Newton.txt");
  Newton.print(ofs);
  ofs.close();
  /* END KIM */

  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeGradient
// Purpose       : Compute the Gradient corresponding to the current
//                 primary solution vector. 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeGradient(const Vector& F, Vector& Gradient)
{
  *xyceFPtr_ = F;
  xyceGradientPtr_->scale(0.0);
  bool status = xyceInterfacePtr_->computeGradient();
  Gradient = *xyceGradientPtr_;
  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::applyJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobian(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool NoTranspose = false;
    xyceJacobianPtr_->matvec(NoTranspose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
    // tscoffe/tmei HB 07/29/08
#ifndef Xyce_NOX_USE_VECTOR_COPY
    const std::string message = "SharedSystem::applyJacobian() - ERROR, Xyce_NOX_USE_VECTOR_COPY required";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
#endif
    bool status = xyceInterfacePtr_->applyJacobian(input.getNativeVectorRef(), result.getNativeVectorRef());
    if (status == false) {
      const std::string message = "SharedSystem::applyJacobian() - apply Jac failed!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::applyJacobianTranspose
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool Transpose = true;
    xyceJacobianPtr_->matvec(Transpose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
      const std::string message = "SharedSystem::applyJacobianTranspose() - Not Supported for Matrix Free Loads!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeDfDpMulti	
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool SharedSystem::computeDfDpMulti	
  (const std::vector< int > & paramIDs, 
   NOX::Abstract::MultiVector & dfdp, 
   bool isValidF)
{
  bool status = xyceInterfacePtr_->computeDfDpMulti	
                   (paramIDs, dfdp, isValidF);
    
  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computePreconditioner
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computePreconditioner()
{
  Epetra_CrsMatrix* crs = dynamic_cast<Epetra_CrsMatrix*>
    (&(xyceJacobianPtr_->epetraObj()));

  if (crs == 0) {
    dout() << "SharedSystem::computePreconditioner() - " 
	 << "Dynamic cast to CRS Matrix failed!" << std::endl;
  }
  
  deletePreconditioner();
  ifpackGraphPtr_ = new Ifpack_IlukGraph(crs->Graph(),
					 1,
					 0);
  ifpackGraphPtr_->ConstructFilledGraph();
  ifpackPreconditionerPtr_ = new Ifpack_CrsRiluk(*ifpackGraphPtr_);
  ifpackPreconditionerPtr_->InitValues(*crs);
  ifpackPreconditionerPtr_->Factor();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::deletePreconditioner
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::deletePreconditioner()
{
  delete ifpackPreconditionerPtr_;
  delete ifpackGraphPtr_;
  ifpackPreconditionerPtr_ = 0;
  ifpackGraphPtr_ = 0;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::applyRightPreconditioning
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyRightPreconditioning(bool useTranspose, 
					     Teuchos::ParameterList& params,
					     const Vector& input, 
					     Vector& result)
{
  if (ifpackPreconditionerPtr_ == 0) {
    const std::string message = "SharedSystem::applyRightPreconditioning - Preconditioner is 0!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  if (useTranspose)
    ifpackPreconditionerPtr_->SetUseTranspose(useTranspose);
  
  Linear::Vector& nonConstInput = 
    const_cast<Linear::Vector&>(input.getNativeVectorRef_());
  Epetra_MultiVector& epVecInput = 
    const_cast<Epetra_MultiVector&>(nonConstInput.epetraObj());
  
  Linear::Vector& nonConstResult = 
    const_cast<Linear::Vector&>(result.getNativeVectorRef_());
  Epetra_MultiVector& epVecResult = 
    const_cast<Epetra_MultiVector&>(nonConstResult.epetraObj());
  
  int errorCode = ifpackPreconditionerPtr_->
    ApplyInverse(epVecInput, epVecResult);
  
  // Unset the transpose call
  if (useTranspose)
    ifpackPreconditionerPtr_->SetUseTranspose(false);    

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Vector& SharedSystem::getSolutionVector()
{
  return *xyceSolnPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const Linear::Matrix& SharedSystem::getJacobian() const
{
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Linear::Matrix& SharedSystem::getJacobian(const Group* grp)
{
  ownerOfJacobian_ = grp;
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getStateVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::getStateVectors(const Group* grp)
{
  ownerOfStateVectors_ = grp;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getLasSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Linear::System* SharedSystem::getLasSystem()
{
  return xyceLasSysPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::cloneSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Vector* SharedSystem::cloneSolutionVector() const
{
  Vector* tmpVectorPtr = 0;
  tmpVectorPtr = 
    dynamic_cast<Vector*>(xyceSolnPtr_->clone(NOX::DeepCopy).release().get());

  if (tmpVectorPtr == 0) {
    const std::string message = 
      "SharedSystem::cloneSolutionVector() - dynamic cast/ memory allocation failure!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return (tmpVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem:: getNewtonVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const Vector & SharedSystem::getNewtonVector() const
{
  return *xyceNewtonPtr_;                                                       
}                                                                               

//-----------------------------------------------------------------------------
// Function      : SharedSystem::debugOutput1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput1
   (Linear::Matrix & jacobian, Linear::Vector & rhs)
{
  xyceInterfacePtr_->debugOutput1(jacobian, rhs);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::debugOutput3
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput3 
   (Linear::Vector & dxVector, Linear::Vector & xVector)
{
  xyceInterfacePtr_->debugOutput3(dxVector, xVector);
}

}}}
