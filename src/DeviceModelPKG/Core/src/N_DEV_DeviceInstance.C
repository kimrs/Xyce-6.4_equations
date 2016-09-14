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
// Filename       : $RCSfile: N_DEV_DeviceInstance.C,v $
//
// Purpose        : Implementation of the base device instance class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/30/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.210 $
//
// Revision Date  : $Date: 2015/07/02 21:29:42 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceInstance.h>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_NumericalJacobian.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>


/* EDIT KIM */
#include <Epetra_CrsMatrix.h>
#include <N_DEV_DeviceMgr.h>
/* DONE KIM */

namespace Xyce {
namespace Device {

namespace {
static const std::vector<std::string> emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::DeviceInstance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::DeviceInstance(
  const InstanceBlock &         instance_block,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceEntity(parametric_data, factory_block.solverState_, factory_block.deviceOptions_, instance_block.getNetlistLocation().getFilename(), instance_block.getNetlistLocation().getLineNumber()),
    name_(instance_block.getInstanceName()), 
    mlData(factory_block.matrixLoadData_),
    extData(factory_block.externData_),
    configuredForLeadCurrent(false),
    cols(factory_block.matrixLoadData_.cols),
    vals(factory_block.matrixLoadData_.vals),
    numJacPtr(NULL),
    psLoaded(false),
    ssLoaded(false),
    rhsLoaded(false),
    origFlag(true),
    numIntVars(0),
    numExtVars(2),
    numStateVars(0),
    numStoreVars(0),
    numLeadCurrentVars(0),
    numLeadCurrentStoreVars(0),
    loadLeadCurrent(false),
    numBranchDataVars(0),
    numBranchDataVarsIfAllocated(0),
    mergeRowColChecked(false)
{
  devConMap.resize(2);
  devConMap[0] = 1;
  devConMap[1] = 1;
  numJacPtr = new NumericalJacobian(factory_block.matrixLoadData_, factory_block.solverState_, factory_block.externData_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::~DeviceInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::~DeviceInstance()
{
  delete numJacPtr;
}

/* EDIT KIM */
int debug_func(int g) {
    while(g < -3) { /* break here */
        --g;
    }
    return g;
}

/* There appears to be a problem with using numrows as length when adding the strings. 
 */
void 
DeviceInstance::print_strings(Linear::Matrix &matrix, EquationStrings* eqStrings, std::string name = "") {
  int n = matrix.epetraObj().NumMyRows();

  for(int i = 0; i < n; ++i) {
      int row = matrix.epetraObj().GRID64(i);
      int n_entries;
      double values[n];
      int indices[n]; 
      matrix.epetraObj().ExtractGlobalRowCopy(row, n, n_entries, values, indices); 
    
      /* Check if nonzero then write to text matrix */
      for(int j = 0; j < n_entries; ++j) {
          if(values[j]) {
//              std::cout << "adding " << values[j] << " to (num " << i <<  ", " << indices[j] << ")\n";
              std::stringstream  * ss = eqStrings->getStringStream(i, indices[j]);
              *ss << "\t+(";
              *ss << values[j];
              *ss << "):";
              *ss << name;
          }
      }
  }
}

void
DeviceInstance::print_strings(Linear::MultiVector &matrix, EquationStrings* eqStrings, std::string name = "") {

}
/* DONE KIM */

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepSolnVars
//
// Purpose       : This function configures a device to an auxiliary F & Q vector
//                 so that lead currents can be calculated for this device.c.
//
// Special Notes : This must be called soon after the constructor call
//                 before the store vector is allocated.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/22/13
//-----------------------------------------------------------------------------
void DeviceInstance::enableLeadCurrentCalc()
{
  if (!configuredForLeadCurrent)
  {
    // indicated that this device is now configured for lead current calculation
    // this avoids claiming too much space in the store vector if this function
    // is called more than once for a device.
    configuredForLeadCurrent = true;

    // set device instance flag to indicate the need to load lead current
    // data into store F & Q vectors
    loadLeadCurrent = true;

    // request additional space in store vector for lead current vars
    numStoreVars = numStoreVars + numLeadCurrentStoreVars;
    
    // migrating the store vector calculation to the branch data vectors
    numBranchDataVars = numBranchDataVars + numBranchDataVarsIfAllocated;
  }
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInstance::getDepSolnVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & DeviceInstance::getDepSolnVars()
{
  return expVarNames;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnLIDs
// Purpose       : Allows registration of LIDs of nodes and instances that
//                 appear in expressions that occur in the device.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/15/05
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnLIDs(
  const std::vector< std::vector<int> > &       depSolnLIDVecRef)
{
  int size = expVarLIDs.size();
  if (size != depSolnLIDVecRef.size())
  {
    DevelFatal0(*this).in("DeviceInstance::registerDepSolnLIDs")
      << "Inconsistent number of LIDs returned from topology";
  }
  for (int i = 0; i < size; ++i)
  {
    if (depSolnLIDVecRef[i].size() != 1)
    {
      UserError0(*this) << "Problem with value for " << expVarNames[i]
                        << ".  This may be an incorrect usage of a lead current in place of a current through a voltage source.";
    }
    expVarLIDs[i] = depSolnLIDVecRef[i][0];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnGIDs(
  const std::vector< std::vector<int> > & varList )
{
  int size = expVarGIDs.size();
  for (int i = 0; i < size; ++i)
  {
    expVarGIDs[i] = varList[i][0];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerJacLIDs( const JacobianStamp & jacLIDVec )
{
  if (getDeviceOptions().numericalJacobianFlag || getDeviceOptions().testJacobianFlag)
  {
    devJacLIDs = jacLIDVec;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepStateVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & DeviceInstance::getDepStateVars()
{
  return emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepStateGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepStateGIDs(
  const std::vector< std::vector<int> > & varList )
{
  if( varList.size() != 0 )
  {
    Report::DevelFatal().in("DeviceInstance::registerDepStateGIDs")
      << "Call to registerDepStateGIDs for a device which doesn't use it.";
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepStoreVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
const std::vector<std::string> & DeviceInstance::getDepStoreVars()
{
  return emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepStoreGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepStoreGIDs(
  const std::vector< std::vector<int> > & varList )
{
  if( varList.size() != 0 )
  {
    Report::DevelFatal().in("DeviceInstance::registerDepStoreGIDs")
      << "Call to registerDepStoreGIDs for a device which doesn't use it.";
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepLeadCurrentVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
const std::vector<std::string> & DeviceInstance::getDepLeadCurrentVars()
{
  return emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepLeadCurrentGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepLeadCurrentGIDs(
  const std::vector< std::vector<int> > & varList )
{
  if( varList.size() != 0 )
  {
    Report::DevelFatal().in("DeviceInstance::registerDepLeadCurrentGIDs")
      << "Call to registerDepLeadCurrentGIDs for a device which doesn't use it.";
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::testDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool DeviceInstance::testDAEMatrices(const std::vector<const std::string *> & nameVec)
{
  bool bsuccess = true;

  // if necessary, consolodate the LIDs vector.
  if (devLIDs.empty())
  {
    devLIDs = extLIDVec;
    devLIDs.insert(devLIDs.end(), intLIDVec.begin(), intLIDVec.end());
    devLIDs.insert(devLIDs.end(), expVarLIDs.begin(), expVarLIDs.end());
  }

  bsuccess = numJacPtr->testDAEMatrices(*this, nameVec);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::trivialStampLoader
// Purpose       : This function contains most of the original
//                 loadTrivialMatrixStamp function.  See comments above.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::trivialStampLoader (Linear::Matrix * matPtr)
{
  std::vector<int>::const_iterator firstVar;
  std::vector<int>::const_iterator lastVar;
  std::vector<int>::const_iterator iterVar;

  int localRows = matPtr->getLocalNumRows();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << "Loading trivial stamp for " << getName() << std::endl;
  }

  if (cols.size() < 1) cols.resize(1);
  if (vals.size() < 1) vals.resize(1);

  for (int i = 0; i < 2; ++i)
  {
    // do external vars first.
    if (i==0)
    {
      firstVar = extLIDVec.begin ();
      lastVar  = extLIDVec.end ();
    }
    // then do internal vars.
    else
    {
      firstVar = intLIDVec.begin ();
      lastVar  = intLIDVec.end ();
    }

    for (iterVar=firstVar; iterVar!=lastVar; ++iterVar)
    {
      int row = *iterVar;

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "matrix row = " << row << std::endl;
      }
      if (row < 0)
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "\tNOT loading this one - too small" << std::endl;
        }
        continue;
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "\tloading this one" << std::endl;
      }

      int count = 1;
      vals[0] = 1.0;
      cols[0] = row;

      //matPtr->putLocalRow(row, count, &vals[0], &cols[0]);
      matPtr->replaceLocalRow(row, count, &vals[0], &cols[0]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::loadTrivialDAE_FMatrixStamp
// Purpose       : See loadTrivialMatrixStamp - this is the same thing, except
//                 for the new-DAE F-Matrix.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::loadTrivialDAE_FMatrixStamp ()
{
  return trivialStampLoader (extData.dFdxMatrixPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::zeroMatrixDiagonal
// Purpose       : puts zeros into on matrix diagonal, but just for this
//                 device.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::zeroMatrixDiagonal (Linear::Matrix * matPtr)
{
  std::vector<int>::const_iterator firstVar;
  std::vector<int>::const_iterator lastVar;
  std::vector<int>::const_iterator iterVar;

  int localRows = matPtr->getLocalNumRows();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << "Zeroing the matrix diagonal for " << getName() << std::endl;
  }

  if (cols.size() < 1) cols.resize(1);
  if (vals.size() < 1) vals.resize(1);

  for (int i = 0; i < 2; ++i)
  {
    // do external vars first.
    if (i == 0)
    {
      firstVar = extLIDVec.begin ();
      lastVar  = extLIDVec.end ();
    }
    // then do internal vars.
    else
    {
      firstVar = intLIDVec.begin ();
      lastVar  = intLIDVec.end ();
    }

    for (iterVar=firstVar; iterVar!=lastVar; ++iterVar)
    {
      int row = *iterVar;

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "matrix row = " << row << std::endl;
      }

      if (row < 0) continue;
      if (row >= localRows) continue;

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "\tloading this one" << std::endl;
      }

      int count = 1;
      vals[0] = 0.0;
      cols[0] = row;

      matPtr->putLocalRow(row, count, &vals[0], &cols[0]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::enablePDEContinuation(int &max_PDE_continuation_steps)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::disablePDEContinuation()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationAlpha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationAlpha (double alpha)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationBeta
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/04
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationBeta (double beta)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::setInitialGuess ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double DeviceInstance::getMaxTimeStepSize  ()
{
  return getDeviceOptions().defaultMaxTimeStep;
}


//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap
// Purpose       : Compute Jacobian Stamp and Map for devices that can have merged nodes
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 2/13/04
//-----------------------------------------------------------------------------
void
DeviceInstance::jacStampMap(
  const JacobianStamp & stamp_parent,
  std::vector<int> &    map_parent,
  JacobianStamp &       map2_parent,
  JacobianStamp &       stamp,
  std::vector<int> &    map,
  JacobianStamp &       map2,
  int                   from,
  int                   to,
  int                   original_size)
{
  if (from <= to)
  {
    Report::DevelFatal().in("DeviceInstance::jacStampMap")
      << "From index " << from << " <= " << " to index " << to;
  }

  if (map_parent.size() == 0)
  {
    map_parent.resize(original_size);
    map2_parent.resize(original_size);
    for (int i = 0; i < original_size; ++i)
    {
      map_parent[i] = i;
      map2_parent[i].resize(stamp_parent[i].size());
      for (int j = 0; j < stamp_parent[i].size(); ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }
  map2.resize(original_size);

  // This is to merge the column that is being eliminated into the column it is merged.
  // If extra elements are present then we must increment the map2 value for later
  // elements.  There are multiple cases depending what is populated.
  for (int i = 0; i < original_size; ++i)
  {
    int f_index = -1;
    int t_index = -1;
    int p_row = map_parent[i];
    for (int j = 0; j < stamp_parent[p_row].size(); ++j)
    {
      if (stamp_parent[p_row][j] == from)
        f_index = j;
      if (stamp_parent[p_row][j] == to)
        t_index = j;
    }
    map2[i].resize(map2_parent[i].size());
    int f_index2 = -1;
    int t_index2 = -1;
    for (int j = 0; j < map2_parent[i].size(); ++j)
    {
      map2[i][j] = map2_parent[i][j];
      if (stamp_parent[p_row][map2[i][j]] == from)
        f_index2 = j;
      if (stamp_parent[p_row][map2[i][j]] == to)
        t_index2 = j;
    }
    if (f_index >= 0)
    {
      if (t_index >= 0)
      {
        for (int j = 0; j < map2[i].size(); ++j)
        {
          if (map2[i][j] > f_index)
            map2[i][j]--;
        }
        if (f_index2 >= 0 && t_index2 >= 0)
          map2[i][f_index2] = map2[i][t_index2];
      }
      else
      {
        for (int j = 0; j < map2[i].size(); ++j)
        {
          if (stamp_parent[p_row][map2[i][j]] > to && stamp_parent[p_row][map2[i][j]] < from)
            ++map2[i][j];
        }
        t_index = 0;
        for (int j = 0; j < stamp_parent[p_row].size(); ++j)
        {
          if (to > stamp_parent[p_row][j])
            ++t_index;
          else
            break;
        }
        if (f_index2 >= 0)
          map2[i][f_index2] = t_index;
      }
    }
  }
  map.resize(original_size);
  int p_size = stamp_parent.size();

  // This is to merge the row that is being eliminated into the row it is merged into
  // if extra elements are present then we must increment the map2 value for later
  // elements in the row
  JacobianStamp map2_tmp = map2;
  for (int i = 0; i < stamp_parent[from].size() && stamp_parent[from][i] < p_size-1; ++i)
  {
    bool new_col = false;
    for (int j = 0; j<stamp_parent[to].size(); ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[from][i] == stamp_parent[to][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (int j = 0; j < map2[to].size(); ++j)
      {
        if (stamp_parent[to][map2_tmp[to][j]] > stamp_parent[from][i])
        {
          ++map2[to][j];
        }
      }
    }
  }

  stamp.resize(p_size-1);

  std::vector<int> dup(original_size);

  int f_mod = from;
  for (int i = 1; i < original_size; ++i)
  {
    dup[i] = -1;
    for (int j = 0; j < i; ++j)
    {
      if (map_parent[i] == map_parent[j])
      {
        dup[i] = j;
        if (i <= f_mod)
          ++f_mod;
        break;
      }
    }
  }

  for (int i = 0; i < f_mod; ++i)
    map[i] = map_parent[i];
  map[f_mod] = map[to];
  
  for (int i = f_mod + 1; i < original_size; ++i)
    map[i] = map_parent[i]-1;

  for (int i = 1; i < original_size; ++i)
  {
    if (dup[i] >= 0)
      map[i] = map[dup[i]];
  }


  // Now that we know where the row originally came from, we can do any renumbering
  // needed in the map2 source row
  map2_tmp = map2;
  int map_2_from = from;
  if (map[from] != map[to]) {
    for (int i = from + 1; i < original_size; ++i) {
      if (map[i] == map[to]) {
        map_2_from = i;
        break;
      }
    }
    if (map_2_from == from) {
      Report::DevelFatal().in("DeviceInstance::jacStampMap") << "Internal Error 2";
    }
  }

  for (int i = 0; i < stamp_parent[to].size() && stamp_parent[to][i] < p_size-1; ++i)
  {
    bool new_col = false;
    for (int j = 0; j < stamp_parent[from].size(); ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[to][i] == stamp_parent[from][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (int j = 0; j < map2[map_2_from].size(); ++j)
      {
        if (stamp_parent[from][map2_tmp[map_2_from][j]] > stamp_parent[to][i])
        {
          ++map2[map_2_from][j];
        }
      }
    }
  }

  JacobianStamp fill(p_size);
  for (int i = 0; i < p_size; ++i)
  {
    fill[i].resize(p_size);
    for (int j = 0 ; j < p_size; ++j)
      fill[i][j] = 0;
    for (int j = 0; j < stamp_parent[i].size(); ++j)
      fill[i][stamp_parent[i][j]] = 1;
  }
  for (int i = 0; i < p_size; ++i)
  {
    fill[to][i] += fill[from][i];
  }
  for (int i = 0; i < p_size; ++i)
  {
    fill[i][to] += fill[i][from];
  }
  for (int i = from; i < p_size - 1; ++i)
  {
    for (int j = 0; j < p_size; ++j)
      fill[i][j] = fill[i+1][j];
    for (int j = 0; j < p_size; ++j)
      fill[j][i] = fill[j][i+1];
  }
  for (int i = 0; i < p_size - 1 ; ++i)
  {
    stamp[i].clear();
    for (int j = 0; j < p_size - 1; ++j)
    {
      if (fill[i][j] > 0)
        stamp[i].push_back(j);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap_fixOrder
//
// Purpose       : This function corrects the compressed row column array so
//                 that the column indices are in ascending order.
//
// Special Notes : The reason for this function is that the
//                 DeviceInstance::jacStampMap function
//                 implicitly requires an ordered jacStamp to work correctly.
//                 Some devices (particularly devices that have meshes),
//                 will start out with an non-ordered stamp, at least in the
//                 column entries.
//
//                 Note, this does require the "map" argument, because,
//                 unlike the function DeviceInstance::jacStampMap,
//                 because this function doesn't change the row ordering,
//                 or remove or merge any rows.
//
//                 This function only changes the column ordering in the
//                 compressed row form of the jacStamp.  It thus requires
//                 modifications to map2, which is essentially a column map.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/01/08
//-----------------------------------------------------------------------------
void
DeviceInstance::jacStampMap_fixOrder(
  const JacobianStamp & stamp_parent,
  JacobianStamp &       map2_parent,
  JacobianStamp &       stamp,
  JacobianStamp &       map2)
{
  int current_size = stamp_parent.size();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_JACSTAMP) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Begin DeviceInstance::jacStampMap_fixOrder." << std::endl
                 << Xyce::section_divider << std::endl;
  }

  // if this is the first time this function is called (for a particular device), then
  // allocate the map and set up their trivial contents.
  if (map2_parent.size() == 0)
  {
    map2_parent.resize(current_size);
    for (int i = 0; i < current_size; ++i)
    {
      map2_parent[i].resize(stamp_parent[i].size());
      for (int j = 0; j < stamp_parent[i].size(); ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }

  stamp.clear();
  map2.clear();

  // To make this simple, start out with a full, dense stamp.
  JacobianStamp denseStamp(current_size);
  for (int i = 0; i < current_size; ++i)
  {
    denseStamp[i].resize(current_size,-1);

    for (int j = 0; j < stamp_parent[i].size(); ++j)
    {
      int denseCol = stamp_parent[i][j];
      denseStamp[i][denseCol] = j;
    }
  }

  // At this point, the denseStamp has been set up.  Now use it to re-create the
  // compressed-row stamp.  By simply looping over the dense stamp, the column order
  // in the compressed row stamp will automatically be ascending.
  // map2 is set up here as well, by pulling the values out that we previously put into
  // dense stamp.
  stamp.resize(current_size);
  map2.resize(current_size);
  for (int i = 0; i < current_size; ++i)
  {
    for (int j = 0; j < current_size; ++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1)
      {
        stamp[i].push_back(j);
      }
    }

    int stampRowLength=stamp[i].size();
    map2[i].resize(stampRowLength, 0);

    for (int j = 0, k = 0; j < current_size; ++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1 && colMapIndex < stampRowLength)
      {
        map2[i][colMapIndex] = k;
        ++k;
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_JACSTAMP) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "From inside of DeviceInstance::jacStampMap_fixOrder:" << std::endl
                 << "The original parent stamp is:" << std::endl;
    outputJacStamp(stamp_parent);
    Xyce::dout() << "The new reduced stamp is:" << std::endl;
    outputJacStamp(stamp);
    Xyce::dout() << "The dense stamp is:" << std::endl;
    outputJacStamp(denseStamp);

    Xyce::dout() << "The new map is:" << std::endl;
    outputJacStamp(map2);

    Xyce::dout() << Xyce::section_divider << std::endl
                 << "End DeviceInstance::jacStampMap_fixOrder."<<std::endl
                 << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacStamp
// Purpose       : Output jacStamp (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/19/06
//-----------------------------------------------------------------------------
void
DeviceInstance::outputJacStamp(
  const JacobianStamp &         jac)
{
  for (int i = 0 ; i < jac.size(); ++i)
  {
    Xyce::dout() << "Row: " << i << " ::";
    for (int j = 0 ; j < jac[i].size(); ++j)
      Xyce::dout() << "  " << jac[i][j];

    Xyce::dout() << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacMaps
// Purpose       : Output jacMap and jacMap2 (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, Electrical & Microsystems Modeling
// Creation Date : 02/20/08
//-----------------------------------------------------------------------------
void
DeviceInstance::outputJacMaps(
  const std::vector<int>  &     jacMap,
  const JacobianStamp &         jacMap2)
{
  for (int i = 0 ; i < jacMap.size(); ++i)
  {
    Xyce::dout() << "Row " << i << ": ";
    for (int j = 0; j < jacMap2[i].size(); j++)
      Xyce::dout() << jacMap[i]<< "," << jacMap2[i][j] << " ";

    Xyce::dout() << std::endl;
  }

  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerGIDData
// Purpose       : Insert GID data into 'indexPairList' object
//
// Special Notes : This information is neccessary for the numerical Jacobian.
//                 It is only called once.
//
//                 The numerical jacobian may get confused by duplicate
//                 matrix entries, in that it might load them 2x.  For that
//                 reason, this function checks for duplicates.
//
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 12/13/04
//-----------------------------------------------------------------------------
void DeviceInstance::registerGIDData(
  const std::vector<int> &      counts,
  const IdVector &              GIDs,
  const JacobianStamp &         jacGIDs )
{
  if (getDeviceOptions().numericalJacobianFlag)
  {
    indexPairList.clear();

    int extSize = counts[0];
    int intSize = counts[1];
    int expSize = counts[2];
    int size = GIDs.size();

    std::map<int,int> testIndexMap;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "DeviceInstance::registerGIDData for " << getName() << std::endl
                   << "  extSize      = " << extSize  << std::endl
                   << "  intSize      = " << intSize  << std::endl
                   << "  expSize      = " << expSize  << std::endl
                   << "  GIDs.size    = " << size << std::endl
                   << "  jacGIDs.size = " << jacGIDs.size () << std::endl
                   << std::endl;
    }

    {
      int gid_index = 0;

      // Copy out external gids:
      extGIDList.clear ();
      for (; gid_index < extSize; ++gid_index )
      {
        if ( testIndexMap.find(GIDs[gid_index]) == testIndexMap.end() )
        {
          extGIDList.push_back( IndexPair( GIDs[gid_index], 1 ) );
          testIndexMap[GIDs[gid_index]] = 1;
        }
      }

      testIndexMap.clear ();

      // Copy out internal gids:
      intGIDList.clear ();
      for (; gid_index < intSize + extSize; ++gid_index)
      {
        if ( testIndexMap.find(GIDs[gid_index]) == testIndexMap.end() )
        {
          intGIDList.push_back( IndexPair( GIDs[gid_index], 1 ) );
          testIndexMap[GIDs[gid_index]] = 1;
        }
      }

      testIndexMap.clear ();

      // Copy out the exp var gid's, if they exist.  These will
      // only exist in devices which depend on expressions, like the Bsrc.
      expVarGIDs.clear ();
      for (; gid_index < intSize + extSize + expSize; ++gid_index)
      {
        if ( testIndexMap.find(GIDs[gid_index]) == testIndexMap.end() )
        {
          expVarGIDs.push_back( GIDs[gid_index] );
          testIndexMap[GIDs[gid_index]] = 1;
        }
      }
    }

    // Now copy the exp var GIDs into the extVarGID's.
    // add the contents of expVarGIDs to the extGIDListRef.
    // This is done because the numerical jacobian treats expression
    // GIDs the same as external (nodal) GIDs.
    int expS = expVarGIDs.size();
    for (int i = 0; i < expS; ++i)
    {
      extGIDList.push_back( IndexPair( expVarGIDs[i], 1 ) );
    }

    testIndexMap.clear();

    // do the index pairs for the jacobian matrix
    indexPairList.clear ();
    for (int i = 0; i < jacGIDs.size () ; ++i )
    {
      if ( testIndexMap.find(GIDs[i]) == testIndexMap.end() )
      {
        testIndexMap[GIDs[i]] = 1;

        std::map<int,int> testJMap;
        testJMap.clear ();
        int length = jacGIDs[i].size();
        for( int j = 0; j < length; ++j )
        {
          if ( testJMap.find(jacGIDs[i][j]) == testJMap.end () )
          {
            indexPairList.push_back( IndexPair( GIDs[i], jacGIDs[i][j] ) );
            testJMap[jacGIDs[i][j]] = 1;
          }
        }
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " Complete GIDs :" << std::endl;
      for (std::vector<int>::const_iterator it = GIDs.begin(), end = GIDs.end(); it != end; ++it)
        Xyce::dout() << "\tgid=" << *it << std::endl;

      Xyce::dout() << std::endl;

      Xyce::dout() << " intGIDList :" << std::endl;
      for (IndexPairVector::const_iterator it = intGIDList.begin(), end  = intGIDList.end(); it != end; ++it)
        Xyce::dout() << "\tgid=" << (*it).row << std::endl;

      Xyce::dout() << std::endl;

      Xyce::dout() << " extGIDList :" << std::endl;
      for (IndexPairVector::const_iterator it = extGIDList.begin(), end  = extGIDList.end(); it != end; ++it)
        Xyce::dout() << "\tgid=" << (*it).row << std::endl;

      Xyce::dout() << std::endl;


      Xyce::dout() << " indexPairList :" << std::endl;
      for (IndexPairVector::const_iterator it = indexPairList.begin(), end  = indexPairList.end(); it != end; ++it)
        Xyce::dout() << "  row=" << (*it).row << "  col=" << (*it).col << std::endl;

      Xyce::dout() << std::endl << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
bool DeviceInstance::updateTemperature(const double & temp_tmp)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool DeviceInstance::processParams ()
{
  Report::DevelFatal0().in("DeviceInstance::processParams")
    << "DeviceInstance::processParams() must be implemented for device " << getName();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setInternalState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
bool DeviceInstance::setInternalState(
  const DeviceState & state )
{
  Report::DevelFatal().in("DeviceInstance::setInternalState") << "does not exist for this device " << getName();

  return false;
}

std::ostream &
DeviceInstance::printName(std::ostream &os) const
{
  return os << "instance " << name_.getEncodedName();
}

} // namespace Device
} // namespace Xyce
