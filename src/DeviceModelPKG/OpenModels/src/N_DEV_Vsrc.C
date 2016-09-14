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
// Filename       : $RCSfile: N_DEV_Vsrc.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.250 $
//
// Revision Date  : $Date: 2015/04/20 20:43:55 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

#include <Epetra_CrsMatrix.h>

namespace Xyce {
namespace Device {

namespace Vsrc {

void Traits::loadInstanceParameters(ParametricData<Vsrc::Instance> &p)
{
    p.addPar ("DCV0", 0.0, &Vsrc::Instance::DCV0)
      .setOriginalValueStored(true)
      .setUnit(U_VOLT)
      .setDescription("DC Voltage")
      .setAnalyticSensitivityAvailable(true)
      .setSensitivityFunctor(&dcv0Sens);

    // Pulse parameters
    p.addPar ("V0", 0.0, &Vsrc::Instance::par0)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Offset Voltage");

    p.addPar ("V1", 0.0, &Vsrc::Instance::par0)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Initial Voltage");

    p.addPar ("V2", 0.0, &Vsrc::Instance::par1)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Pulsed Voltage");

    p.addPar ("TD", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Delay");

    p.addPar ("TR", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Time");

    p.addPar ("TF", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Time");

    p.addPar ("PW", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Pulse Width");

    p.addPar ("PER", 0.0, &Vsrc::Instance::par6)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Period");

    p.addPar ("SF", 0.0, &Vsrc::Instance::par7)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Scale Factor -- smooth pulse only");

    // Sin parameters
    p.addPar ("VA", 0.0, &Vsrc::Instance::par1)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Amplitude");

    p.addPar ("FREQ", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Frequency");

    p.addPar ("THETA", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Theta");

    p.addPar ("PHASE", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Phase");

    // Exp parameters
    p.addPar ("TD1", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Delay Time");

    p.addPar ("TAU1", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Time Constant");

    p.addPar ("TD2", 0.0, &Vsrc::Instance::par4)
      .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Delay Time");

    p.addPar ("TAU2", 0.0, &Vsrc::Instance::par5)
      .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Time Constant");

    // AC parameters
    p.addPar ("ACMAG", 0.0, &Vsrc::Instance::ACMAG)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Amplitude");

    p.addPar ("ACPHASE", 0.0, &Vsrc::Instance::ACPHASE)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Phase");

    // SFFM parameters
    p.addPar ("FC", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Carrier Frequency");

    p.addPar ("FS", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Signal Frequency");

    p.addPar ("MDI", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Modulation Index");

    // PWL params
    p.addPar ("R", 0.0, &Vsrc::Instance::REPEATTIME)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Repeat Time");

    p.addPar ("T", 0.0, &Vsrc::Instance::T)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Time");  // time-voltage pairs

    p.addPar ("V", 0.0, &Vsrc::Instance::V)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Voltage"); // time-voltage pairs

    // Set up exceptions (ie variables that are not doubles):
    p.addPar ("TRANSIENTSOURCETYPE", (int) _DC_DATA, &Vsrc::Instance::TRANSIENTSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::TRANSIENTSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("ACSOURCETYPE", (int) _AC_DATA, &Vsrc::Instance::ACSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::ACSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("DCSOURCETYPE", (int) _DC_DATA, &Vsrc::Instance::DCSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::DCSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("NUM", 0, &Vsrc::Instance::NUM)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );
}

void Traits::loadModelParameters(ParametricData<Vsrc::Model> &p)
{
}


std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : SourceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    srcCurrent(0.0),
    srcVoltage(0.0),
    srcDrop(0.0),
    srcBC(0.0),

    scale(1.0),
    nlstep(-1),
    ACSpecified_(factory_block.deviceManager_.getACSpecified()),
    HBSpecified_(factory_block.deviceManager_.getHBSpecified()),
    DCV0(0.0),
    par0(0.0),
    par1(0.0),
    par2(0.0),
    par3(0.0),
    par4(0.0),
    par5(0.0),
    par6(0.0),
    par7(0.0),
    REPEATTIME(),
    T(0.0),
    V(0.0),
    ACMAG(1.0),
    ACPHASE(0.0),
    NUM(0),
    REPEAT(false),
    TRANSIENTSOURCETYPE(_DC_DATA),
    TRANSIENTSOURCETYPEgiven(false),
    ACSOURCETYPE(_AC_DATA),
    ACSOURCETYPEgiven(false),
    DCSOURCETYPE(_AC_DATA),
    DCSOURCETYPEgiven(false),

    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),

    gotParams(false),

    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    ABraEquBraVarOffset(-1),
    APosEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0),
    fPosEquPosNodePtr(0),
    fNegEquNegNodePtr(0),
    fBraEquBraVarPtr(0),
#endif

    source(0.0),
    v_pos(0.0),
    v_neg(0.0),
    i_bra(0.0),
    li_branch_data(0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(3);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0] = 0;
    jacStampPDE[0][1] = 2;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0] = 1;
    jacStampPDE[1][1] = 2;
    jacStampPDE[2].resize(3);
    jacStampPDE[2][0] = 0;
    jacStampPDE[2][1] = 1;
    jacStampPDE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  const SolverState &solver_state = factory_block.solverState_;
  const DeviceOptions &device_options = factory_block.deviceOptions_;

  // Set any non-constant parameter defaults:
  if (ACSpecified_ && ACSOURCETYPEgiven)
  {
    acSourceData_ = new ACData(*this, IB.params, solver_state, device_options);
  }

  if (DCSOURCETYPEgiven) // this will always be given, if the source spec was valid.
  {
    dcSourceData_ = new ConstData(*this, IB.params, solver_state, device_options);
  }

  if (HBSpecified_ || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        tranSourceData_ = new SinData(*this, IB.params, solver_state, device_options);
        break;

      case _EXP_DATA:
        tranSourceData_ = new ExpData(*this, IB.params, solver_state, device_options);
        break;

      case _PULSE_DATA:
        tranSourceData_ = new PulseData(*this, IB.params, solver_state, device_options);
        break;

      case _PWL_DATA:
        tranSourceData_ = new PWLinData(*this, IB.params, solver_state, device_options);
        break;

      case _SFFM_DATA:
        tranSourceData_ = new SFFMData(*this, IB.params, solver_state, device_options);
        break;

      case _DC_DATA:
        tranSourceData_ = 0; // this forces  us to use the dcSourceData_ object instead
        break;

      case _SMOOTH_PULSE_DATA:
        tranSourceData_ = new SmoothPulseData(*this, IB.params, solver_state, device_options);
        break;

      default:
        UserFatal0(*this) << "Cannot identify source data type for " << getName();
        break;
    }
  }

  processParams();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();
  processParams();

  // calculate dependent (ie computed) params and check for errors:

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/26/03
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  if (gotParams)
  {
    if (dcSourceData_ != 0)
    {
      dcSourceData_->setParams (&DCV0);
    }
    if (acSourceData_ != 0)
    {
      acSourceData_->setParams (&ACMAG);
    }
    if (tranSourceData_ != 0)
    {
      tranSourceData_->setParams(&par0);
    }
  }
  else
  {
    if (dcSourceData_ != 0)
    {
      dcSourceData_->getParams (&DCV0);
    }
    if (acSourceData_ != 0)
    {
      acSourceData_->getParams (&ACMAG);
    }
    if (tranSourceData_ != 0)
    {
      tranSourceData_->getParams(&par0);
    }
    gotParams = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
  delete tranSourceData_;
  delete acSourceData_;
  delete dcSourceData_;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  VsrcInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal0(*this).in("Instance::registerLIDs") << "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal0(*this).in("Instance::registerLIDs") << "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The Voltage Source device uses exactly one "branch data vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(V1)" for the current through resistor V1. and a junction voltage.
///
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_Bra, getName(), "branch");
  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().isPDESystem_)
    return jacStampPDE;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (getSolverState().isPDESystem_)
  {
    APosEquBraVarOffset  = jacLIDVec[0][1];
    ANegEquBraVarOffset  = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }
  else
  {
    APosEquBraVarOffset  = jacLIDVec[0][0];
    ANegEquBraVarOffset  = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // Get the value for the source.
  SourceData *dataPtr = dcSourceData_; // by default assume the DC value.
  if ((HBSpecified_ || getSolverState().tranopFlag || getSolverState().transientFlag) && tranSourceData_ != 0 )
  {
    dataPtr = tranSourceData_;
  }

  if (dataPtr != 0)
  {
    source = dataPtr->returnSource();
  }
  else
  {
    source = 0.0;
  }

  // get the value for v_pos, v_neg, i_bra
  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];
  i_bra = solVec[li_Bra];

  srcCurrent = i_bra;
  srcDrop    = (v_pos-v_neg);
  srcBC      = source;
  srcVoltage = srcDrop-srcBC;

  if( getDeviceOptions().scale_src != 0.0 )
  {
    srcCurrent *= scale;

    // first newton step, generate new scaling
    //if( nlsMgrPtr->getNonLinearIter() != nlstep )
    if( getSolverState().newtonIter != nlstep )
    {
      //nlstep = nlsMgrPtr->getNonLinearIter();
      nlstep = getSolverState().newtonIter;
      double new_scale = fabs(i_bra) * scale;

      scale = std::max( new_scale, getDeviceOptions().scale_src );
    }

    srcVoltage *= scale;
    srcDrop    *= scale;
    srcBC      *= scale;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       : Loads the F-vector contributions for a single vsrc instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;

  fVec[li_Pos] += srcCurrent;
  fVec[li_Neg] += -srcCurrent;
  fVec[li_Bra] += srcDrop;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data] = srcCurrent;
    junctionV[li_branch_data] = srcDrop;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  double * bVec = extData.daeBVectorRawPtr;
  bVec[li_Bra] += srcBC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBVectorsforAC
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool Instance::loadBVectorsforAC(double * bVecReal, double * bVecImag )
{
  if (acSourceData_ != 0)
  {
    bool flag = true;
    acSourceData_->setRealFlag(flag);

    acSourceData_->updateSource ();
    source = acSourceData_->returnSource();
    srcBC = source;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      srcBC *= scale;
    }

    bVecReal[li_Bra] += srcBC;

    flag = false;
    acSourceData_->setRealFlag(flag);

    acSourceData_->updateSource ();
    source = acSourceData_->returnSource();
    srcBC = source;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      srcBC *= scale;
    }

    bVecImag[li_Bra] += srcBC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;
  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize  ()
{
  double maxStep = 1.0e+100;
  if (tranSourceData_ != 0)
  {
    maxStep = tranSourceData_->getMaxTimeStepSize ();
  }
  return maxStep;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    DC_TRAN (0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
///
/// @param[in] op Operator to apply to all instances.
///
///
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}



//-----------------------------------------------------------------------------
// Vsrc Master functions:
//-----------------------------------------------------------------------------

Master::Master(
  const Configuration & configuration,
  const FactoryBlock &  factory_block,
  const SolverState &   solver_state,
  const DeviceOptions & device_options)
  : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options),
    HBSpecified_(factory_block.deviceManager_.getHBSpecified())
  {}

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi        = *(*it);
    // Get the value for the source.
    SourceData *dataPtr  = vi.dcSourceData_; // by default assume the DC value.
    if ((HBSpecified_ || getSolverState().tranopFlag || getSolverState().transientFlag) && vi.tranSourceData_ != 0 )
    {
      dataPtr            = vi.tranSourceData_;
    }

    if (dataPtr != 0)
    {
      vi.source          = dataPtr->returnSource();
    }
    else
    {
      vi.source          = 0.0;
    }

    // get the value for v_pos, v_neg, i_bra
    vi.v_pos             = solVec[vi.li_Pos];
    vi.v_neg             = solVec[vi.li_Neg];
    vi.i_bra             = solVec[vi.li_Bra];

    vi.srcCurrent        = vi.i_bra;
    vi.srcDrop           = (vi.v_pos-vi.v_neg);
    vi.srcBC             = vi.source;
    vi.srcVoltage        = vi.srcDrop-vi.srcBC;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      vi.srcCurrent     *= vi.scale;

      // first newton step, generate new scaling
      if( getSolverState().newtonIter != vi.nlstep )
      {
        vi.nlstep        = getSolverState().newtonIter;
        double new_scale = fabs(vi.i_bra) * vi.scale;

        vi.scale         = std::max( new_scale, getDeviceOptions().scale_src );
      }

      vi.srcVoltage     *= vi.scale;
      vi.srcDrop        *= vi.scale;
      vi.srcBC          *= vi.scale;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * storeLeadF, double * storeLeadQ, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    (*it)->F_Strings = F_strings;
    (*it)->Q_Strings = Q_strings;
    (*it)->B_Strings = B_strings;

    /* EDIT KIM */
    Epetra_MultiVector old_epetra_f(externData.daeFVectorPtr->epetraObj());
    Linear::MultiVector old_F(&old_epetra_f, false); /* SIGABRT if set to true you fucker */

    Epetra_MultiVector old_epetra_b(externData.daeBVectorPtr->epetraObj());
    Linear::MultiVector old_B(&old_epetra_b, false); /* SIGABRT if set to true you fucker */
    /* END KIM */

    Instance & vi = *(*it);
    fVec[vi.li_Pos] += vi.srcCurrent;
    fVec[vi.li_Neg] += -vi.srcCurrent;
    fVec[vi.li_Bra] += vi.srcDrop;

    bVec[vi.li_Bra] += vi.srcBC;

    if( vi.loadLeadCurrent )
    {
      leadF[vi.li_branch_data] = vi.srcCurrent;
      junctionV[vi.li_branch_data] = vi.srcDrop;
    }

    /* EDIT KIM */
    /* Old F is set in the beginning of the for loop */
    Epetra_MultiVector tmp_epetra_f(externData.daeFVectorPtr->epetraObj());
    Linear::MultiVector tmp_F(&tmp_epetra_f, false);

    /* tmp_F = F - old_F */
    tmp_F.linearCombo(1.0, tmp_F, -1.0, old_F);
    
    Epetra_CrsMatrix crs_mat_f = ExternData::multivec2crsMat(tmp_F.epetraObj());
    Linear::Matrix mat_f(&crs_mat_f, false);
    vi.print_strings(mat_f, vi.F_Strings, vi.getName().getEncodedName());

    Epetra_MultiVector tmp_epetra_b(externData.daeBVectorPtr->epetraObj());
    Linear::MultiVector tmp_B(&tmp_epetra_b, false);

    /* tmp_F = F - old_F */
    tmp_B.linearCombo(1.0, tmp_B, -1.0, old_B);
    
    Epetra_CrsMatrix crs_mat_b = ExternData::multivec2crsMat(tmp_B.epetraObj());
    Linear::Matrix mat_b(&crs_mat_b, false);
    vi.print_strings(mat_b, vi.B_Strings, vi.getName().getEncodedName());

    /* END KIM */
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);
    /* EDIT KIM */
    Epetra_CrsMatrix old(dFdx.epetraObj());
    Linear::Matrix old_dFdx(&old, false); /* SIGABRT if set to true you fucker */
    /* END KIM */

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(vi.fPosEquBraVarPtr) += 1.0;
    *(vi.fNegEquBraVarPtr) -= 1.0;
    *(vi.fBraEquPosNodePtr) += 1.0;
    *(vi.fBraEquNegNodePtr) -= 1.0;
#else
    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;
    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;
    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;
    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
#endif
    /* EDIT KIM */
    /* Old dQdx is set in the beginning of the for loop */
    Epetra_CrsMatrix tmp(dFdx.epetraObj());
    Linear::Matrix tmp_dFdx(&tmp, false);

    /* tmp_dQdx = dQdx - old_dQdx */
    tmp_dFdx.linearCombo(1.0, tmp_dFdx, -1.0, old_dFdx);
    vi.print_strings(tmp_dFdx, vi.dFdx_Strings, vi.getName().getEncodedName());
    /* END KIM */

  }
  return true;
}

Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("v", 1);
}


//-----------------------------------------------------------------------------
// Function      : dcVsrcSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=DCV0.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/18/2014
//-----------------------------------------------------------------------------
void dcVsrcSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp,
    std::vector<double> & dqdp,
    std::vector<double> & dbdp,
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  dbdp.resize(1);
  dbdp[0] += 1.0;
  Bindices.resize(1);
  Bindices[0] = in->li_Bra;
}

} // namespace Vsrc
} // namespace Device
} // namespace Xyce
