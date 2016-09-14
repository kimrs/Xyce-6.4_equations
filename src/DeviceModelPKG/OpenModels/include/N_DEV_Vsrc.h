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
// Filename       : $RCSfile: N_DEV_Vsrc.h,v $
//
// Purpose        : Independent voltage source classes.
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
// Revision Number: $Revision: 1.124 $
//
// Revision Date  : $Date: 2015/04/20 20:43:55 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Vsrc_h
#define Xyce_N_DEV_Vsrc_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <Epetra_CrsMatrix.h>


namespace Xyce {
namespace Device {
namespace Vsrc {

class Model;
class Instance;

/// sensitivity functor
class dcVsrcSensitivity :  public baseSensitivity
{
  public:
  dcVsrcSensitivity() : 
    baseSensitivity() {};

  virtual ~dcVsrcSensitivity() {};

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const ;
};


static dcVsrcSensitivity dcv0Sens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Independent Voltage Source";}
  static const char *deviceTypeName() {return "V level 1";}
  static ParametricData<Instance> &instanceParameters();

  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "DCV0";}
  static const char *instanceDefaultParameter() {return "DCV0";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
  static ParametricData<Model> &modelParameters();
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Instance : public SourceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
  friend class Master;
  friend class dcVsrcSensitivity;

public:
/* EDIT KIM */
 // EquationStrings * dFdx_Strings;
 // EquationStrings * dQdx_Strings;
 // void print_strings(Linear::Matrix &matrix, EquationStrings * strings, std::string name); 
 // int m1; int n1; 
 // int m2; int n2; 
 // int m3; int n3; 
 // int m4; int n4; 
/* END KIM */

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Viter,
     const FactoryBlock &        factory_block);

  Instance(const Instance & right);
  ~Instance();

  bool isLinearDevice() const { return true; }

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  // load functions, residual:
  bool loadBVectorsforAC (double * bVecReal, double * bVecImag);
  bool loadDAEQVector () { return true; }
  bool loadDAEFVector ();
  bool loadDAEBVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () { return true; }
  bool loadDAEdFdx ();

  void setupPointers();

  double getMaxTimeStepSize ();

  void varTypes( std::vector<char> & varTypeVec );

  void getLIDs(int & lpos, int & lneg,int & lbra)
  {lpos = li_Pos; lneg = li_Neg; lbra = li_Bra;}

public:
  // iterator reference to the vsrc model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // state variables:
  double srcCurrent;
  double srcVoltage;
  double srcDrop;
  double srcBC;

  // scale factor
  double scale;
  int nlstep;

  bool          HBSpecified_;
  bool          ACSpecified_;

  // Parameters
  double DCV0;
  double par0;
  double par1;
  double par2;
  double par3;
  double par4;
  double par5;
  double par6;
  double par7;
  double REPEATTIME;
  double T;
  double V;
  double ACMAG;
  double ACPHASE;

  int NUM;
  bool REPEAT;
  int TRANSIENTSOURCETYPE;
  bool TRANSIENTSOURCETYPEgiven;
  int ACSOURCETYPE;
  bool ACSOURCETYPEgiven;
  int DCSOURCETYPE;
  bool DCSOURCETYPEgiven;
  bool gotParams;

  // load variables
  double source, v_pos, v_neg, i_bra;

  // indices into state vector:
  int istate_I;  // index for i0;

  // Matrix equation index variables:

  //local indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;
  int li_branch_data;         ///< Index for Lead Current and junction voltage (for power calculations)

  // Jacobian matrix indices:
  //Locally indexed offsets for jacobian
  int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
  // branch current equ.

  int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
  // branch current equ.

  int APosEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the pos node

  int ANegEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the neg node

  //  The following jacobian offsets are only neccessary
  // for 2-level newton.
  int APosEquPosNodeOffset;  // Offset, positive node variable
  // contribution, positive node KCL.

  int ANegEquNegNodeOffset;  // Offset, negative node variable
  // contribution, negative node KCL.

  int ABraEquBraVarOffset;  // Offset, branch current variable
  // contribution, branch current equation.

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Jacobian matrix pointers:
  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fPosEquBraVarPtr;
  double * fNegEquBraVarPtr;

  //  The following jacobian pointers are only neccessary for 2-level newton.
  double * fPosEquPosNodePtr;
  double * fNegEquNegNodePtr;
  double * fBraEquBraVarPtr;
#endif

  static std::vector< std::vector<int> > jacStamp;
  static std::vector< std::vector<int> > jacStampPDE;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;
  friend class Master;
  friend class dcVsrcSensitivity;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
     const FactoryBlock &      factory_block);
  ~Model ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;
  virtual bool processParams() 
  {
    return true;
  }

  virtual bool processInstanceParams() 
  {
    return true;
  }


public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

private:
  std::vector<Instance*> instanceContainer;

private:

  // This is the dc and transient analysis value of the source.
  double DC_TRAN;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  /*EDIT KIM */
  Master(
    const Configuration &       configuration,
    const FactoryBlock &        factory_block,
    const SolverState &         solver_state,
    const DeviceOptions &       device_options);

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * storeLeadF, double * storeLeadQ, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);

private:
  bool          HBSpecified_;
};

void registerDevice();

} // namespace Vsrc
} // namespace Device
} // namespace Xyce

#endif
