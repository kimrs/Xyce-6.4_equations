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
// Filename       : $RCSfile: N_LAS_AmesosSolver.C,v $
//
// Purpose        : Amesos direct solver wrapper
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.73 $
//
// Revision Date  : $Date: 2015/08/03 21:12:39 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <Amesos.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_AmesosSolver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_TransformTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::AmesosSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
AmesosSolver::AmesosSolver(
  const std::string &   type,
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(false),
    type_(type),
    lasProblem_(problem),
    problem_(problem.epetraObj()),
    solver_(0),
    repivot_(true),
    reindex_(false),
    outputLS_(0),
    outputBaseLS_(0),
    outputFailedLS_(0),
    tProblem_(0),
    optProb_(0),
    optMat_(0),
    origMat_(0),
    optExporter_(0),
    options_( new Util::OptionBlock( options ) ),
    timer_( new Util::Timer() )
{
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::~AmesosSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
AmesosSolver::~AmesosSolver()
{
  delete solver_;
  delete timer_;
  delete options_;
  delete optProb_;
  delete optMat_;
  delete optExporter_;
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool AmesosSolver::setOptions( const Util::OptionBlock & OB )
{
  for( Util::ParamList::const_iterator it_tpL = OB.begin();
         it_tpL != OB.end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();

    if( tag == "KLU_REPIVOT" ) repivot_ = static_cast<bool>(it_tpL->getImmutableValue<int>());
    
    if( tag == "KLU_REINDEX" ) reindex_ = static_cast<bool>(it_tpL->getImmutableValue<int>());

    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

#ifdef Xyce_PARALLEL_MPI
  options_->addParam(Util::Param("TR_reindex", 1));

  // Turn off partitioning and AMD if we're doing a parallel load serial solve
  if (type_ == "KLU") {
    options_->addParam(Util::Param("TR_partition", 0));
    options_->addParam(Util::Param("TR_amd", 0));
  }
#endif

  if( !transform_.get() ) transform_ = TransformTool()( *options_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool AmesosSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::setDefaultOption
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool AmesosSolver::setDefaultOption( const std::string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool AmesosSolver::setParam( const Util::Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool AmesosSolver::getInfo( Util::Param & info )
{
  return true;
}

//void print(Epetra_MultiVector * vec, char * filename) {
//    FILE * file = fopen(filename, "w");
//    int n = vec->GlobalLength();
//    double values[n];
//    vec->ExtractCopy(&values[0], 0);
//    for(int i = 0; i < n; ++i) {
//        fprintf(file, "%d,%d,%0.16f\n", i, 0, values[i]);
//    }
//    fclose(file);
//}

//-----------------------------------------------------------------------------
// Function      : AmesosSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
//
/* EDIT KIM */
void print(const Epetra_MultiVector *vec, const char * filename) {
    FILE * file = fopen(filename, "w");
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        fprintf(file, "%d,%d,%0.32e\n", 0, i, values[i]);
    
    fclose(file);
}

void print(const Epetra_RowMatrix * A) {
   int n_rows = A->NumMyRows();
   int n_maxCols = A->NumGlobalCols();
   for(int i = 0; i < n_rows; ++i) {
       int n_entries;
       double * values = new double[n_maxCols];
       int * indices = new int[n_maxCols];
       A->ExtractMyRowCopy(i, n_rows, n_entries, values, indices); 
       for(int j = 0; j < n_entries; ++j) 
         printf("%d,%d,%0.32e\n", i, indices[j], values[j]); 

       delete values;
       delete indices;
   }
}
void print(Epetra_RowMatrix * A, std::string syntax, int disambiguation) {
   int n_rows = A->NumMyRows();
   int n_maxCols = A->NumGlobalCols();
   for(int i = 0; i < n_rows; ++i) {
       int n_entries;
       double * values = new double[n_maxCols];
       int * indices = new int[n_maxCols];
       A->ExtractMyRowCopy(i, n_rows, n_entries, values, indices); 
       for(int j = 0; j < n_entries; ++j) 
         printf(syntax.c_str(), i, indices[j], values[j]); 

       delete values;
       delete indices;
   }
}

void print(Epetra_MultiVector *vec, std::string output, int disambiguation) {
    int n = vec->GlobalLength();
    double values[n];
    vec->ExtractCopy(&values[0], 0);
    for(int i = 0; i < n; ++i) 
        printf(output.c_str(), i, 0, values[i]);
}

void print(const Epetra_RowMatrix * A, const char * filename) {
   int n_rows = A->NumMyRows();
   int n_maxCols = A->NumGlobalCols();
   FILE * file = fopen(filename, "w");
   for(int i = 0; i < n_rows; ++i) {
       int n_entries;
       double * values = new double[n_maxCols];
       int * indices = new int[n_maxCols];
       A->ExtractMyRowCopy(i, n_rows, n_entries, values, indices); 
       for(int j = 0; j < n_entries; ++j) 
         fprintf(file, "%d,%d,%0.32e\n", i, indices[j], values[j]); 

       delete values;
       delete indices;
   }
   fclose(file);
}

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

int AmesosSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();

  int linearStatus = 0;

  Epetra_LinearProblem * prob = &problem_;

  if( transform_.get() )
  {
    if( !tProblem_ )
      tProblem_ = &((*transform_)( problem_ ));
    prob = tProblem_;
    transform_->fwd();
  }

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int failure_number = 0, file_number = 1, base_file_number = 1;
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      if (!reuse_factors) {
        if (file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Transformed_BlockMap.mm", (prob->GetMatrix())->Map() );
        }
        sprintf( file_name, "Transformed_Matrix%d.mm", file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Transformed_RHS%d.mm", file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
      }
    }
    // file_number++;  This will be incremented after the solution vector is written to file.
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (!reuse_factors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (problem_.GetMatrix())->Map() );
        }
        sprintf( file_name, "Base_Matrix%d.mm", base_file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_.GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Base_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetRHS()) );
      }
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

  // Set the traceback mode in Epetra so it prints out warnings
  if (DEBUG_LINEAR)
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 2 );
  else
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 0 );


  Amesos localAmesosObject;
  if( !solver_ )
  {
    // the Query() function expects a string
    // in lower case with the first letter in upper case
    // So, our "KLU" must become "Klu"
    std::string solverType( type_ );
    if( type_ == "KLU" )
    {
      solverType = "Amesos_Klu";
    }
    else if( type_ == "SUPERLU" )
    {
      solverType = "Amesos_Superlu";
    }
    else if( type_ == "SUPERLUDIST" )
    {
      solverType = "Amesos_Superludist";
    }
    else if( type_ == "PARAKLETE" )
    {
      solverType = "Amesos_Paraklete";
    }
    else if( type_ == "PARDISO" )
    {
      solverType = "Amesos_Pardiso";
    }
    else if( type_ == "LAPACK" )
    {
      solverType = "Amesos_Lapack";
    }
    else if( type_ == "SCALAPACK" )
    {
      solverType = "Amesos_Scalapack";
    }
    else if( type_ == "MUMPS" )
    {
      solverType = "Amesos_Mumps";
    }

    if( !localAmesosObject.Query( solverType ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
                              "Unknown or Unavailable Linear Solver: " + type_ );


#ifndef Xyce_PARALLEL_MPI
    //setup optimized storage version of problem for serial
    //only do this if the linear system is nontrivial (not a single equation)
    origMat_ = dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix());
    if (origMat_->NumGlobalRows() > 1) {
      Epetra_Map const& rowMap = origMat_->RowMap();
      Epetra_BlockMap const& blockRowMap = dynamic_cast<Epetra_BlockMap const&>(rowMap);
      optMat_ = new Epetra_CrsMatrix( Copy, rowMap, 0 );
      optExporter_ = new Epetra_Export( blockRowMap, blockRowMap );
      optMat_->Export( *origMat_, *optExporter_, Insert );
      optMat_->FillComplete();
      optMat_->OptimizeStorage();

      optProb_ = new Epetra_LinearProblem( optMat_, prob->GetLHS(), prob->GetRHS() );
      prob = optProb_;
    }
#endif

    solver_ = localAmesosObject.Create( solverType, *prob );

    Teuchos::ParameterList params;

#ifndef Xyce_PARALLEL_MPI
    // Inform solver not to check inputs to reduce overhead.
    params.set( "TrustMe", true );
    // If repivot == true (default), recompute the pivot order each numeric factorization,
    // else try to re-use pivot order to expedite numeric factorization.
    params.set( "Refactorize", !repivot_ );
#else
    if (type_ == "SUPERLUDIST") {
      Teuchos::ParameterList& sludistParams = params.sublist("Superludist");
      sludistParams.set("ReuseSymbolic", true );
    }
#endif

    // Let Amesos reindex the linear problem.
    // NOTE:  This is used by MPDE and HB since the map indices are not continguous.
    if (reindex_) {
      params.set( "Reindex", reindex_ );
    }

    if (VERBOSE_LINEAR)
      Xyce::dout() << "AmesosSolver::solve() setting solver : " << type_ << "\n"
                   << "AmesosSolver::solve() setting parameters : " << params << std::endl;

    solver_->SetParameters( params );

    double begSymTime = timer_->elapsedTime();

    // Perform symbolic factorization and check return value for failure
    linearStatus = solver_->SymbolicFactorization();
    if (linearStatus != 0)
      return linearStatus;

    if (VERBOSE_LINEAR)
    {
      double endSymTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos (" << type_ << ") Symbolic Factorization Time: "
                   << (endSymTime - begSymTime) << std::endl;
    }
  }

  if( optMat_ ) optMat_->Export( *origMat_, *optExporter_, Insert );

  // Set the transpose flag only if that has changed since the last solve.
  if ( solver_->UseTranspose() != transpose )
  {
    solver_->SetUseTranspose( transpose );
  }

  // Perform numeric factorization and check return value for failure
  if( !reuse_factors ) {

    double begNumTime = timer_->elapsedTime();
    //Amesos_BaseSolver bs = *solver_; /* EDIT KIM */
    solver_->PrintStatus();
    linearStatus = solver_->NumericFactorization();
    std::cout << "linearStatus: " << linearStatus << "\n";

    if (VERBOSE_LINEAR)
    {
      double endNumTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos (" << type_ << ") Numeric Factorization Time: "
                   << (endNumTime - begNumTime) << std::endl;
    }
    
    if (linearStatus != 0) {

      // Inform user that singular matrix was found and linear solve has failed.
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,
                              "Numerically singular matrix found by Amesos, returning zero solution to nonlinear solver!");

      // Put zeros in the solution since Amesos was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );
      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        char file_name[40];
        if (failure_number== 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Failed_BlockMap.mm", (prob->GetMatrix())->Map() );
        }
        sprintf( file_name, "Failed_Matrix%d.mm", failure_number );
        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Failed_RHS%d.mm", failure_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
      }

      return linearStatus;  // return the actual status (see bug 414 SON)
    }
  }

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  solver_->Solve();
  /* EDIT KIM */
  print(prob->GetMatrix(),  genFileName("Jac").c_str());
  print(prob->GetRHS(),     genFileName("rhs").c_str());
  print(prob->GetLHS(),     genFileName("lhs").c_str());

  std::cout << "jac\n";
  print(prob->GetMatrix()   ,"%d,%d,%0.64e\n", 0);
  std::cout << "rhs\n";
  print(prob->GetRHS()      ,"%d,%d,%0.64e\n", 0);
  std::cout << "lhs\n";
  print(prob->GetLHS()      ,"%d,%d,%0.64e\n", 0);
  std::cout << "done\n";
  /* END KIM */

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::dout() << "  Amesos (" << type_ << ") Solve Time: "
                 << (endSolveTime - begSolveTime) << std::endl;
  }

  if (DEBUG_LINEAR) {
    double resNorm = 0.0, bNorm = 0.0;
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm );
    prob->GetRHS()->Norm2( &bNorm );
    Xyce::lout() << "Linear System Residual (AMESOS_" << type_ << "): "
                 << (resNorm/bNorm) << std::endl;
  }

  if( transform_.get() ) transform_->rvs();

  // Output computed solution vectors, if requested.
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetLHS()) );
    }
    file_number++;
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      sprintf( file_name, "Base_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetLHS()) );
    }
    base_file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::dout() << "Total Linear Solution Time (Amesos " << type_ << "): "
                 << solutionTime_ << std::endl;

  return 0;
}

} // namespace Linear
} // namespace Xyce
