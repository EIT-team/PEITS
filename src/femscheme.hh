// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef ELLIPT_FEMSCHEME_HH
#define ELLIPT_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
//#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/lagrange.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptationmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>

// lagrange interpolation 
#include <dune/fem/space/lagrange/interpolation.hh>

#include <dune/fem/misc/femtimer.hh>
/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#endif
#if HAVE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif

// local includes
#include "probleminterface.hh" 

#include "model.hh"

#include "rhs.hh"
#include "electrode_helpers.hh"
#include "elliptic.hh"

#include <time.h>

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

// FemScheme 
//----------

/*******************************************************************************
 * template arguments are:
 * - GridPsrt: the part of the grid used to tesselate the 
 *             computational domain
 * - Model: description of the data functions and methods required for the
 *          elliptic operator (massFlux, diffusionFlux)
 *     Model::ProblemType boundary data, exact solution, 
 *                        and the type of the function space
 *******************************************************************************/

enum SolverType
{
  fem,         // use the matrix based version of the dune-fem solvers
  femoem,      // use the matrix based version of the dune-fem solvers with blas
  istl,        // use the dune-istl solvers
  umfpack,     // use the direct solver umfpack
  petsc        // use the petsc package
};

template <class DFSpace, SolverType solver>
struct Solvers;
template <class DFSpace>
struct Solvers<DFSpace,fem>
{
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
};
template <class DFSpace>
struct Solvers<DFSpace,femoem>
{
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::OEMCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
};
#if HAVE_DUNE_ISTL
template <class DFSpace>
struct Solvers<DFSpace,istl>
{
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
};
#endif // HAVE_ISTL
#if HAVE_UMFPACK
template <class DFSpace>
struct Solvers<DFSpace,umfpack>
{
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType, true > LinearInverseOperatorType;
};
#endif
#if HAVE_PETSC
template <class DFSpace>
struct Solvers<DFSpace,petsc>
{
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
};
#endif
template < class Model, SolverType solver >
class FemScheme 
{
public:   
  //! type of the mathematical model 
  typedef Model ModelType ;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, f: \Omega -> R) 
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

  typedef Solvers<DiscreteFunctionSpaceType,solver> SolverType;

  // choose type of discrete function, Matrix implementation and solver implementation 
  typedef typename SolverType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename SolverType::LinearOperatorType LinearOperatorType;
  typedef typename SolverType::LinearInverseOperatorType LinearInverseOperatorType;
  /*********************************************************/

  //! define Laplace operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > EllipticOperatorType;

  // define entities
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::EntitySeed     EntitySeed;

  FemScheme( GridPartType &gridPart, 
             const ModelType& implicitModel, ElectrodePositions &electrodes, CurrentProtocol &current_protocol )
    : implicitModel_( implicitModel ),
      electrodes_(electrodes),
      current_protocol_(current_protocol),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_( implicitModel_, discreteSpace_, electrodes ), 

/*********************************************************/
/***                 NEW FOR LESSON 2                  ***/
/*********************************************************/
      // create linear operator (domainSpace,rangeSpace)
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_),// implicitOperator_.stencil() ),
/*********************************************************/

      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ), //second value is default
      solverIter_( Dune::Fem::Parameter::getValue< unsigned int >( "poisson.solveriter", std::numeric_limits< unsigned int >::max() ) ), //second value is default
	  invOp_( 0 )
  {
    // set all DoF to zero 
    solution_.clear();

    // get the time to put into name of save files later on
    time_t now;
    time(&now);
    strftime(time_buf_, 21, "%Y-%m-%d_%H:%M:%S", gmtime(&now));
  }
  ~FemScheme()
  {
    delete invOp_;
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  void prepare(int pattern_no)
  {
    solution_.clear();
    implicitOperator_.prepare( implicitModel_.dirichletBoundary(), solution_ );

    const std::vector< std::vector<EntitySeed> > entitySeeds = implicitOperator_.entitySeedVector();
    const std::vector< double > electrodeAreas = implicitOperator_.electrodeAreas();

    // assemble rhs
    assembleRHS ( implicitModel_.rightHandSide(), rhs_, pattern_no, electrodes_, current_protocol_, entitySeeds, electrodeAreas );
    /*{ // test that l(1) = 0
      typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType ; 
      const DofIteratorType dend = rhs_.dend();
      double sum = 0;
      for( DofIteratorType it = rhs_.dbegin(); it != dend; ++it )
         sum += *it;
      MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      assert(std::abs(sum)<1e-12);
    }*/
   
    // apply constraints, e.g. Dirichlet contraints, to the result
    implicitOperator_.prepare( solution_, rhs_ );
  }

  void assemble ( ) 
  {
    // assemble linear operator (i.e. setup matrix)
    implicitOperator_.jacobian( solution_ , linearOperator_ );

    // inverse operator using linear operator 
    if (!invOp_)
      invOp_ = new LinearInverseOperatorType( linearOperator_, solverEps_, solverEps_, solverIter_ );
  }

  //! solve the system 
  void solve ( int protocol_no )
  { 
    // solve system 
    (*invOp_)( rhs_, solution_ );
    if (Dune::Fem::Parameter::getValue< bool >( "fem.solver.shift_average" ))
      shiftAverage(solution_);
  }

  void calculateElecPot(const DiscreteFunctionType &result, const int pattern_index, vector<double> &electrode_potentials)
  {
    const std::vector< std::vector<EntitySeed> > entitySeeds = implicitOperator_.entitySeedVector();
    const std::vector< double > electrodeAreas = implicitOperator_.electrodeAreas();
    ::calculateElecPot( result, pattern_index, electrodes_, current_protocol_, entitySeeds, electrodeAreas, electrode_potentials );
  }

  void writeElecPot( const vector<double> &electrode_potentials, const int measure_index, const double drive_factor, const double measure_factor, const int protocol_step)
  {
    const std::vector< double > electrodeAreas = implicitOperator_.electrodeAreas();
    ::writeElecPot( electrode_potentials, measure_index, time_buf_, current_protocol_, drive_factor, measure_factor, electrodeAreas, protocol_step );
  }

  const std::vector< std::vector<EntitySeed> > electrodeElements() const
  {
    return implicitOperator_.entitySeedVector(); 
  }

  const std::vector< double > electrodeAreas() const
  {
    return implicitOperator_.electrodeAreas();
  }

  const LinearOperatorType &matrix() const
  {
    return linearOperator_;
  }
  const DiscreteFunctionType &rhs() const
  {
    return rhs_;
  }
protected:  
  const ModelType& implicitModel_;   // the mathematical model 

  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with 

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space 
  DiscreteFunctionType solution_;   // the unknown 
  DiscreteFunctionType rhs_;        // the right hand side 

  EllipticOperatorType implicitOperator_; // the implicit operator 

/*********************************************************/
/***                 NEW FOR LESSON 2                  ***/
/*********************************************************/
  LinearOperatorType linearOperator_;  // the linear operator (i.e. jacobian of the implicit)

  //LinearInverseOperatorType invOp_;
/*********************************************************/

  const double solverEps_ ; // eps for linear solver
  const int solverIter_;
  char time_buf_[21];
  LinearInverseOperatorType *invOp_;
  ElectrodePositions electrodes_;
  CurrentProtocol current_protocol_;
};

#endif // end #if ELLIPT_FEMSCHEME_HH
