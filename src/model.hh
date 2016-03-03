#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

#include "probleminterface.hh"
#include "rhs.hh"

// DiffusionModel
// --------------

template< class FunctionSpace, class GridPart, class SigmaFunctionType >
struct DiffusionModel
{
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef ProblemInterface< FunctionSpaceType > ProblemType ;

protected:
  enum FunctionId { rhs, bnd };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bnd>, GridPartType > DirichletBoundaryType;

  //! constructor taking problem reference 
  DiffusionModel( const ProblemType& problem, const GridPart &gridPart,
		  const SigmaFunctionType &sigma, const SigmaFunctionType &elementID )
    : problem_( problem ),
      gridPart_(gridPart),
      rhs_(problem_),
      bnd_(problem_),
      sigma_(sigma),
  	  elementID_(elementID),
      uniformCond_(Dune::Fem::Parameter::getValue< bool >( "fem.uniform_conductivity" )),
      uniformCondValue_(Dune::Fem::Parameter::getValue< double >( "fem.uniform_conductivity_value" )),
      assignCond_(Dune::Fem::Parameter::getValue< bool >( "fem.assign_conductivities" )),
      havePerturbation_(Dune::Fem::Parameter::getValue< bool >( "mesh.perturbation" )),
      pertRadius_(Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.radius" )),
      multORabs_(Dune::Fem::Parameter::getValue< bool >( "mesh.perturbation.multORabs" )),
      value_(Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.value" )),
      pertcenter_(3)
  {
	  pertcenter_[0]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_x" );
	  pertcenter_[1]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_y" );
	  pertcenter_[2]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_z" );

    // load the conductivity values of the different tissues
    const string filename = Dune::Fem::Parameter::getValue< string >( "conductivities" );
    error_ = 0;
    FILE *F = 0;
	  F = fopen(filename.c_str(), "r");
    if (F == NULL)
    {
      cout << "Conductivities file could not be loaded!" << endl;
      error_ = 1;
      return;
    }

	  while(!feof(F))
	  {
		  double c;
		  int error = fscanf(F, "%lf\n", &c);
		  if (!error)
		    cout << "File could not be read: " << filename << endl;
		  conductivities_.push_back(c);
	  }
	  fclose(F);
  }

  template< class Entity >
  void sigmaValue( const Entity &entity,
                   RangeType &sigmavalue ) const
  {
    if (uniformCond_)
    {
      sigmavalue = uniformCondValue_;
      if (havePerturbation_)
        applyPerturbation(entity.geometry().center(), sigmavalue);
      return;
    }
    
    sigmavalue = sigma_.localFunction(entity)[0];
    
    if (!assignCond_)
    {
      if (havePerturbation_)
      applyPerturbation(entity.geometry().center(), sigmavalue);
      return;
    }

    sigmavalue = conductivities_[sigmavalue-1];

    if (havePerturbation_)
      applyPerturbation(entity.geometry().center(), sigmavalue);
  }

  template< class CenterType >
  void applyPerturbation( const CenterType &center,
                   RangeType &sigmavalue ) const
  {
    CenterType electrodecenter(center);

	    /*electrodecenter[0]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_x" );
	    electrodecenter[1]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_y" );
	    electrodecenter[2]=Dune::Fem::Parameter::getValue< double >( "mesh.perturbation.pos_z" );*/
    
    double radius = pertRadius_/1000;

    //electrodecenter -= center;
    electrodecenter[0] -= pertcenter_[0]; electrodecenter[1] -= pertcenter_[1]; electrodecenter[2] -= pertcenter_[2];
    if( electrodecenter.two_norm() < radius )
    {
      if (multORabs_)
        sigmavalue *= value_;
      else
        sigmavalue = value_;
    } 
  }

  template< class Entity >
  void elementID( const Entity &entity,
                   RangeType &elementID ) const
  {
    elementID = elementID_.localFunction(entity)[0] - 1;
  }


/*********************************************************/

  //! exact some methods from the problem class
  bool hasDirichletBoundary () const 
  {
    return problem_.hasDirichletBoundary() ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  bool isDirichletPoint( const DomainType& x ) const 
  {
    return problem_.isDirichletPoint(x) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar, 
          const Entity &entity, 
          const Point &x,
          RangeType &u ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    problem_.g( xGlobal, u );
  }

  // return Fem :: Function for Dirichlet boundary values 
  DirichletBoundaryType dirichletBoundary( ) const 
  {
    return DirichletBoundaryType( "boundary function", bnd_, gridPart_, 5 );  
  }

  // return Fem :: Function for right hand side 
  RightHandSideType rightHandSide(  ) const 
  {
    return RightHandSideType( "right hand side", rhs_, gridPart_, 5 );  
  }

  int error_;
   
protected:
  template <FunctionId id>
  class FunctionWrapper : public Dune::Fem::Function< FunctionSpaceType, FunctionWrapper< id > >
  {
    const ProblemInterface<FunctionSpaceType>& impl_;
    public:   
    FunctionWrapper( const ProblemInterface<FunctionSpaceType>& impl )
    : impl_( impl ) {}
 
    //! evaluate function 
    void evaluate( const DomainType& x, RangeType& ret ) const 
    {
      if( id == rhs ) 
      {
        // call right hand side of implementation 
        impl_.f( x, ret );
      }
      else if( id == bnd ) 
      {
        // call dirichlet boudary data of implementation 
        impl_.g( x, ret );
      }
      else 
      {
        DUNE_THROW(Dune::NotImplemented,"FunctionId not implemented"); 
      }
    }
  };
   
  const ProblemType& problem_;
  const GridPart &gridPart_;
  FunctionWrapper<rhs> rhs_;
  FunctionWrapper<bnd> bnd_;
  const SigmaFunctionType sigma_;
  const SigmaFunctionType elementID_;
  const bool uniformCond_;
  const double uniformCondValue_;
  const bool assignCond_;
  const bool havePerturbation_;
  const double pertRadius_;
  const bool multORabs_;
  const double value_;
  std::vector<double> pertcenter_;
  std::vector<double> conductivities_;
};

#endif // #ifndef ELLIPTC_MODEL_HH
