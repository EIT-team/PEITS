#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include "probleminterface.hh"

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

  static const int dimRange = FunctionSpaceType::dimRange;

  static const bool isLinear = true;
  static const bool isSymmetric = true;

protected:
  enum FunctionId { rhs, bndD, bndN };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndD>, GridPartType > DirichletBoundaryType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndN>, GridPartType > NeumanBoundaryType;

  //! constructor taking problem reference
  DiffusionModel( const ProblemType& problem, const GridPart &gridPart,
		  const SigmaFunctionType &sigma, const SigmaFunctionType &elementID )
    : problem_( problem ),
      gridPart_(gridPart),
      rhs_(problem_),
      bndD_(problem_),
      bndN_(problem_),
      penalty_(Dune::Fem::Parameter::getValue<double>("dg.penalty")),
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

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &flux ) const
  {
    linSource( value, gradient, entity, x, value, gradient, flux );
  }

  // the linearization of the source function
  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const JacobianRangeType &gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   const JacobianRangeType &gradient,
                   RangeType &flux ) const
  {
    const DomainType xGlobal = entity.geometry().global( Dune::Fem::coordinate( x ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<flux.size();++i)
      flux[i] = m[i]*value[i];
  }
  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  // linearization of diffusiveFlux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar,
                          const JacobianRangeType& gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
    // the flux is simply the identity
    flux = gradient;
  }

  template< class Entity, class Point >
  void alpha(const Entity &entity, const Point &x,
             const RangeType &value,
             RangeType &val) const
  {
    linAlpha(value,entity,x,value,val);
  }
  template< class Entity, class Point >
  void linAlpha(const RangeType &uBar,
                const Entity &entity, const Point &x,
                const RangeType &value,
                RangeType &val) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType alpha;
    problem_.alpha(xGlobal,alpha);
    for (unsigned int i=0;i<val.size();++i)
      val[i] = alpha[i]*value[i];
  }

  // extract some methods from the problem class

  //! return true if Dirichlet boundary is set
  bool hasDirichletBoundary () const
  {
    return problem_.hasDirichletBoundary() ;
  }

  //! return true if Neumann boundary is set
  bool hasNeumanBoundary () const
  {
    return problem_.hasNeumanBoundary() ;
  }

  //! return true if given intersection belongs to the Dirichlet boundary -
  //! we test here if the center is a Dirichlet point
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return problem_.isDirichletPoint( inter.geometry().center() ) ;
  }

  // return Fem :: Function for Dirichlet boundary values
  DirichletBoundaryType dirichletBoundary( ) const
  {
    return DirichletBoundaryType( "boundary function", bndD_, gridPart_, 5 );
  }
  NeumanBoundaryType neumanBoundary( ) const
  {
    return NeumanBoundaryType( "boundary function", bndN_, gridPart_, 5 );
  }

  // return Fem :: Function for right hand side
  RightHandSideType rightHandSide(  ) const
  {
    return RightHandSideType( "right hand side", rhs_, gridPart_, 5 );
  }

  //! penalty parameter for DG methods
  double penalty() const
  {
    return penalty_;
  }

  // Public custom attributes for the EIT model
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
      else if( id == bndD )
      {
        // call dirichlet boudary data of implementation
        impl_.g( x, ret );
      }
      else if( id == bndN )
      {
        // call dirichlet boudary data of implementation
        impl_.n( x, ret );
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
  FunctionWrapper<bndD> bndD_;
  FunctionWrapper<bndN> bndN_;
  double penalty_;
  // Protected custom attributes for the EIT model
  const SigmaFunctionType sigma_;
  const SigmaFunctionType elementID_;
  const bool uniformCond_;
  const bool uniformCondValue_;
  const bool assignCond_;
  const bool havePerturbation_;
  const double pertRadius_;
  const bool multORabs_;
  const double value_;
  std::vector<double> pertcenter_;
  std::vector<double> conductivities_;
};

#endif // #ifndef ELLIPTC_MODEL_HH
