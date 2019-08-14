// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#include <iostream>

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include "dirichletconstraints.hh"

#include <dune/fem/function/localfunction/temporary.hh>
// EllipticOperator
// ----------------

template< class DiscreteFunction, class Model >
struct EllipticOperator
: public virtual Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  //! type of Dirichlet constraints 
  typedef Dune::DirichletConstraints< ModelType, DiscreteFunctionSpaceType > ConstraintsType;

public:
  //! contructor 
  EllipticOperator ( const ModelType &model, const DiscreteFunctionSpaceType &space )
  : model_( model )
  , constraints_( model, space )
  {}
      
  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
    // set boundary values for solution 
    constraints()( func, u );
  }

  //! application operator
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

protected:
  const ModelType &model () const { return model_; }
  const ConstraintsType &constraints () const { return constraints_; }

private:
  ModelType model_;
  ConstraintsType constraints_;
};

/*********************************************************/
/***                 NEW FOR LESSON 2                  ***/
/*********************************************************/
// DifferentiableEllipticOperator
// ------------------------------

template< class JacobianOperator, class Model >
struct DifferentiableEllipticOperator
: public EllipticOperator< typename JacobianOperator::DomainFunctionType, Model >,
  public Dune::Fem::DifferentiableOperator< JacobianOperator >
{
  typedef EllipticOperator< typename JacobianOperator::DomainFunctionType, Model > BaseType;

  typedef JacobianOperator JacobianOperatorType;

public:
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename EntityType::EntitySeed EntitySeed;
  // typedef typename EntityType::EntityPointer EntityPointer;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename GridType::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType IdType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef Dune::Fem::Stencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;

public:
  //! contructor 
  DifferentiableEllipticOperator ( const ModelType &model, const DiscreteFunctionSpaceType &space, ElectrodePositions &electrodes )
  : BaseType( model, space ),
    electrodes_(electrodes),
    areaE_(electrodes_.electrode_positions_.size(),0),
    electrodeElements_(electrodes_.electrode_positions_.size()),
    stencil_(space,space)
  {
    const GridPartType& gridPart = space.gridPart();
    const GlobalIdSet &globalIdSet = space.grid().globalIdSet();

    const IteratorType end = space.end();
    for( IteratorType it = space.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      bool element_already_added_to_list = false;

      if( !entity.hasBoundaryIntersections() )
      {
        stencil_.fill(entity,entity);
        continue;
      }

      const IntersectionIteratorType endiit = gridPart.iend( entity );
      int intersection_counter = 0;
      for ( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != endiit ; ++iit, ++intersection_counter )
      {
        const IntersectionType& intersection = *iit ;

        if( intersection.boundary() )
        {
          typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
          const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

          int electrodenumber = 0;
          if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
          {
            std::vector<int> globalIndices(4);
            IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
            abc.getKey().extractKey(globalIndices);
            electrodenumber = electrodes_.electrodeNumber_node(globalIndices);
          }
          else
            electrodenumber = electrodes_.electrodeNumber(intersectionGeometry.center());

          if (electrodenumber)
          {
            areaE_[electrodenumber-1] += intersectionGeometry.volume();
            if (!element_already_added_to_list)
            {
              electrodeElements_[electrodenumber-1].push_back(entity.seed()); // remember that this entity belongs to electrode
              element_already_added_to_list = true;
            }
            
            // fill the stencil
            for( unsigned int electrodeiterator = 0; electrodeiterator < electrodeElements_[electrodenumber-1].size(); ++electrodeiterator)
            {
              // EntityPointer ep = gridPart.grid().entityPointer(electrodeElements_[electrodenumber-1][electrodeiterator]);
              // const EntityType &entitySecond = *ep;
              const EntityType &entitySecond = gridPart.grid().entity(electrodeElements_[electrodenumber-1][electrodeiterator]);
              // fill in the 'diagonal' part
              stencil_.fill( entity, entitySecond );
              stencil_.fill( entitySecond, entity);
            }
          }
          else
            stencil_.fill( entity, entity);
        }
        else
          stencil_.fill( entity, entity);
      }
    }
    // communicate the electrode areas. not necessary, since each electrode is on one process only.
    MPI_Allreduce(MPI_IN_PLACE, &areaE_[0], areaE_.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  }
      
  //! method to setup the jacobian of the operator for storage in a matrix 
  void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

  const std::vector< std::vector<EntitySeed> > &entitySeedVector() const
  { return electrodeElements_; }
  
  const std::vector< double > &electrodeAreas() const
  { return areaE_; }

protected:
  using BaseType::model;
  using BaseType::constraints;
  std::string electrodes_string_;
  ElectrodePositions electrodes_;
  std::vector<double> areaE_;
  std::vector< std::vector<EntitySeed> > electrodeElements_;
  StencilType stencil_;
};
/*********************************************************/

// Implementation of EllipticOperator
// ----------------------------------

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  // have a look at the first poisson howto if you want to use this operator

}

// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template< class JacobianOperator, class Model >
void DifferentiableEllipticOperator< JacobianOperator, Model >
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BaseFunctionSetType;

  const DiscreteFunctionSpaceType &dfSpace = u.space();
  const GridPartType& gridPart = dfSpace.gridPart();
  const GlobalIdSet &globalIdSet = u.space().grid().globalIdSet();

  // read in if we are grounding by an epsilon weighted mass term (if yes, need to specify h^2)
  const bool epsmassgrounding = Dune::Fem::Parameter::getValue< bool >( "ground.epsmass" );
  const double hsquared = Dune::Fem::Parameter::getValue< double >( "ground.hsquared" );

  jOp.reserve(stencil_);
  jOp.clear();

  std::vector< typename LocalFunctionType::RangeType > phi( dfSpace.blockMapper().maxNumDofs() );
  std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace.blockMapper().maxNumDofs() );

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    // std::cout << "assembling next element" << std::endl;
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const BaseFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    // find sigma for this entity
    RangeType sigmavalue(0);
    model().sigmaValue( entity, sigmavalue );

    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      // evaluate all basis functions at given quadrature point 
      baseSet.evaluateAll( quadrature[ pt ], phi );

      // evaluate jacobians of all basis functions at given quadrature point
      baseSet.jacobianAll( quadrature[pt], dphi );

      RangeType aphi( 0 );
      JacobianRangeType adphi( 0 );
      for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
      {
    	adphi = dphi[ localCol ];
        adphi *= sigmavalue;
        if (epsmassgrounding)
          aphi = phi[localCol]*hsquared;
        else
          aphi = 0.0;
        // get column object and call axpy method 
        jLocal.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
      }
    }

    //////// BEGIN SURFACE INTEGRALS
    if( !entity.hasBoundaryIntersections() )
      continue;

    const IntersectionIteratorType endiit = gridPart.iend( entity );
    int intersection_counter=0;
    for ( IntersectionIteratorType iit = gridPart.ibegin( entity );
          iit != endiit ; ++iit, ++intersection_counter )
    {
      const IntersectionType& intersection = *iit ;

      if( intersection.boundary() )
      {
        typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
        const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

        int electrodenumber = 0;
        if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
        {
          std::vector<int> globalIndices(4);
          IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
          abc.getKey().extractKey(globalIndices);
          electrodenumber = electrodes_.electrodeNumber_node(globalIndices);
        }
        else
          electrodenumber = electrodes_.electrodeNumber(intersectionGeometry.center());

        if (!electrodenumber)
          continue;

        FaceQuadratureType faceQuadInside(gridPart, intersection, 2*dfSpace.order(),
                                          FaceQuadratureType::INSIDE);

        const size_t numFaceQuadPoints = faceQuadInside.nop();
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const typename FaceQuadratureType::LocalCoordinateType &x = faceQuadInside.localPoint( pt );

          const double weight = faceQuadInside.weight( pt ) * intersectionGeometry.integrationElement(x);

          // evaluate all basis functions for quadrature point pt
          baseSet.evaluateAll( faceQuadInside[ pt ], phi );

          for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
          {
            jLocal.column( localCol ).axpy( phi, 
                RangeType(phi[localCol]/electrodes_.contactImpedance(electrodenumber)/areaE_[electrodenumber-1]),
                weight );
          }
        }

        /////////// Second iteration over all elements with surfaces under the selected electrode
        for( unsigned int electrodeiterator = 0; electrodeiterator < electrodeElements_[electrodenumber-1].size(); ++electrodeiterator)
        {
          // EntityPointer ep = gridPart.grid().entityPointer(electrodeElements_[electrodenumber-1][electrodeiterator]);
          // const EntityType &entitySecond = *ep;
          const EntityType &entitySecond = gridPart.grid().entity(electrodeElements_[electrodenumber-1][electrodeiterator]);
          
          LocalMatrixType jLocalToGlobal = jOp.localMatrix( entitySecond, entity );

          const IntersectionIteratorType endiitSecond = gridPart.iend( entitySecond );
          int intersectionSecond_counter=0;
          for ( IntersectionIteratorType iitSecond = gridPart.ibegin( entitySecond );
                iitSecond != endiitSecond ; ++iitSecond, ++intersectionSecond_counter )
          {
            const IntersectionType& intersectionSecond = *iitSecond ;

            if( intersectionSecond.boundary() )
            {
              typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
              const IntersectionGeometryType &intersectionGeometrySecond = intersectionSecond.geometry();

              int electrodenumberSecond = 0;
              if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
              {
                std::vector<int> globalIndices(4);
                IdType abc = globalIdSet.subId(entitySecond, intersectionSecond_counter, 1);
                abc.getKey().extractKey(globalIndices);
                electrodenumberSecond = electrodes_.electrodeNumber_node(globalIndices);
              }
              else
                electrodenumberSecond = electrodes_.electrodeNumber(intersectionGeometrySecond.center());

              if (electrodenumberSecond != electrodenumber)
                continue;

              typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TempLocalFunction;
              TempLocalFunction sumUlocal(dfSpace,entitySecond);
              sumUlocal.clear();

              //// calculate sumU contribution locally
              FaceQuadratureType faceQuadInsideSecond(gridPart, intersectionSecond, 2*dfSpace.order(),
                                                      FaceQuadratureType::INSIDE);

              const size_t numFaceQuadPointsSecond = faceQuadInsideSecond.nop();
              for( size_t ptSecond = 0; ptSecond < numFaceQuadPointsSecond; ++ptSecond )
              {
                const typename FaceQuadratureType::LocalCoordinateType &xSecond = faceQuadInsideSecond.localPoint( ptSecond );

                const double weightSecond = faceQuadInsideSecond.weight( ptSecond )*intersectionGeometrySecond.integrationElement(xSecond);

                // evaluate all basis functions at given quadrature point
                sumUlocal.axpy(faceQuadInsideSecond[ptSecond],RangeType(weightSecond));
              }

              // write facet-to-facet part of third term
              for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
              {
                const typename FaceQuadratureType::LocalCoordinateType &x = faceQuadInside.localPoint( pt );

                const double weight = faceQuadInside.weight( pt ) * intersectionGeometry.integrationElement(x);

                // evaluate all basis functions for quadrature point pt
                baseSet.evaluateAll( faceQuadInside[ pt ], phi );

                for( unsigned int localColSecond = 0; localColSecond < numBaseFunctions; ++localColSecond )
                {
                  jLocalToGlobal.column(localColSecond).axpy(phi,
                      RangeType(sumUlocal[localColSecond]/electrodes_.contactImpedance(electrodenumber)/areaE_[electrodenumber-1]/areaE_[electrodenumber-1]), 
                      -weight);
                }
              }
            }
          } // end iteration over facets of element
        } // end iteration of elements under same electrode
      }
    } // end surface integral
  } // end grid traversal

  // apply constraints to matrix operator 
  constraints().applyToOperator( jOp );
  jOp.communicate();
  const bool writeMatrixToMatlab = Dune::Fem::Parameter::getValue< bool >( "output.write_matrix" );
#if HAVE_PETSC && WANT_PETSC
  if (writeMatrixToMatlab)
  {
    if (Dune::Fem::MPIManager::size() == 1) { jOp.view(); }     // write matrix to a binary file which can be read with matlab with function PetscBinaryRead in petsc/bin/matlab
    else { cout << "You can only write the matrix in serial" << endl; }
  }
#endif
}
/*********************************************************/

#endif // #ifndef ELLIPTIC_HH
