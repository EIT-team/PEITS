// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef RHS_HH
#define RHS_HH

#include <dune/common/dynvector.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include "matlabhelpers.hh"
#include "electrode_helpers.hh"

#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/function/localfunction/temporary.hh>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

using namespace std;

// used in calcElecPot. simply returns the signum of val
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


// assembleRHS
// -----------
template< class Function, class DiscreteFunction, class EntitySeed >
void assembleRHS ( const Function &function, DiscreteFunction &rhs, int pattern_no, ElectrodePositions &electrodes, CurrentProtocol &current_protocol, const std::vector< std::vector<EntitySeed> > &electrodeElements, const std::vector<double> &areaE )
{
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  // typedef typename EntityType::EntityPointer EntityPointer;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename GridType::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType IdType;

  rhs.clear();

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();
  const GlobalIdSet &globalIdSet = rhs.space().grid().globalIdSet();

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;

    if ( !entity.hasBoundaryIntersections() )
      continue;

    LocalFunctionType rhsLocal = rhs.localFunction( entity );

    const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
    double factor;
    int intersection_counter=0;
    for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit, ++intersection_counter ) // looping over intersections
    {
      const IntersectionType &intersection = *iit;

      if ( ! intersection.boundary() )
        continue;

      typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
      const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

      int number = 0;
      if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
      {
        std::vector<int> globalIndices(4);
        IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
        abc.getKey().extractKey(globalIndices);
        number = electrodes.electrodeNumber_node(globalIndices);
      }
      else
        number = electrodes.electrodeNumber(intersectionGeometry.center());

      if( number )
        factor = current_protocol.current(number, pattern_no) / areaE[number-1];
      else
      {
        factor = 0;
        continue;
      }

      FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, dfSpace.order(), FaceQuadratureType::INSIDE );
      const size_t numQuadraturePoints = quadInside.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        // evaluate f
        typename Function::RangeType
        f( quadInside.weight( pt ) * intersectionGeometry.integrationElement(x) * factor );
        rhsLocal.axpy( quadInside[ pt ], f );
      }
    }
  }
  rhs.communicate();
}

// calculate electrode potentials
// -----------
template< class DiscreteFunction, class EntitySeed >
void calculateElecPot ( const DiscreteFunction &solution, const int pattern_index, ElectrodePositions &electrodes, CurrentProtocol &current_protocol, const std::vector< std::vector<EntitySeed> > &electrodeElements, const std::vector<double> &areaE, vector<double> &electrode_potentials )
{
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename EntityType::EntityPointer EntityPointer;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename GridType::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType IdType;

  const DiscreteFunctionSpaceType &dfSpace = solution.space();
  const GlobalIdSet &globalIdSet = solution.space().grid().globalIdSet();

  if( Dune::Fem::MPIManager::rank() == 0 )
  {
    for( unsigned int i = 0 ; i<electrodes.electrode_positions_.size() ; ++i )
  	  electrode_potentials[i] = current_protocol.current(i+1, pattern_index) * electrodes.contactImpedance(i+1);
  }

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;

    if ( !entity.hasBoundaryIntersections() )
      continue;

    LocalFunctionType solutionLocal = solution.localFunction( entity );

    const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
    double factor;
    int intersection_counter=0;
    for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit, ++intersection_counter ) // looping over intersections
    {
      const IntersectionType &intersection = *iit;
      if ( ! intersection.boundary() )
        continue;

      typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
      const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

      unsigned int number = 0;
      if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
      {
        std::vector<int> globalIndices(4);
        IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
        abc.getKey().extractKey(globalIndices);
        number = electrodes.electrodeNumber_node(globalIndices);
      }
      else
        number = electrodes.electrodeNumber(intersectionGeometry.center());

      if( number )
        factor = 1 / areaE[number-1];
      else
      {
        factor = 0;
        continue;
      }

      FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, 2*dfSpace.order(), FaceQuadratureType::INSIDE );
      const size_t numQuadraturePoints = quadInside.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

        // evaluate f
        RangeType f(0);
        f= 0;
        solutionLocal.evaluate( quadInside[ pt ], f );

        // multiply by quadrature weight
        f *= quadInside.weight(pt) * intersectionGeometry.integrationElement(x) * factor;

        electrode_potentials[number-1] += f;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &electrode_potentials[0], electrode_potentials.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}


// calculate and write measured voltages
// -----------
void writeElecPot ( const vector<double> &electrode_potentials, const int measure_index, const char time_buf[], CurrentProtocol &current_protocol, const double drive_factor, const double measure_factor, const std::vector<double> &areaE, const int protocol_step )
{
  vector<double> el_pots (electrode_potentials);
 
  if( Dune::Fem::MPIManager::rank() == 0 )
  {
	cout.precision(20);
	double temp = 0;
    for( unsigned int i=0 ; i<el_pots.size() ; ++i )
    {
      el_pots[i] *= drive_factor;
   	  if (Dune::Fem::Parameter::getValue< bool >( "output.print_electrode_voltages" ))
        cout << "    Electrode " << i+1 << " = " << el_pots[i] << "V ; Area = " << areaE[i] << endl;
   	  temp += el_pots[i]*current_protocol.current(i+1, measure_index)*measure_factor;
    }
    cout << "    Measured voltage = " << temp << endl;
    cout.precision(10);

    int noProtocols = current_protocol.current_protocol_.size();
    int mode = 0;
    if (!protocol_step)
       mode = noProtocols;
    else if (protocol_step == noProtocols-1)
       mode = -1;

    // write electrode voltages to matlab readable file
    Dune::Fem::MatlabHelper mhelp;
    std::string filename("electrodevoltages.bin");
    filename.insert(17,time_buf);
    const string outputpath = Dune::Fem::Parameter::getValue< string >( "fem.io.outputpath" );
    filename.insert(0,outputpath);
    if (Dune::Fem::Parameter::getValue< bool >("fem.io.write_only_measured_voltage"))
    { 
      vector<double> measuredV {temp};
      mhelp.saveElectrodeVoltagesBinary(filename.c_str(), measuredV, mode);
    }
    else
      mhelp.saveElectrodeVoltagesBinary(filename.c_str(), el_pots, mode);
  }
}


// structure for sides
struct Side {
  vector<int> globalID;
  vector<int> IDmapper;
  vector<Dune::FieldVector<double,3>> coords;
};


// function to sort sides
bool sortMySides (Side i, Side j)
{
  if (i.globalID[0] < j.globalID[0])
    return true;
  else
  {
    if ( (i.globalID[0] == j.globalID[0]) && (i.globalID[1] < j.globalID[1]) )
      return true;
    else
      return false;
  }
}


template< class NodalFunction, class PiecewiseFunction >
struct JacobianRowCalculator
{
  typedef typename NodalFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename NodalFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename EntityType::EntitySeed EntitySeed;
  // typedef typename EntityType::EntityPointer EntityPointer;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename GridType::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType IdType;

  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > TemporaryLocalMatrix;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BaseFunctionSetType;
  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TempLocalFunction;

  JacobianRowCalculator( const NodalFunction &dummyFunction, ElectrodePositions &electrodes, const vector<vector<EntitySeed>> &electrodeElements, const std::vector<double> &areaE )
    : dfSpace_(dummyFunction.space().gridPart()),
      electrodes_(electrodes)
  {
    std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace_.blockMapper().maxNumDofs() );
  
    const GridPartType& gridPart = dfSpace_.gridPart();

    const GlobalIdSet &globalIdSet = dummyFunction.space().grid().globalIdSet();

    // container for all sides that are part of an electrode
    int no_elecs = electrodeElements.size();
    do_electrode_jacobian_ = Dune::Fem::Parameter::getValue< bool >("fem.io.do_electrode_jacobian");
    vector<vector<Side>> sides_containers(no_elecs);
    vector<Dune::FieldVector<double,3>> surface_coordinates;
    int dim_srf_coord_x, dim_srf_coord_y;
    vector<Dune::FieldVector<double,3>> field_vector_x, field_vector_y;

    // obtain required size for sparse boost matrix
    int matrix_size = dummyFunction.size();

    if ( do_electrode_jacobian_ )
    {
      int error = 0;
      const string srf_coord_file = Dune::Fem::Parameter::getValue< string >("surface.coordinates");
      FILE *F = 0;
      F = fopen(srf_coord_file.c_str(), "r");
      if (F == NULL)
      {
        cout << "Surface coordinate file could not be loaded!" << endl;
        error = 1;
        return;
      }
      if (!fscanf(F, "%d, %d \n", &dim_srf_coord_x, &dim_srf_coord_y))
        cout << "Size of surface coordinate system is not specified in file" << endl;
	  while(!std::feof(F))
	  {
		double c1, c2, c3;
		error = fscanf(F, "%lf,%lf,%lf\n", &c1, &c2, &c3);
		if (!error)
		  cout << "File could not be read: " << F << endl;
		Dune::FieldVector<double,3> coord;
		coord[0]=c1;
		coord[1]=c2;
		coord[2]=c3;
		surface_coordinates.push_back(coord);
	  }
	  fclose(F);

      // find the field vectors for all electrodes
      for (int i=0; i<no_elecs; ++i)
      {
        Dune::FieldVector<double,3> elec_coord;
        elec_coord[0] = electrodes_.electrode_positions_[i][0]; elec_coord[1] = electrodes_.electrode_positions_[i][1]; elec_coord[2] = electrodes_.electrode_positions_[i][2];
        int minimum = 0;
        double value = 1e10;
        for (int j=0; j<dim_srf_coord_x*dim_srf_coord_y/3; ++j)
        {
          double dist = (elec_coord-surface_coordinates[j]).two_norm();
          if (dist<value)
          {
            value=dist;
            minimum=j;
          }
        }

        Dune::FieldVector<double,3> pos_x = elec_coord-surface_coordinates[minimum+dim_srf_coord_x];
        Dune::FieldVector<double,3> neg_x = elec_coord-surface_coordinates[minimum-dim_srf_coord_x];
        Dune::FieldVector<double,3> pos_y = elec_coord-surface_coordinates[minimum+1];
        Dune::FieldVector<double,3> neg_y = elec_coord-surface_coordinates[minimum-1];

        if ( pos_x.two_norm() < neg_x.two_norm() )
          field_vector_x.push_back(surface_coordinates[minimum+dim_srf_coord_x] - surface_coordinates[minimum]);
        else
          field_vector_x.push_back(surface_coordinates[minimum] - surface_coordinates[minimum-dim_srf_coord_x]);
        if ( pos_y.two_norm() < neg_y.two_norm() )
          field_vector_y.push_back(surface_coordinates[minimum+1] - surface_coordinates[minimum]);
        else
          field_vector_y.push_back(surface_coordinates[minimum] - surface_coordinates[minimum-1]);
        field_vector_x[i] /= field_vector_x[i].two_norm();
        field_vector_y[i] /= field_vector_y[i].two_norm();

        // create template electrode jacobian for this electrode
        elec_jacs_x_.push_back(make_shared<boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major>> (matrix_size+no_elecs,matrix_size+no_elecs));
        elec_jacs_y_.push_back(make_shared<boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major>> (matrix_size+no_elecs,matrix_size+no_elecs));
      }
    }

    // ASSEMBLE THE LOCAL STIFFNESS MATRIX FOR EACH ELEMENT
    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      stiffnessMatrices_.push_back( new TemporaryLocalMatrix (dfSpace_, dfSpace_, entity, entity) );
      stiffnessMatrices_.back().clear();

      const BaseFunctionSetType &baseSet = stiffnessMatrices_.back().domainBasisFunctionSet();
      const unsigned int numBaseFunctions = baseSet.size();

      mapper_.push_back( std::vector<int>( numBaseFunctions ) );
      dfSpace_.blockMapper().mapEach(entity, Dune::Fem::AssignFunctor< std::vector<int> >( mapper_.back() ) );

      // If element is part of an electrode, we have to compute the electrode boundary Jacobian
      // For this, we first save all sides in order to find the boundary of the electrode:
      if( entity.hasBoundaryIntersections() && do_electrode_jacobian_ )
      {
        int intersection_counter = 0;
        const IntersectionIteratorType endiit = gridPart.iend( entity );
        for ( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != endiit ; ++iit, ++intersection_counter  )
        {
          const IntersectionType& intersection = *iit ;

          if( intersection.boundary() )
          {
            typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
            const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

            int electrodenumber = electrodes_.electrodeNumber(intersectionGeometry.center());
            if (!electrodenumber)
              continue;

            std::vector<int> globalIndices(4), elementGID(4), side_mapper(3);
            std::vector<Dune::FieldVector<double,3>> coordinate_mapper(3);
            IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
            IdType abcd = globalIdSet.id(entity);
            abc.getKey().extractKey(globalIndices);
            abcd.getKey().extractKey(elementGID);

            for (int i=0; i<3; i++)
            {
              for (int j=0; j<4; j++)
              {
                if( globalIndices[i] == elementGID[j] )
                {
                  side_mapper[i] = mapper_.back()[j];
                  coordinate_mapper[i] = geometry.corner(j);
                }
              }
            }

            // write all three edges into the sides_container with their corresponding "local global index"
            // (we use the fact that the indices of the triangle are already sorted - otherwise a globalIndices[0]<globalIndices[1] sorting would be necessary)
            sides_containers[electrodenumber-1].push_back( Side { vector<int> {globalIndices[0], globalIndices[1]} , vector<int> {side_mapper[0], side_mapper[1]} , 
                                                                  vector<Dune::FieldVector<double,3>> {coordinate_mapper[0], coordinate_mapper[1], coordinate_mapper[2]} } );
            sides_containers[electrodenumber-1].push_back( Side { vector<int> {globalIndices[0], globalIndices[2]} , vector<int> {side_mapper[0], side_mapper[2]} , 
                                                                  vector<Dune::FieldVector<double,3>> {coordinate_mapper[0], coordinate_mapper[2], coordinate_mapper[1]} } );
            sides_containers[electrodenumber-1].push_back( Side { vector<int> {globalIndices[1], globalIndices[2]} , vector<int> {side_mapper[1], side_mapper[2]} , 
                                                                  vector<Dune::FieldVector<double,3>> {coordinate_mapper[1], coordinate_mapper[2], coordinate_mapper[0]} } );

            // write the mapper ID's into the corresponding vector
            for (int i=0; i<3; i++)
            {
              if (find(elec_nodes_.begin(), elec_nodes_.end(), side_mapper[i]) == elec_nodes_.end())
                elec_nodes_.push_back(side_mapper[i]);
            }
          }
        }
      }

      QuadratureType quadrature( entity, 2*dfSpace_.order() );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        // evaluate jacobians of all basis functions at given quadrature point
        baseSet.jacobianAll( quadrature[pt], dphi );

        JacobianRangeType adphi( 0 );
        for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
        {
      	  adphi = dphi[ localCol ];

          // get column object and call axpy method 
          stiffnessMatrices_.back().column( localCol ).axpy( dphi, adphi, weight );
        }
      }
    }

    if (!do_electrode_jacobian_)
      return; // the "normal" jacobian matrix is now set up

    // let us now find out which sides belong the the electrode boundary
    for (int electrodenumber = 0; electrodenumber<no_elecs; electrodenumber++)
    {
      if ( sides_containers[electrodenumber].empty() )
        continue;

      sort(sides_containers[electrodenumber].begin(),sides_containers[electrodenumber].end(),sortMySides);

      int shift = 0;
      for (unsigned int side = 0; side-shift<sides_containers[electrodenumber].size()-1; side++)
      {
        if (sides_containers[electrodenumber][side-shift].globalID == sides_containers[electrodenumber][side+1-shift].globalID)
        {
          // remove these two entries since they are not part of electrode boundary
          sides_containers[electrodenumber].erase(sides_containers[electrodenumber].begin() + side-shift, sides_containers[electrodenumber].begin() + side-shift + 2);
          shift = shift + 2; side++;
        }
      }
    }

    // now to the second iteration, where we actually set up the electrode Jacobians
    for( int electrode = 0; electrode<no_elecs; electrode++)
    {
      for( unsigned int electrodeiterator = 0; electrodeiterator < electrodeElements[electrode].size(); ++electrodeiterator)
      {
        // EntityPointer ep = gridPart.grid().entityPointer(electrodeElements[electrode][electrodeiterator]);        
        // const EntityType &entity = *ep;
        const EntityType &entity = gridPart.grid().entityPointer(electrodeElements[electrode][electrodeiterator]);

        int intersection_counter = 0;
        const IntersectionIteratorType endiit = gridPart.iend( entity );
        for ( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != endiit ; ++iit, ++intersection_counter  )
        {
          const IntersectionType& intersection = *iit ;

          if( intersection.boundary() )
          {
            typedef typename IntersectionType :: Geometry  IntersectionGeometryType;
            const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

            int electrodenumber = electrodes_.electrodeNumber(intersectionGeometry.center());
            if (electrodenumber-1!=electrode)
              continue;

            std::vector<int> globalIndices(4);
            IdType abc = globalIdSet.subId(entity, intersection_counter, 1);
            abc.getKey().extractKey(globalIndices);

            vector<int> delete_me;
            // test if surface is part of electrode boundary
            for (unsigned int side = 0; side < sides_containers[electrodenumber-1].size(); side++)
            {
              int count = 0;
              Dune::FieldVector<double,3> internal_node(0), side_node_1(0), side_node_2(0), outward_vector(0), side_vector(0), perpendicular_vector(0), normal_vector(0);

              // analyse first node
              if ( (sides_containers[electrodenumber-1][side].globalID[0] != globalIndices[0]) && (sides_containers[electrodenumber-1][side].globalID[1] != globalIndices[0]) )
              {
                count = 1;
                internal_node = sides_containers[electrodenumber-1][side].coords[2];
              }
              else
              {
                if ( globalIndices[0] == sides_containers[electrodenumber-1][side].globalID[0] )
                  side_node_1 = sides_containers[electrodenumber-1][side].coords[0];
                else
                  side_node_2 = sides_containers[electrodenumber-1][side].coords[1];
              }

              // analyse second node
              if ( (sides_containers[electrodenumber-1][side].globalID[0] != globalIndices[1]) && (sides_containers[electrodenumber-1][side].globalID[1] != globalIndices[1]) )
              {
                if (count == 1) 
                  continue;
                else
                {
                  count = 1;
                  internal_node = sides_containers[electrodenumber-1][side].coords[2];
                }
              }
              else
              {
                if ( globalIndices[1] == sides_containers[electrodenumber-1][side].globalID[0] )
                  side_node_1 = sides_containers[electrodenumber-1][side].coords[0];
                else
                  side_node_2 = sides_containers[electrodenumber-1][side].coords[1];
              }

              // analyse third node
              if ( (sides_containers[electrodenumber-1][side].globalID[0] != globalIndices[2]) && (sides_containers[electrodenumber-1][side].globalID[1] != globalIndices[2]) )
              {
                if (count == 1) 
                  continue;
                else
                {
                  internal_node = sides_containers[electrodenumber-1][side].coords[2];
                }
              }
              else
              {
                if ( globalIndices[2] == sides_containers[electrodenumber-1][side].globalID[0] )
                  side_node_1 = sides_containers[electrodenumber-1][side].coords[0];
                else
                  side_node_2 = sides_containers[electrodenumber-1][side].coords[1];
              }

              // store this side for deletion from sides_containers
              delete_me.push_back(side);

              // set up the matrix for this side
              outward_vector = (side_node_1-internal_node) + (side_node_2-internal_node);
              side_vector = side_node_1 - side_node_2;
              perpendicular_vector = outward_vector; perpendicular_vector.axpy(-(side_vector.dot(outward_vector)),side_vector);
              normal_vector = perpendicular_vector; normal_vector /= perpendicular_vector.two_norm();

              vector<int> nodes = {sides_containers[electrodenumber-1][side].IDmapper[0], sides_containers[electrodenumber-1][side].IDmapper[1], matrix_size+electrodenumber-1};

              double multiplier_x = -1/(electrodes_.contactImpedance(electrodenumber) * areaE[electrodenumber-1]) * side_vector.two_norm() * (normal_vector*field_vector_x[electrodenumber-1]);
              double multiplier_y = -1/(electrodes_.contactImpedance(electrodenumber) * areaE[electrodenumber-1]) * side_vector.two_norm() * (normal_vector*field_vector_y[electrodenumber-1]);

              // now we write the matrix: multiplier* [1/3, 1/6, -1/2; 1/6, 1/3, -1/2; -1/2, -1/2, 1]
              (*elec_jacs_x_[electrodenumber-1])(nodes[0],nodes[0]) += multiplier_x/3; (*elec_jacs_y_[electrodenumber-1])(nodes[0],nodes[0]) += multiplier_y/3;
              (*elec_jacs_x_[electrodenumber-1])(nodes[0],nodes[1]) += multiplier_x/6; (*elec_jacs_y_[electrodenumber-1])(nodes[0],nodes[1]) += multiplier_y/6;
              (*elec_jacs_x_[electrodenumber-1])(nodes[0],nodes[2]) -= multiplier_x/2; (*elec_jacs_y_[electrodenumber-1])(nodes[0],nodes[2]) -= multiplier_y/2;
              (*elec_jacs_x_[electrodenumber-1])(nodes[1],nodes[0]) += multiplier_x/6; (*elec_jacs_y_[electrodenumber-1])(nodes[1],nodes[0]) += multiplier_y/6;
              (*elec_jacs_x_[electrodenumber-1])(nodes[1],nodes[1]) += multiplier_x/3; (*elec_jacs_y_[electrodenumber-1])(nodes[1],nodes[1]) += multiplier_y/3;
              (*elec_jacs_x_[electrodenumber-1])(nodes[1],nodes[2]) -= multiplier_x/2; (*elec_jacs_y_[electrodenumber-1])(nodes[1],nodes[2]) -= multiplier_y/2;
              (*elec_jacs_x_[electrodenumber-1])(nodes[2],nodes[0]) -= multiplier_x/2; (*elec_jacs_y_[electrodenumber-1])(nodes[2],nodes[0]) -= multiplier_y/2;
              (*elec_jacs_x_[electrodenumber-1])(nodes[2],nodes[1]) -= multiplier_x/2; (*elec_jacs_y_[electrodenumber-1])(nodes[2],nodes[1]) -= multiplier_y/2;
              (*elec_jacs_x_[electrodenumber-1])(nodes[2],nodes[2]) += multiplier_x;   (*elec_jacs_y_[electrodenumber-1])(nodes[2],nodes[2]) += multiplier_y;

            }
            // delete the used entries
            for (unsigned int side = 0; side<delete_me.size(); side++)
              sides_containers[electrodenumber-1].erase(sides_containers[electrodenumber-1].begin() + delete_me[side] - side);
          }
        }
      }
    }
  }


  void getJacobianRow ( const NodalFunction &drive, const NodalFunction &measure, PiecewiseFunction &jacobianRow,
                          const double factor )
  {
    jacobianRow.clear();
    
    for (unsigned int index=0;index<mapper_.size();++index)
    {
      Dune::DynamicVector<double> dlocal( mapper_[0].size(),0 ),
                  tmp( mapper_[0].size(),0 );
      for (unsigned int i=0;i<mapper_[index].size();++i)
        dlocal[i] = (*drive.block(mapper_[index][i]))[0];
      stiffnessMatrices_[index].multiplyAdd(dlocal, tmp);
      for (unsigned int i=0;i<mapper_[index].size();++i)
        (*jacobianRow.block(index))[0] -= tmp[i] * (*measure.block(mapper_[index][i]))[0]*factor;
    }

    jacobianRow.communicate();
  }


  void getElecJacRow( const NodalFunction &drive, const vector<double> &elec_drive, const NodalFunction &measure, const vector<double> &elec_meas, vector<double> &elecJacRow, const double factor)
  {
    // fill a sparse vector with the values we actually require
    boost::numeric::ublas::compressed_vector<double> boost_drive ( drive.size()+elec_drive.size() , elec_nodes_.size()+elec_drive.size() );
    boost::numeric::ublas::compressed_vector<double> boost_meas ( measure.size()+elec_meas.size() , elec_nodes_.size()+elec_meas.size() );
    for (unsigned int i=0; i<elec_nodes_.size(); i++)
    {
      boost_drive(elec_nodes_[i]) = (*drive.block(elec_nodes_[i]))[0];
      boost_meas(elec_nodes_[i]) = (*measure.block(elec_nodes_[i]))[0];
    }
    for (unsigned int i=0; i<elec_drive.size(); i++)
    {
      boost_drive(drive.size()+i) = elec_drive[i];
      boost_meas(measure.size()+i) = elec_meas[i];
    }

    double entry_x, entry_y;
    for (unsigned int electrode=0; electrode<elec_drive.size(); ++electrode)
    {
      // compute the jacobian entries by left and right multiplication of the template jacobian
      entry_x = factor * boost::numeric::ublas::prec_inner_prod(boost_drive, boost::numeric::ublas::prec_prod(*elec_jacs_x_[electrode], boost_meas));
      entry_y = factor * boost::numeric::ublas::prec_inner_prod(boost_drive, boost::numeric::ublas::prec_prod(*elec_jacs_y_[electrode], boost_meas));

      // MPI communication (additive reduction - since only one process actually has the value)
      MPI_Allreduce(MPI_IN_PLACE, &entry_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &entry_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
      // add the entries to the jacobian matrix
      elecJacRow.push_back(entry_x); elecJacRow.push_back(entry_y);
    }
  }


  boost::ptr_vector< TemporaryLocalMatrix > stiffnessMatrices_;
  std::vector< std::vector<int> > mapper_;
  const DiscreteFunctionSpaceType dfSpace_;
  ElectrodePositions electrodes_;
  bool do_electrode_jacobian_;
  vector<shared_ptr<boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major>>> elec_jacs_x_, elec_jacs_y_; // template electrode position matrices
  vector<int> elec_nodes_; // all nodes involved in the electrode position jacobian matrix
};

template< class Function >
void shiftAverage ( Function &solution )
{
  typedef typename Function::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename Function::LocalFunctionType LocalFunctionType;
  typedef typename Function::DofIteratorType DofIteratorType ; 

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  const DiscreteFunctionSpaceType &dfSpace = solution.space();

  double average[] = {0,0}; // average,volume

  IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();
    average[1] += geometry.volume();

    LocalFunctionType localSol = solution.localFunction( entity);

    QuadratureType quadrature( entity, 2*dfSpace.order() );

    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
    	const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
    	const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );
        RangeType value;
        localSol.evaluate( quadrature[pt], value );
        average[0] += value[0] * weight;
    }
  }

  dfSpace.gridPart().grid().comm().sum( average,2 );
  average[0] /= average[1];
  if (dfSpace.gridPart().grid().comm().rank() == 0)
    std::cout << "average=" << average[0] << std::endl;
  const DofIteratorType dend = solution.dend();
  for( DofIteratorType it = solution.dbegin(); it != dend; ++it )
    *it -= average[0];
}

#endif // #ifndef RHS_HH
