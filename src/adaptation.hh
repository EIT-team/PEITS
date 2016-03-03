/******************************************************************************

  ALUGrid - a library providing a mesh manager supporting simplicial
  and hexahedral meshes, local grid adaptivity for use in parallel
  computations including dynamic load balancing.

  Copyright (C) 1998 - 2002 Bernhard Schupp
  Copyright (C) 1998 - 2002 Mario Ohlberger
  Copyright (C) 2004 - 2012 Robert Kloefkorn
  Copyright (C) 2005 - 2012 Andreas Dedner
  Copyright (C) 2010 - 2012 Martin Nolte

  The DUNE ALUGrid module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The ALUGrid library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

******************************************************************************/

#ifndef ADAPTATION_HH 
#define ADAPTATION_HH 

/** include the grid capabilities 
 ** to distiguish grids without local adaptation **/
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include "datamap.hh"
#include "loadbalance.hh"

#include <time.h>


// GridMarker
// ----------

/** \class GridMarker
 *  \brief class for marking entities for adaptation.
 *
 *  This class provides some additional strategies for marking elements
 *  which are not so much based on an indicator but more on mere
 *  geometrical and run time considerations. If based on some indicator
 *  an entity is to be marked, this class additionally tests for example
 *  that a maximal or minimal level will not be exceeded. 
 */
template< class Grid >
struct GridMarker 
{
  typedef typename Grid::template Codim< 0 >::Entity Entity;

  /** \brief constructor
   *  \param grid     the grid. Here we can not use a grid view since they only
   *                  a constant reference to the grid and the
   *                  mark method can not be called.
   *  \param minLevel the minimum refinement level
   *  \param maxLevel the maximum refinement level
   */
  GridMarker( Grid &grid, int minLevel, int maxLevel )
  : grid_(grid),
    minLevel_( minLevel ),
    maxLevel_( maxLevel ),
    wasMarked_( 0 )
  {}

  /** \brief mark an element for refinement 
   *  \param entity  the entity to mark; it will only be marked if its level is below maxLevel.
   */
  void refine ( const Entity &entity )
  {
    if( entity.level() < maxLevel_ )
    {
      grid_.mark( 1, entity );
      wasMarked_ = 1;
    }
  }

  /** \brief mark all neighbors of a given entity for refinement 
   *  \param gridView the grid view from which to take the intersection iterator 
   *  \param entity the corresponding entity
   */
  template< class GridView >
  void refineNeighbors ( const GridView &gridView, const Entity &entity )
  {
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    const IntersectionIterator end = gridView.iend( entity );
    for( IntersectionIterator it = gridView.ibegin( entity ); it != end; ++it )
    {
      const Intersection &intersection = *it;
      if( intersection.neighbor() )
        refine( *intersection.outside() );
    }
  }

  /** \brief mark an element for coarsening
   *  \param entity  the entity to mark; it will only be marked if its level is above minLevel.
   */
  void coarsen ( const Entity &entity )
  {
    if( (get( entity ) <= 0) && (entity.level() > minLevel_) )
    {
      grid_.mark( -1, entity );
      wasMarked_ = 1;
    }
  }

  /** \brief get the refinement marker 
   *  \param entity entity for which the marker is required
   *  \return value of the marker
   */
  int get ( const Entity &entity ) const
  {
    return grid_.getMark( entity );
  }

  /** \brief returns true if any entity was marked for refinement 
   */
  bool marked() 
  {
    wasMarked_ = grid_.comm().max (wasMarked_);
    return (wasMarked_ != 0);
  }

private:
  Grid &grid_;
  int minLevel_;
  int maxLevel_;
  int wasMarked_;
};

// LeafAdaptation
// --------------

/** \class LeafAdaptation
 *  \brief class used the adaptation procedure.
 *
 *  \tparam Grid     is the type of the underlying grid
 */
template< class Grid>
struct LeafAdaptation
{
  // dimensions of grid and world
  static const int dimGrid = Grid::dimension;
  static const int dimWorld = Grid::dimensionworld;

  // type used for coordinates in the grid
  typedef typename Grid::ctype ctype;

  // type of id set
  typedef typename Grid::Traits::LocalIdSet IdSet;

  // types of grid's level and hierarchic iterator
  static const Dune::PartitionIteratorType partition = Dune::Interior_Partition;
  typedef typename Grid::template Codim< 0 >::template Partition< partition >::LevelIterator LevelIterator;
  typedef typename Grid::Traits::HierarchicIterator HierarchicIterator;

  // types of entity, entity pointer and geometry
  typedef typename Grid::template Codim< 0 >::Entity Entity;

public:
  /** \brief constructor
   *  \param grid   the grid to be adapted
   */
  LeafAdaptation ( Grid &grid, ZOLTAN_PARTITIONING new_partitioning )
  : grid_( grid ),
	new_partitioning_( new_partitioning ),
    adaptTime_( 0.0 ),
    lbTime_( 0.0 ),
    commTime_( 0.0 )
  {}

  /** \brief main method performing the adaptation and
             perserving the data.
      \param solution the data vector to perserve during 
             adaptation. This class must conform with the
             parameter class V in the DataMap class and additional
             provide a resize and communicate method.
  **/
  template< class Vector >
  void operator() ( Vector &solution );

  //! return time spent for the last adapation in sec 
  double adaptationTime() const { return adaptTime_; }
  //! return time spent for the last load balancing in sec
  double loadBalanceTime() const { return lbTime_; }
  //! return time spent for the last communication in sec
  double communicationTime() const { return commTime_; }

private:
  /** \brief do restriction of data on leafs which might vanish
   *         in the grid hierarchy below a given entity
   *  \param entity   the entity from where to start the restriction process
   *  \param dataMap  map containing the leaf data to be used to store the
   *                  restricted data.
   **/
  template< class Vector, class DataMap >
  void hierarchicRestrict ( const Entity &entity, DataMap &dataMap ) const;

  /** \brief do prolongation of data to new elements below the given entity
   *  \param entity  the entity from where to start the prolongation
   *  \param dataMap the map containing the data and used to store the
   *                 data prolongt to the new elements
   **/
  template< class Vector, class DataMap >
  void hierarchicProlong ( const Entity &entity, DataMap &dataMap ) const;

  Grid &grid_;
  ZOLTAN_PARTITIONING new_partitioning_;

  double adaptTime_;
  double lbTime_;
  double commTime_;
};

template< class Grid >
template< class Vector >
inline void LeafAdaptation< Grid >::operator() ( Vector &conductivity )
{
  if (Dune :: Capabilities :: isCartesian<Grid> :: v)
    return;

  adaptTime_ = 0.0;
  lbTime_    = 0.0;
  commTime_  = 0.0;
  Dune :: Timer adaptTimer ; 

  // copy complete solution vector to map
  typedef typename Vector::GridView GridView;
  typedef typename GridView
    ::template Codim< 0 >::template Partition< partition >::Iterator 
    Iterator;
  const GridView &gridView = conductivity.gridView();

  // container to keep data save during adaptation and load balancing
  typedef Dune::PersistentContainer<Grid,typename Vector::LocalDofVector> Container;
  Container container(grid_,0);
  const Container & ccontainer = container;

  // first store all leaf data in container
  {
    const Iterator &end = gridView.template end< 0, partition >();
    for( Iterator it = gridView.template begin< 0, partition >(); it != end; ++it )
    {
      const Entity &entity = *it;
      conductivity.getLocalDofVector( entity, container[ entity ] );
    }
  }

  // check if elements might be removed in next adaptation cycle 
  const bool mightCoarsen = grid_.preAdapt();

  // if elements might be removed
  if( mightCoarsen )
  {
    // restrict data and save leaf level 
    const LevelIterator end = grid_.template lend< 0, partition >( 0 );
    for( LevelIterator it = grid_.template lbegin< 0, partition >( 0 ); it != end; ++it )
	{
      hierarchicRestrict<Vector>( *it, container );
	}
  }

  // adapt grid, returns true if new elements were created 
  const bool refined = grid_.adapt();

  // interpolate all new cells to leaf level 
  if( refined )
  {
    container.resize();
    const LevelIterator end = grid_.template lend< 0, partition >( 0 );
    for( LevelIterator it = grid_.template lbegin< 0, partition >( 0 ); it != end; ++it )
	{
      hierarchicProlong<Vector>( *it, container );
	}
  }

  adaptTime_ = adaptTimer.elapsed();

  Dune :: Timer lbTimer ;
  // re-balance grid 
  DataHandle<Grid,Container> dataHandle( grid_, container );
  LoadBalanceHandle<Grid> ldb( grid_, new_partitioning_ );

  timeval secondstart,secondstop;
	if( Dune::Fem::MPIManager::rank() == 0 )
  {
	  std::cout << "    Second loadbalance... ";
    gettimeofday(&secondstart, NULL);
  }

  grid_.loadBalance( ldb, dataHandle );

	if( Dune::Fem::MPIManager::rank() == 0 )
  {
    gettimeofday(&secondstop, NULL);
      double seconddiff = secondstop.tv_sec-secondstart.tv_sec + (secondstop.tv_usec-secondstart.tv_usec)/1000000.0;
	  std::cout << "done!                   It took " << seconddiff << " seconds" << endl;
  }

  // cleanup adaptation markers 
  grid_.postAdapt();

  // resize solution vector if elements might have been removed 
  // or were created 
  if( refined || mightCoarsen )
  {
    conductivity.resize();
  }

  // retrieve data from container and store on new leaf grid
  {
    const Iterator &end = gridView.template end< 0, partition >();
    for( Iterator it = gridView.template begin< 0, partition >(); it != end; ++it )
    {
      const Entity &entity = *it;

      typename Vector::LocalDofVector dat( ccontainer[ entity ] );
      conductivity.setLocalDofVector( entity, dat );
    }
  }
  lbTime_ = lbTimer.elapsed();

  Dune::Timer commTimer ;
  // copy data to ghost entities
  conductivity.communicate();
  commTime_ = commTimer.elapsed();
}

template< class Grid >
template< class Vector, class DataMap >
inline void LeafAdaptation< Grid >
  ::hierarchicRestrict ( const Entity &entity, DataMap &dataMap ) const
{
  // for leaf entities just copy the data to the data map
  if( !entity.isLeaf() )
  {
    // check all children first 
    bool doRestrict = true;
    const int childLevel = entity.level() + 1;
    const HierarchicIterator hend = entity.hend( childLevel );
    for( HierarchicIterator hit = entity.hbegin( childLevel ); hit != hend; ++hit )
    {
      const Entity &child = *hit;
      hierarchicRestrict<Vector>( child, dataMap );
      doRestrict &= child.mightVanish();
    }

    // if there is a child that does not vanish, this entity may not vanish
    assert( doRestrict || !entity.mightVanish() );

    // if( doRestrict )
      Vector::restrictLocal( entity, dataMap );
  }
}

template< class Grid >
template< class Vector, class DataMap >
inline void LeafAdaptation< Grid >
  ::hierarchicProlong ( const Entity &entity, DataMap &dataMap ) const
{
  if ( !entity.isLeaf() )
  {
    const int childLevel = entity.level() + 1;
    const HierarchicIterator hend = entity.hend( childLevel );
    HierarchicIterator hit = entity.hbegin( childLevel );

    const bool doProlong = hit->isNew();
    if( doProlong )
      Vector::prolongLocal( entity, dataMap );

    // if the children have children then we have to go deeper 
    for( ; hit != hend; ++hit ) 
    {
      assert(doProlong == hit->isNew());
      hierarchicProlong<Vector>( *hit, dataMap );
    }
  }
}

#endif 
