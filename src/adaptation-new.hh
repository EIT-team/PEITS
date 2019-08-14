#ifndef ADAPTATION_HH
#define ADAPTATION_HH

/** include the grid capabilities
 ** to distiguish grids without local adaptation **/
#include <dune/common/version.hh>
#include <dune/common/timer.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/utility/persistentcontainer.hh>

// global counter of adaptation cyclces
static int adaptationSequenceNumber = 0;

#include "datamap.hh"
#include "loadbalance_zoltan.hh"

// interface class for callback adaptation
#include <dune/grid/common/adaptcallback.hh>

typedef struct
{
  int changes; // 1 if partitioning was changed, 0 otherwise
  int numGidEntries;  // Number of integers used for a global ID
  int numExport;      // Number of vertices I must send to other processes
  unsigned int *exportGlobalGids; // Global IDs of the vertices I must send
  int *exportProcs;    // Process to which I send each of the vertices
  
} ZOLTAN_PARTITIONING;

// LeafAdaptation
// --------------

/** \class LeafAdaptation
 *  \brief class used the adaptation procedure.
 *
 *  \tparam Grid     is the type of the underlying grid
 *  \tparam Vector   is the type of the solution vector
 */
template< class Grid, class Vector, class LoadBalanceHandle >
class LeafAdaptation : public Dune::AdaptDataHandle< Grid, LeafAdaptation< Grid, Vector, LoadBalanceHandle > >
{
  typedef LeafAdaptation<Grid,Vector,LoadBalanceHandle> ThisType;
public:
  // dimensions of grid and world
  static const int dimGrid = Grid::dimension;
  static const int dimWorld = Grid::dimensionworld;

  // type used for coordinates in the grid
  typedef typename Grid::ctype ctype;

  // types of grid's level and hierarchic iterator
  static const Dune::PartitionIteratorType partition = Dune::Interior_Partition;
  typedef typename Grid::template Codim< 0 >::template Partition< partition >::LevelIterator LevelIterator;
  typedef typename Grid::Traits::HierarchicIterator HierarchicIterator;

  // types of entity, entity pointer and geometry
  typedef typename Grid::template Codim< 0 >::Entity Entity;

#ifdef USE_VECTOR_FOR_PWF
  // container to keep data save during adaptation and load balancing
  typedef Dune::PersistentContainer<Grid,typename Vector::LocalDofVector> Container;
#else
  typedef typename Vector::VectorType Container;
#endif

  // type of grid view used
  typedef typename Vector :: GridView  GridView;

  typedef typename GridView
      ::template Codim< 0 >::template Partition< partition >::Iterator
      Iterator;
public:
  /** \brief constructor
   *  \param grid   the grid to be adapted
   */
  LeafAdaptation ( Grid &grid, LoadBalanceHandle &ldb )
  : grid_( grid ),
    ldb_( ldb ),
#ifdef USE_VECTOR_FOR_PWF
    // create persistent container for codimension 0
    container_( grid_, 0 ),     // in this version we need to provide extra storage for the dofs
#endif
    solution_( 0 ),
    adaptTimer_(),
    adaptTime_( 0.0 ),
    restProlTime_( 0.0 ),
    lbTime_( 0.0 ),
    commTime_( 0.0 )
  {}

  /** \brief main method performing the adaptation and
             perserving the data.
      \param solution  the data vector to perserve during
                       adaptation. This class must conform with the
                       parameter class V in the DataMap class and additional
                       provide a resize and communicate method.
  **/
  void operator() ( Vector &solution );

  //! return time spent for the last adapation in sec
  double adaptationTime() const { return adaptTime_; }
  //! return time spent for the last adapation in sec
  double restProlTime() const { return restProlTime_; }
  //! return time spent for the last load balancing in sec
  double loadBalanceTime() const { return lbTime_; }
  //! return time spent for the last communication in sec
  double communicationTime() const { return commTime_; }

  // this is called before the adaptation process starts
  void initialize ();

  // this is called after the adaptation process is finished
  void finalize ();

  //--------------------------------------------------
  //  Interface methods for callback adaptation
  //--------------------------------------------------
  // called when children of father are going to vanish
  void preCoarsening ( const Entity &father )
  {
#ifndef USE_VECTOR_FOR_PWF
    Container &container_ = getSolution().container();
#endif
    Vector::restrictLocal( father, container_ );
  }

  // called when children of father where newly created
  void postRefinement ( const Entity &father )
  {
#ifndef USE_VECTOR_FOR_PWF
    Container &container_ = getSolution().container();
#endif
    container_.resize();
    Vector::prolongLocal( father, container_ );
  }

private:
  Vector& getSolution()             { assert( solution_ ); return *solution_; }
  const Vector& getSolution() const { assert( solution_ ); return *solution_; }

  Grid&              grid_;
  LoadBalanceHandle& ldb_;
#ifdef USE_VECTOR_FOR_PWF
  mutable Container  container_;
#endif
  Vector*            solution_;

  Dune :: Timer      adaptTimer_ ;

  double adaptTime_;
  double restProlTime_;
  double lbTime_;
  double commTime_;
};

template< class Grid, class Vector, class LoadBalanceHandle >
inline void LeafAdaptation< Grid, Vector,LoadBalanceHandle >::operator() ( Vector &solution )
{
  if (Dune :: Capabilities :: isCartesian<Grid> :: v)
    return;

  // set pointer to solution
  solution_ = & solution ;

  adaptTime_ = 0.0;
  lbTime_    = 0.0;
  commTime_  = 0.0;
  restProlTime_ = 0.0;

  // reset timer
  adaptTimer_.reset() ;

  // copy solution to PersistentContainer if necessary
  initialize();

  // callback adaptation, see interface methods above
  grid_.adapt( *this );

  // copy solution from PersistentContainer if necessary
  finalize();

  // increase adaptation secuence number
  ++adaptationSequenceNumber;
}

template< class Grid, class Vector, class LoadBalanceHandle >
inline void LeafAdaptation< Grid, Vector, LoadBalanceHandle >
  ::initialize()
{
#ifdef USE_VECTOR_FOR_PWF
  const Vector& solution = getSolution();
  const GridView &gridView = solution.gridView();

  // first store all leave data in container
  const Iterator end = gridView.template end  < 0, partition >();
  for(  Iterator it  = gridView.template begin< 0, partition >(); it != end; ++it )
  {
    const Entity &entity = *it;
    solution.getLocalDofVector( entity, container_[ entity ] );
  }
#endif
}

template< class Grid, class Vector, class LoadBalanceHandle >
inline void LeafAdaptation< Grid, Vector, LoadBalanceHandle > ::finalize()
{
  Vector& solution = getSolution();
#ifndef USE_VECTOR_FOR_PWF
  Container &container_ = solution.container();
#endif

  adaptTime_ = adaptTimer_.elapsed();

  Dune :: Timer lbTimer ;
  DataHandle<Grid,Container> dataHandle( grid_, container_ ) ;
#if USE_ZOLTANLB || USE_SIMPLELB
  if ( ldb_.repartition() )
  {
    grid_.repartition( ldb_, dataHandle );
  }
#else
  grid_.loadBalance( ldb_, dataHandle );
#endif
  lbTime_ = lbTimer.elapsed();

  // reduce size of container, if possible
  container_.resize();

#ifdef USE_VECTOR_FOR_PWF
  // reset timer to count again
  adaptTimer_.reset();

  // resize to current grid size
  solution.resize();

  // retrieve data from container and store on new leaf grid
  const GridView &gridView = solution.gridView();
  const Iterator end = gridView.template end  < 0, partition >();
  for(  Iterator it  = gridView.template begin< 0, partition >(); it != end; ++it )
  {
    const Entity &entity = *it;
    solution.setLocalDofVector( entity, container_[ entity ] );
  }

  // store adaptation time
  adaptTime_ += adaptTimer_.elapsed();
#endif

  Dune::Timer commTimer ;
  // copy data to ghost entities
  solution.communicate();
  commTime_ = commTimer.elapsed();

  // reset pointer
  solution_ = 0;
}

#endif
