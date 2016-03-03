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

#ifndef PIECEWISEFUNCTION_HH
#define PIECEWISEFUNCTION_HH

#include <vector>
#include <fstream>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/common/datahandleif.hh>

// PiecewiseFunction
// ----------
/** \class PiecewiseFunction
 *  \brief a piecewise constant function
 *
 *  \tparam  View   grid view on which the function is defined
 *  \tparam  Range  type for the range vector of the function
 *                  This class must have a vector like structure, i.e.,
 *                  a double operator[](int) method.
 *                  Furthermore we expect an axpy method on the class.
 */
template< class View, class Range >
class PiecewiseFunction      
{ 
  typedef PiecewiseFunction< View, Range > This;

public:
  typedef View GridView;
  typedef Range RangeType;

  static const unsigned int order = 0;

  typedef std::vector< RangeType > VectorType; 

  typedef typename GridView::template Codim< 0 >::Entity Entity;
  typedef typename Entity::Geometry Geometry;
  typedef typename Geometry::LocalCoordinate DomainType;
  typedef typename Geometry::GlobalCoordinate GlobalType;

  /* vector of all dofs on an codim 0 entity */
  typedef RangeType LocalDofVector;

private:
  struct CommDataHandle;

public:
  /** \brief constructor
   *
   *  \param[in]  gridView  grid view the function lives on
   */
  PiecewiseFunction( const GridView &gridView )
  : gridView_( gridView ),
    dof_( gridView.indexSet().size( 0 ) )
  {}

  const RangeType &operator[] ( const size_t &index ) const
  {
    return dof_[ index ];
  }

  RangeType &operator[] ( const size_t &index )
  {
    return dof_[ index ];
  }

  /** \brief access the DoFs on one entity (of codimension 0)
   *
   *  \param[in]  entity  entity whose DoFs to access
   *
   *  \returns a reference to the DoFs of the entity
   */
  const RangeType &operator[] ( const Entity &entity ) const;

  /** \brief access the DoFs on one entity (of codimension 0)
   *
   *  \param[in]  entity  entity whose DoFs to access
   *
   *  \returns a reference to the DoFs of the entity
   */
  RangeType &operator[] ( const Entity &entity );

  /** \brief evaluate the discrete function at an arbitrary point on the
   *         reference element 
   *  \param entity the entity on which the function is to be evaluated
   *  \param x      point on entity in local coordinate system
   */
  RangeType evaluate ( const Entity &entity, const DomainType &x ) const
  {
    return (*this)[ entity ];
  }

  /** \brief obtain the size of the dof vector;
             might be larger that number of elements in grid view. */
  size_t size() const 
  {
    return dof_.size();
  }

  /** \brief make a copy of the of he vector **/
  void assign( const This &other )
  {
    this->dof_ = other.dof_;
  }

  /** \brief obtain the grid view, this function lives on */
  const GridView &gridView () const
  {
    return gridView_;
  }

  /** \brief method to initialize data
   *  \param problemData instance of a class having a method initial taking
   *         points \c x in the global coordinate system
   */
  template <class ProblemData>
  void initialize ( const ProblemData &problemData );

  /** \brief add a multiple of another function to this one
   *
   *  The classic BLAS level 1 function. In pseudo code it does the following:
   *  \code
   *  (*this) += lambda * other
   *  \endcode
   *
   *  \param[in]  lambda  factor to multiply the other function by
   *  \param[in]  other   function to add to this one
   */
  void axpy ( const double &lambda, const This &other );

  /** \brief add a multiple of another function to this one
   *
   *  \code
   *  (*this) = mu * (*this) + lambda * other
   *  \endcode
   *
   *  \param[in]  mu      factor to multiply the this function by
   *  \param[in]  lambda  factor to multiply the other function by
   *  \param[in]  other   function to add to this one
   */
  void addAndScale ( const double &mu, const double &lambda, const This &other );

  /** \brief set this function to zero */
  void clear ();

  /** \brief adapt the DoF vector to the size of the grid view */
  void resize();

  /** \brief extract all DoFs on an entity 
   *  \param[in] entity the entity
   *  \param[in] localDofs insert the DoFs into this storage vector
   */
  void getLocalDofVector ( const Entity &entity, LocalDofVector &localDofs ) const;
  /** \brief copy all DoFs onto an entity from a local storage vector
   *  \param[in] entity the entity
   *  \param[in] localDofs extract the DoFs from here
   */
  void setLocalDofVector ( const Entity &entity, const LocalDofVector &localDofs );

  /** \brief method to restrict the data from the children of an entity
   *  \param[in] father    the father entity
   *  \param[in] localDofs the local storage vector where the restricted
   *                       data is to be stored
   *                       (both the father data storage space and
   *                       the data on the children can be accessed by
   *                       an \ operator[](const Entity&) method)
   */
  template< class LocalDofMap >
  static void restrictLocal ( const Entity &father, LocalDofMap &localDofs );
  /** \brief method to prolong the data to all children entities from a
   *  \param[in] father    the father entity
   *  \param[in] localDofs the storage vector where the prolonged
   *                       data is to be stored 
   *                       (both the father data and
   *                       the child data stroage space can be accessed by
   *                       an \ operator[](const Entity&) method)
   */
  template< class LocalDofMap >
  static void prolongLocal ( const Entity &father, LocalDofMap &localDofs );

  /** \brief copy interior cell data to overlap of ghost cells */
  void communicate ();

private:  
  GridView gridView_;
  /* storage for dofs */
  VectorType dof_;
};

template< class View, class Range >
inline const typename PiecewiseFunction< View, Range >::RangeType &
PiecewiseFunction< View, Range >::operator[] ( const Entity &entity ) const
{
  assert( gridView().indexSet().contains( entity ) );
  return dof_[ gridView().indexSet().index( entity ) ];
}

template< class View, class Range >
inline typename PiecewiseFunction< View, Range >::RangeType &
PiecewiseFunction< View, Range >::operator[] ( const Entity &entity )
{
  assert( gridView().indexSet().contains( entity ) );
  return dof_[ gridView().indexSet().index( entity ) ];
}

template< class View, class Range >
template< class ProblemData >
inline void
PiecewiseFunction< View, Range >::initialize ( const ProblemData &problemData )
{
  // resize data vector and set all data to zero
  resize();

  /* first we extract the dimensions of the grid 
  static const int dim = GridView::dimension;
  static const int dimworld = GridView::dimensionworld;
  */

  /* get type of iterator over leaf entities of codimension zero */
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;

  /* types of entity and geometry */
  typedef typename GridView::template Codim< 0 >::Entity Entity;
  typedef typename GridView::template Codim< 0 >::Geometry Geometry;

  /* types of vectors */
  typedef typename Geometry::LocalCoordinate DomainType;
  typedef typename Geometry::GlobalCoordinate GlobalType;

// MJE
  typedef typename GridView::IndexSet IndexSetType;
  const IndexSetType &indexSet = gridView_.indexSet();
// end MJE

  // loop over all entities
  const Iterator end = gridView().template end< 0 >();
  for( Iterator it = gridView().template begin< 0 >(); it != end; ++it )
  {   
    const Entity &entity = *it;
    (*this)[ entity ] = problemData[indexSet.index(*it)];//problemData.initial( baryCenter );  // MJE
  }

  // copy data to ghosts
  communicate();
}

template< class View, class Range >
inline void
PiecewiseFunction< View, Range >::axpy ( const double &lambda, const This &other )
{
  const size_t size = dof_.size();
  for( size_t i = 0; i < size; ++i )
    (*this)[ i ].axpy( lambda, other[ i ] );
}

template< class View, class Range >
inline void
PiecewiseFunction< View, Range >::addAndScale ( const double &mu, const double &lambda, const This &other )
{
  const size_t size = dof_.size();
  for( size_t i = 0; i < size; ++i )
  {
    (*this)[ i ] *= mu;
    (*this)[ i ].axpy( lambda, other[ i ] );
  }
}

template< class View, class Range >
inline void PiecewiseFunction< View, Range >::clear ()
{
  const typename VectorType::iterator end = dof_.end();
  for( typename VectorType::iterator it = dof_.begin(); it != end; ++it )
    *it = RangeType( 0 );
}

template< class View, class Range >
inline void PiecewiseFunction< View, Range >::resize ()
{
  const size_t newSize = gridView().indexSet().size( 0 );
  if( dof_.capacity() < newSize ) 
  {
    //const size_t reserveSize = ((size_t) ((double) 1.1 * newSize));
    const size_t reserveSize = newSize * (1024 + 128) / 1024;
    // reserve some more  memory to avoid to much re-allocation 
    // (possibly not required for std::vector)
    dof_.reserve( reserveSize );
  }
  // resize to current size 
  dof_.resize( newSize );
}

template< class View, class Range >
inline void PiecewiseFunction< View, Range >
  ::getLocalDofVector ( const Entity &entity, LocalDofVector &localDofs ) const
{
  localDofs = (*this)[ entity ];
}

template< class View, class Range >
inline void PiecewiseFunction< View, Range >
  ::setLocalDofVector ( const Entity &entity, const LocalDofVector &localDofs )
{
  (*this)[ entity ] = localDofs;
}

template< class View, class Range >
template< class LocalDofMap >
inline void PiecewiseFunction< View, Range >
  ::restrictLocal ( const Entity &father, LocalDofMap &localDofs )
{
  const LocalDofMap &clocalDofs = localDofs;
  typedef typename Entity::HierarchicIterator HierarchicIterator;

  LocalDofVector &fatherValue = localDofs[ father ];

  double volume = 0;

  const int childLevel = father.level() + 1;
  const HierarchicIterator hend = father.hend( childLevel );
  HierarchicIterator hit = father.hbegin( childLevel );
  const Entity &child = *hit;
  const double childVolume = child.geometry().volume();
  fatherValue = clocalDofs[ child ];
  fatherValue *= childVolume;
  volume += childVolume;
  ++hit;
  for( ; hit != hend; ++hit )
  {
    const Entity &child = *hit;
    const double childVolume = child.geometry().volume();
    fatherValue.axpy( childVolume, clocalDofs[ child ] );
    volume += childVolume;
  }

  // Note: only true if we are not on a surface
  assert( std::abs( volume - father.geometry().volume() ) < 1e-12 );
  fatherValue /= volume;
}

template< class View, class Range >
template< class LocalDofMap >
inline void PiecewiseFunction< View, Range >
  ::prolongLocal ( const Entity &father, LocalDofMap &localDofs )
{
  const LocalDofMap &clocalDofs = localDofs;
  typedef typename Entity::HierarchicIterator HierarchicIterator;

  const LocalDofVector &fatherValue = clocalDofs[ father ];

  const int childLevel = father.level() + 1;
  const HierarchicIterator hend = father.hend( childLevel );
  for( HierarchicIterator hit = father.hbegin( childLevel ); hit != hend; ++hit )
  {
    const Entity &child = *hit;
    localDofs[ child ] = fatherValue;
  }
}

template< class View, class Range >
inline void PiecewiseFunction< View, Range >::communicate ()
{
  const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface;
  const Dune::CommunicationDirection direction = Dune::ForwardCommunication;

  CommDataHandle handle( *this );
  gridView().communicate( handle, interface, direction );
}

// PiecewiseFunction::CommDataHandle
// ----------------------
/** \brief the communication data handle for the PiecewiseFunction class
 */
template< class View, class Range >
struct PiecewiseFunction< View, Range >::CommDataHandle
: public Dune::CommDataHandleIF< CommDataHandle, Range >
{
  // type of communicated data (i.e. double) 
  typedef RangeType DataType;

  //! constructor taking a PiecewiseFunction instance
  CommDataHandle ( PiecewiseFunction< View, Range > &data )
  : data_( data )
  {}

  //! see documentation in Dune::CommDataHandleIF 
  bool contains ( int dim, int codim ) const
  {
    return (codim == 0);
  }
  
  //! see documentation in Dune::CommDataHandleIF 
  bool fixedsize ( int dim, int codim ) const
  {
    return true;
  }
  
  //! see documentation in Dune::CommDataHandleIF 
  template< class E >
  size_t size ( const E &e ) const
  {
    return (Entity::codimension == 0 ? 1 : 0); 
  }
  
  //! see documentation in Dune::CommDataHandleIF (method for codim 0 entities)
  template< class Buffer >
  void gather ( Buffer &buffer, const Entity &entity ) const
  { 
    buffer.write( data_[ entity ] );
  }

  //! see documentation in Dune::CommDataHandleIF (method for codim 0 entities)
  template< class Buffer >
  void scatter ( Buffer &buffer, const Entity &entity, size_t n )
  { 
    // communication size is fixed, so we can check the sizes
    assert( n == size( entity ) );
    
    // read data from buffer 
    buffer.read( data_[ entity ] );
  }

  //! see documentation in Dune::CommDataHandleIF (method for general entities)
  template< class Buffer, class E >
  void gather ( Buffer &buffer, const E &entity ) const
  {}

  //! see documentation in Dune::CommDataHandleIF (method for general entities)
  template< class Buffer, class E >
  void scatter ( Buffer &buffer, const E &entity, size_t n)
  {
    // communication size is fixed, so we can check the sizes
    assert( n == size( entity ) );
  }
private:
  // data vector 
  PiecewiseFunction< View, Range > &data_;
};

template< class Data >
struct VTKData;

template <class GridViewType>
class PartitioningData
  : public Dune::VTKFunction< GridViewType >
{
  typedef PartitioningData   ThisType;

public:
  typedef typename GridViewType :: template Codim< 0 >::Entity EntityType;
  typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

  //! constructor taking discrete function 
  PartitioningData( const int rank ) : rank_( rank ) {}

  //! virtual destructor
  virtual ~PartitioningData () {}

  //! return number of components
  virtual int ncomps () const { return 1; }

  //! evaluate single component comp in
  //! the entity
  virtual double evaluate ( int comp, const EntityType &e, const LocalCoordinateType &xi ) const
  {
    return double( rank_ );
  }

  //! get name
  virtual std::string name () const
  {
    return std::string( "rank" );
  }

private:
  const int rank_;
};



/**
 * \brief a class for vtk output of a PiecewiseFunction instance
 */
template< class GridView, int dimRange >
struct VTKData< PiecewiseFunction< GridView, Dune::FieldVector<double,dimRange> > >
: public Dune::VTKWriter< GridView >::VTKFunction
{
  typedef PiecewiseFunction< GridView, Dune::FieldVector<double,dimRange> > Data;
  typedef VTKData< Data > This;

  static const int dim = GridView::dimension;
  typedef typename GridView::template Codim< 0 >::Entity Entity;
  typedef Dune::FieldVector< double, dim > DomainType;

  //! number of components for scalar use 1 and for vector 2 or 3
  int ncomps () const
  {
    return 1;
  }

  //! evaluate function (comp<ncomps) on entity for local coordinate xi
  double evaluate ( int comp, const Entity &e, const DomainType &xi ) const
  {
    int index = data_.gridView().indexSet().index(e);
    return data_[ index ][ comp_ ];
  }

  //! name for this function
  std::string name () const 
  {
    std::stringstream ret;
    ret << name_ << "_" << comp_;
    return ret.str();
  }

  //! add a PiecewiseFunction to the VTKWriter 
  static void addTo ( const Data &data, Dune::VTKWriter< GridView > &vtkWriter )
  {
    /* the vtk-Writer class takes ownership of the function added - VTKData
     * is merely a wrapper for the Data class */
    for( int i = 0; i < dimRange; ++i )
      vtkWriter.addCellData( new This( data, i, "data" ) );
  }
  //! add a PiecewiseFunction to the VTKWriter 
  static void addTo ( const Data &data, const std::string &name, Dune::VTKWriter< GridView > &vtkWriter )
  {
    /* the vtk-Writer class takes ownership of the function added - VTKData
     * is merely a wrapper for the Data class */
    for( int i = 0; i < dimRange; ++i )
      vtkWriter.addCellData( new This( data, i,name ) );
  }
  //! add rank function for visualization of the partitioning 
  static void addPartitioningData( const int rank, Dune::VTKWriter< GridView > &vtkWriter ) 
  {
    vtkWriter.addCellData( new PartitioningData< GridView >(rank) );
  }

private:
  VTKData ( const Data &data, unsigned int comp, const std::string &name )
  : data_( data ),
    comp_( comp ),
    name_(name)
  {}
  const Data &data_;
  unsigned int comp_;
  const std::string name_;
};

#endif // #ifndef PIECEWISEFUNCTION_HH
