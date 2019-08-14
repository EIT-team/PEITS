#ifndef DATAMAP_HH
#define DATAMAP_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

// DataMap::LoadBalanceHandle
// --------------------------
 /** \brief the communication data handle for load balancing
 */
template< class Grid, class Container >
class DataHandle
: public Dune::CommDataHandleIF< DataHandle<Grid,Container>, Container >
{
  typedef DataHandle This;
  typedef Dune::CommDataHandleIF< This, Container > Base;


public:
  // type of data transported in the stream
  typedef typename Base::DataType DataType;
  typedef typename Container::Grid::template Codim< 0 >::Entity Entity;

protected:
  // data map
  Container &data_;
  const Grid &grid_;

public:
  //! create DiscreteOperator with a LocalOperator
  DataHandle ( const Grid &grid, Container &data )
  : data_( data )
  , grid_(grid)
  {}

  //! see documentation in Dune::CommDataHandleIF
  bool contains ( int dim, int codim ) const
  {
    return (codim == 0);
  }

  //! see documentation in Dune::CommDataHandleIF
  bool fixedsize ( int dim, int codim ) const
  {
    return (codim > 0);
  }

  //! see documentation in Dune::CommDataHandleIF (note that isLeaf available only on codimension 0 entity
  size_t size ( const Entity &entity ) const
  {
    return entity.isLeaf() ? 1 : 0 ;
  }

  //! see documentation in Dune::CommDataHandleIF (method for codim 0 entities)
  template< class Buffer >
  void gather ( Buffer &buffer, const Entity &entity ) const
  {
    // we only have data on the leaf level
    if( entity.isLeaf() )
    {
      // write data to stream
      buffer.write( data_[ entity ] );
    }
  }

  //! see documentation in Dune::CommDataHandleIF (method for codim 0 entities)
  template< class Buffer >
  void scatter ( Buffer &buffer, const Entity &entity, size_t n )
  {
    assert( n == size( entity ) );

    // we only have data on the leaf level
    if( entity.isLeaf() )
    {
      data_.resize();
      buffer.read( data_[ entity ] );
    }
  }

  //! see documentation in Dune::CommDataHandleIF (method for general entities)
  template <class E>
  size_t size ( const E &entity ) const
  {
    return 0;
  }

  //! see documentation in Dune::CommDataHandleIF (method for general entities)
  template< class Buffer, class E >
  void gather ( Buffer &buffer, const E &entity ) const
  {}

  //! see documentation in Dune::CommDataHandleIF (method for general entities)
  template< class Buffer, class E >
  void scatter ( Buffer &buffer, E &entity, size_t n )
  {
    assert( n == size( entity ) );
  }
};

#endif // #ifndef DATAMAP_HH
