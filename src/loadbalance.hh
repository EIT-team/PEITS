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

#ifndef LOADBALANCE_HH
#define LOADBALANCE_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

typedef struct{

  int changes; // 1 if partitioning was changed, 0 otherwise 
  int numGidEntries;  // Number of integers used for a global ID 
//  int numLidEntries;  // Number of integers used for a local ID 
//  int numImport;      // Number of vertices to be sent to me MJE: SHOULD NOT NEED THAT
//  int *importGlobalGids;  // Global IDs of vertices to be sent to me MJE: SHOULD NOT NEED THAT
//  int *importLocalGids;   // Local IDs of vertices to be sent to me MJE: DO NOT NEED THAT
//  int *importProcs;    // Process rank for source of each incoming vertex MJE: SHOULD NOT NEED THAT
//  int *importToPart;   // New partition for each incoming vertex   MJE: DO NOT NEED THAT
  int numExport;      // Number of vertices I must send to other processes
  unsigned int *exportGlobalGids;  // Global IDs of the vertices I must send 
//  int *exportLocalGids;   // Local IDs of the vertices I must send MJE: DO NOT NEED THAT
  int *exportProcs;    // Process to which I send each of the vertices 
//  int *exportToPart;  // Partition to which each vertex will belong MJE: DO NOT NEED THAT

} ZOLTAN_PARTITIONING;



template< class Grid >
class LoadBalanceHandle
: public Dune::LoadBalanceHandleIF< LoadBalanceHandle<Grid> >
{
  typedef LoadBalanceHandle This;
  typedef Dune::LoadBalanceHandleIF< This > Base;

  typedef typename Grid::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType gIdType;

public:
  typedef typename Grid :: ObjectStreamType ObjectStream;

  static const int dimension = Grid :: dimension;

  template< int codim >
  struct Codim
  {
    typedef typename Grid :: Traits :: template Codim< codim > :: Entity Entity;
    typedef typename Grid :: Traits :: template Codim< codim > :: EntityPointer
      EntityPointer;
  };
  typedef typename Codim< 0 > :: Entity Element;

private:
  const Grid &grid_;

public:
  LoadBalanceHandle ( const Grid &grid, ZOLTAN_PARTITIONING new_partitioning)
  : grid_( grid ),
	globalIdSet_( grid.globalIdSet() ),
	new_partitioning_( new_partitioning )
  {}

  bool userDefinedPartitioning () const
  {
    return true;
  }
  // return true if user defined load balancing weights are provided
  bool userDefinedLoadWeights () const
  {
    return false;
  }

  // returns true if user defined partitioning needs to be readjusted 
  bool repartition () const 
  { 
    angle_ += 2.*M_PI/50.;
    return true;
  }
  // return load weight of given element 
  int loadWeight( const Element &element ) const 
  { 
    return 1;
  }
  // return destination (i.e. rank) where the given element should be moved to 
  // this needs the methods userDefinedPartitioning to return true
  int destination( const Element &element ) const 
  { 
	  gIdType bla = globalIdSet_.id(element);
	  std::vector<int> elementGID(4); // because we have 4 vertices
	  bla.getKey().extractKey(elementGID);

	  // add one to the GIDs, so that they match the ones from Zoltan
	  transform(elementGID.begin(), elementGID.end(), elementGID.begin(), bind2nd(std::plus<int>(), 1));

	  int p = int(this->grid_.comm().rank());

	  for (int i = 0; i<new_partitioning_.numExport; ++i)
	  {
	    if (std::equal(elementGID.begin(),elementGID.end(), &new_partitioning_.exportGlobalGids[i*new_partitioning_.numGidEntries]) )
		  p = new_partitioning_.exportProcs[i];
	  }

    return p;
  }

private:
  static double angle_;
  const GlobalIdSet &globalIdSet_;
  const ZOLTAN_PARTITIONING new_partitioning_;
};

template< class Grid >
double LoadBalanceHandle<Grid>::angle_ = 0;


#endif // #ifndef LOADBALNCE_HH
