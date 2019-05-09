#ifndef LOADBALANCE_ZOLTAN_HH
#define LOADBALANCE_ZOLTAN_HH

#include <iostream>
#include <fstream>
#include <vector>

#include <dune/common/version.hh>
#include <dune/grid/common/gridenums.hh>

#if HAVE_ZOLTAN && HAVE_MPI
#include <dune/alugrid/impl/parallel/aluzoltan.hh>

template< class Grid >
class ZoltanLoadBalanceHandle
{
  typedef ZoltanLoadBalanceHandle This;

private:
  typedef typename Grid::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType GIdType;
  static const int dimension = Grid :: dimension;
  static const int NUM_GID_ENTRIES = 1;
  template< int codim >
  struct Codim
  {
    typedef typename Grid :: Traits :: template Codim< codim > :: Entity Entity;
  };

  struct ZoltanPartitioning{
    int changes; // 1 if partitioning was changed, 0 otherwise
    int numGidEntries;  // Number of integers used for a global ID
    int numLidEntries;  // Number of integers used for a global ID
    int numExport;      // Number of vertices I must send to other processes
    int numImport;      // Number of vertices I must send to other processes
    unsigned int *importLocalGids;  // Global IDs of the vertices I must send
    unsigned int *importGlobalGids; // Global IDs of the vertices I must send
    unsigned int *exportLocalGids;  // Global IDs of the vertices I must send
    unsigned int *exportGlobalGids; // Global IDs of the vertices I must send
    int *importProcs;    // Process to which I send each of the vertices
    int *exportProcs;    // Process to which I send each of the vertices
    int *importToPart;
    int *exportToPart;
  };
  struct FixedElements {
    int fixed_entities;
    std::vector<int> fixed_GID;
    std::vector<int> fixed_Process;
    FixedElements() : fixed_GID(0), fixed_Process(0) {}
  };
  struct HGraphData {              /* Zoltan will partition vertices, while minimizing edge cuts */
    int numMyVertices;             /* number of vertices that I own initially */
    ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */
    float *vtxWEIGHT;              /* weight for each node in the graph */
    int numMyHEdges;               /* number of my hyperedges */
    int numAllNbors;               /* number of vertices in my hyperedges */
    ZOLTAN_ID_TYPE *edgeGID;       /* global ID of each of my hyperedges */
    int *nborIndex;                /* index into nborGID array of edge's vertices */
    ZOLTAN_ID_TYPE *nborGID;       /* Vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]] */
    float *nborWEIGHT;             /* weight for each graph edge */
    int *nborPROC;                 /* processor for each edge node */
    FixedElements fixed_elmts;
    HGraphData() : vtxGID(0), vtxWEIGHT(0), edgeGID(0), nborIndex(0), nborGID(0), nborWEIGHT(0), nborPROC(0) {}
    ~HGraphData() { freeMemory();}
    void freeMemory()
    {
      if (!nborWEIGHT)
        free(nborWEIGHT);
      if (!nborPROC)
        free(nborPROC);
      if (!nborGID)
        free(nborGID);
      if (!nborIndex)
        free(nborIndex);
      if (!edgeGID)
        free(edgeGID);
      if (!vtxWEIGHT)
        free(vtxWEIGHT);
      if (!vtxGID)
        free(vtxGID);
    }
  };
public:
  typedef typename Codim< 0 > :: Entity Element;
  ZoltanLoadBalanceHandle ( const Grid &grid);
  ~ZoltanLoadBalanceHandle();

  // this method is called before invoking the repartition method on the
  // grid, to check if the user defined partitioning needs to be readjusted
  // So it is not part of the interface
  bool repartition ()
  {
    int elements = grid_.size(0);
    size_t sumElements = grid_.comm().sum( elements );
    size_t minElements = grid_.comm().min( elements );
    size_t maxElements = grid_.comm().max( elements );
    double mean = sumElements / grid_.comm().size();
    const bool repartition = ((maxElements > (ldbOver_ * mean)) || (minElements < (ldbUnder_ * mean) ))
                             ? true : false ;
    if (!repartition)
      return false;
    if (!first_)
    {
      Zoltan_LB_Free_Part(&(new_partitioning_.importGlobalGids),
                   &(new_partitioning_.importLocalGids),
                   &(new_partitioning_.importProcs),
                   &(new_partitioning_.importToPart) );
      Zoltan_LB_Free_Part(&(new_partitioning_.exportGlobalGids),
                   &(new_partitioning_.exportLocalGids),
                   &(new_partitioning_.exportProcs),
                   &(new_partitioning_.exportToPart) );
    }
    generateHypergraph();
    /******************************************************************
    ** Zoltan can now partition the vertices of hypergraph.
    ** In this simple example, we assume the number of partitions is
    ** equal to the number of processes.  Process rank 0 will own
    ** partition 0, process rank 1 will own partition 1, and so on.
    ******************************************************************/
    Zoltan_LB_Partition(zz_, // input (all remaining fields are output)
          &new_partitioning_.changes,        // 1 if partitioning was changed, 0 otherwise
          &new_partitioning_.numGidEntries,  // Number of integers used for a global ID
          &new_partitioning_.numLidEntries,  // Number of integers used for a local ID
          &new_partitioning_.numImport,      // Number of vertices to be sent to me
          &new_partitioning_.importGlobalGids,  // Global IDs of vertices to be sent to me
          &new_partitioning_.importLocalGids,   // Local IDs of vertices to be sent to me
          &new_partitioning_.importProcs,    // Process rank for source of each incoming vertex
          &new_partitioning_.importToPart,   // New partition for each incoming vertex
          &new_partitioning_.numExport,      // Number of vertices I must send to other processes
          &new_partitioning_.exportGlobalGids,  // Global IDs of the vertices I must send
          &new_partitioning_.exportLocalGids,   // Local IDs of the vertices I must send
          &new_partitioning_.exportProcs,    // Process to which I send each of the vertices
          &new_partitioning_.exportToPart);  // Partition to which each vertex will belong
    first_ = false;
    return (new_partitioning_.changes == 1);
  }

  // return destination (i.e. rank) where the given element should be moved to
  int operator()( const Element &element ) const
  {
	  std::vector<int> elementGID(NUM_GID_ENTRIES);
    // GIdType id = globalIdSet_.id(element);
    // id.getKey().extractKey(elementGID);
	  elementGID[0] = grid_.macroGridView().macroId(element); //   element.impl().macroID();

    // add one to the GIDs, so that they match the ones from Zoltan
    transform(elementGID.begin(), elementGID.end(), elementGID.begin(), bind2nd(std::plus<int>(), 1));

    int p = int(grid_.comm().rank());

    for (int i = 0; i<new_partitioning_.numExport; ++i)
    {
      if (std::equal(elementGID.begin(),elementGID.end(), &new_partitioning_.exportGlobalGids[i*new_partitioning_.numGidEntries]) )
      {
        p = new_partitioning_.exportProcs[i];
        break;
      }
    }
    return p;
  }
  // This method can simply return false, in which case ALUGrid will
  // internally compute the required information through some global
  // communication. To avoid this overhead the user can provide the ranks
  // of particians from which elements will be moved to the calling partitian.
  bool importRanks( std::set<int> &ranks ) const
  {
    ranks.insert( new_partitioning_.importProcs,
                  new_partitioning_.importProcs+new_partitioning_.numImport );
    // std::ostream_iterator< int > output( cout, " " );
    // cout << "Import ranks: ";
    // std::copy( ranks.begin(), ranks.end(), output );
    return true;
  }
private:
  void generateHypergraph();

  // ZOLTAN query functions
  static int get_number_of_vertices(void *data, int *ierr);
  static void get_vertex_list(void *data, int sizeGID, int sizeLID,
              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
              int wgt_dim, float *obj_wgts, int *ierr);
  static int get_num_fixed_obj(void *data, int *ierr);
  static void get_fixed_obj_list(void *data, int num_fixed_obj,
                                 int num_gid_entries, ZOLTAN_ID_PTR fixed_gids, int *fixed_part, int *ierr);
  static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                                 int num_obj,
                                 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                 int *numEdges, int *ierr);
  static void get_edge_list(void *data, int sizeGID, int sizeLID,
                            int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int *num_edges,
                            ZOLTAN_ID_PTR nborGID, int *nborProc,
                            int wgt_dim, float *ewgts, int *ierr);

  const Grid &grid_;
  const GlobalIdSet &globalIdSet_;

  Zoltan_Struct *zz_;
  HGraphData hg_;
  ZoltanPartitioning new_partitioning_;
  bool first_;
  double ldbUnder_, ldbOver_;
  const bool fix_bnd_;
};

template< class Grid >
ZoltanLoadBalanceHandle<Grid>::
ZoltanLoadBalanceHandle(const Grid &grid)
: grid_( grid )
, globalIdSet_( grid.globalIdSet() )
, first_(true)
, ldbUnder_(0), ldbOver_(1.2)
, fix_bnd_(false)
{
  zz_ = Zoltan_Create(MPI_COMM_WORLD);

  // General parameters
  Zoltan_Set_Param(zz_, "DEBUG_LEVEL", "0");
  if (!fix_bnd_) // fixing element requires using hypergraph partitioning (which is perhaps better anyway?)
    Zoltan_Set_Param(zz_, "LB_METHOD", "GRAPH");             /* partitioning method */
  else
    Zoltan_Set_Param(zz_, "LB_METHOD", "HYPERGRAPH");        /* partitioning method */
  Zoltan_Set_Param(zz_, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
  Zoltan_Set_Param(zz_, "NUM_GID_ENTRIES", "1");      /* global IDs are 1 integers */
  Zoltan_Set_Param(zz_, "NUM_LID_ENTRIES", "1");      /* local IDs are 1 integers */
  Zoltan_Set_Param(zz_, "RETURN_LISTS", "ALL");       /* export AND import lists */
  Zoltan_Set_Param(zz_, "OBJ_WEIGHT_DIM", "1");       /* provide a weight for graph nodes */
  Zoltan_Set_Param(zz_, "EDGE_WEIGHT_DIM", "1");      /* provide a weight for graph edge */
  Zoltan_Set_Param(zz_, "GRAPH_SYM_WEIGHT","MAX");
  Zoltan_Set_Param(zz_, "GRAPH_SYMMETRIZE","NONE" );
  Zoltan_Set_Param(zz_, "PHG_EDGE_SIZE_THRESHOLD", ".25");
  Zoltan_Set_Param(zz_, "CHECK_HYPERGRAPH", "0");
  Zoltan_Set_Param(zz_, "CHECK_GRAPH", "0");

  /* PHG parameters  - see the Zoltan User's Guide for many more
  */
  Zoltan_Set_Param(zz_, "LB_APPROACH", "REPARTITION");
  Zoltan_Set_Param(zz_, "REMAP", "0");
  Zoltan_Set_Param(zz_, "IMBALANCE_TOL", "1.05" );

  /* Application defined query functions */
  Zoltan_Set_Num_Obj_Fn(zz_, get_number_of_vertices, &hg_);
  Zoltan_Set_Obj_List_Fn(zz_, get_vertex_list, &hg_);
  Zoltan_Set_Num_Edges_Multi_Fn(zz_, get_num_edges_list, &hg_);
  Zoltan_Set_Edge_List_Multi_Fn(zz_, get_edge_list, &hg_);

  /* Register fixed object callback functions */
  if (Zoltan_Set_Fn(zz_, ZOLTAN_NUM_FIXED_OBJ_FN_TYPE,
        (void (*)()) get_num_fixed_obj,
        (void *) &hg_) == ZOLTAN_FATAL) {
    return;
  }
  if (Zoltan_Set_Fn(zz_, ZOLTAN_FIXED_OBJ_LIST_FN_TYPE,
        (void (*)()) get_fixed_obj_list,
        (void *) &hg_) == ZOLTAN_FATAL) {
    return;
  }

  /******************************************************************
  ** Zoltan can now partition the vertices of hypergraph.
  ** In this simple example, we assume the number of partitions is
  ** equal to the number of processes.  Process rank 0 will own
  ** partition 0, process rank 1 will own partition 1, and so on.
  ******************************************************************/

  // read config file if avaialble
  std::ifstream in( "alugrid.cfg" );
  if( in )
  {
    in >> ldbUnder_;
    in >> ldbOver_;
  }
}
template <class Grid>
ZoltanLoadBalanceHandle<Grid>::
~ZoltanLoadBalanceHandle()
{
  if (!first_)
  {
    Zoltan_LB_Free_Part(&(new_partitioning_.importGlobalGids),
                 &(new_partitioning_.importLocalGids),
                 &(new_partitioning_.importProcs),
                 &(new_partitioning_.importToPart) );
    Zoltan_LB_Free_Part(&(new_partitioning_.exportGlobalGids),
                 &(new_partitioning_.exportLocalGids),
                 &(new_partitioning_.exportProcs),
                 &(new_partitioning_.exportToPart) );
  }
  Zoltan_Destroy(&zz_);
}

/********************************************************************************
 * This method is called in the repartition callback function.
 * It needs to setup the graph/hypergraph and call the zoltan partition function
 ********************************************************************************/
template< class Grid >
void ZoltanLoadBalanceHandle<Grid>::
generateHypergraph()
{
  hg_.freeMemory();
  // setup the hypergraph by iterating over the macro level
  // (ALU can only partition on the macro level)
  const Dune::PartitionIteratorType partition = Dune::Interior_Partition;
  typedef typename Grid::MacroGridView GridView;
  const GridView &gridView = grid_.macroGridView();
  typedef typename GridView::template Codim< 0 >::template Partition< partition >::Iterator Iterator;
  typedef typename Codim< 0 >::Entity Entity;

  typedef typename GridView::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  int tempNumMyVertices = gridView.size(0);
  ZOLTAN_ID_TYPE *tempVtxGID  = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * NUM_GID_ENTRIES * tempNumMyVertices);
  ZOLTAN_ID_TYPE *tempEdgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * NUM_GID_ENTRIES * tempNumMyVertices);
  float* tempVtxWeight = (float*)malloc(sizeof(float)*tempNumMyVertices);

  std::vector<ZOLTAN_ID_TYPE> tempNborGID(0);
  std::vector<float> tempNborWeight(0);
  std::vector<int> tempNborProc(0);
  std::vector<int> fixedProcVector(0), fixedElmtVector(0);
  int *tempNborIndex = (int *)malloc(sizeof(int) * (tempNumMyVertices + 1));
  tempNborIndex[0] = 0;

  unsigned int element_count = 0;
  const Iterator &end = gridView.template end< 0, partition >();
  for( Iterator it = gridView.template begin< 0, partition >(); it != end; ++it )
  {
	  const Entity &entity = *it;
	  std::vector<int> elementGID(NUM_GID_ENTRIES);
    // use special ALU method that returns a pure integer tuple which is a
    // unique id on the macrolevel
	  elementGID[0] = gridView.macroId(entity);

	  for (int i=0; i<NUM_GID_ENTRIES; ++i)
	  {
	    tempVtxGID[element_count*NUM_GID_ENTRIES + i] = (ZOLTAN_ID_TYPE)elementGID[i] + 1;   // global identifier of element    ADD ONE BECAUSE WE START COUNTING AT 0 AND ZOLTAN DOES NOT LIKE IT
	    tempEdgeGID[element_count*NUM_GID_ENTRIES + i] = (ZOLTAN_ID_TYPE)elementGID[i] + 1;  // global identifier of hyperedge
  	}
    // get weight associated with entity using ALU specific function
    tempVtxWeight[element_count] = gridView.weight(entity);

    // now setup the edges
    const IntersectionIterator iend = gridView.iend( entity );
	  int num_of_neighbors = 0;
    float weight = 0;
    for( IntersectionIterator iit = gridView.ibegin( entity ); iit != iend; ++iit )
    {
      const Intersection &intersection = *iit;
      if( intersection.neighbor() )
	    {
		    const Entity &neighbor = intersection.outside();

		    std::vector<int> neighborGID(NUM_GID_ENTRIES);
        // use special ALU method that returns a pure integer tuple which is a
        // unique id on the macrolevel
	      neighborGID[0] = gridView.macroId(neighbor);
        // use the alu specific weight function between neighboring elements
        weight += gridView.weight( intersection );

		    for (int i=0; i<NUM_GID_ENTRIES; ++i)
		    {
		      tempNborGID.push_back((ZOLTAN_ID_TYPE)neighborGID[i] + 1);
		    }
        tempNborWeight.push_back( gridView.weight( intersection ) );
        tempNborProc.push_back( gridView.master( neighbor ) );

		    num_of_neighbors++;
	    }
      else
      {
        // Find if element is candidate for user-defined partitioning:
        // we keep the left boundary on process 0
        if ( intersection.centerUnitOuterNormal()[0]<-0.9)
        {
          for (int i=0; i<NUM_GID_ENTRIES; ++i)
          {
            fixedElmtVector.push_back((ZOLTAN_ID_TYPE)elementGID[i]+1);
          }
          fixedProcVector.push_back(0);
        }
      }

    }
    // add one because not only neighbors are used in graph, but also entity itself
	  tempNborIndex[element_count+1] = tempNborIndex[element_count] + num_of_neighbors;

	  element_count++;
  }

  assert( tempNumMyVertices >= (int)element_count );

  // now copy into hypergraph structure
  hg_.numMyVertices = element_count;    // How many global elements there are
  hg_.numMyHEdges = element_count;	   // We have the same amount of Hyperedges
  std::swap(tempVtxGID, hg_.vtxGID);
  std::swap(tempVtxWeight, hg_.vtxWEIGHT);
  std::swap(tempEdgeGID, hg_.edgeGID);
  std::swap(tempNborIndex, hg_.nborIndex);

  hg_.numAllNbors = tempNborGID.size()/NUM_GID_ENTRIES;
  hg_.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * tempNborGID.size());
  std::copy(tempNborGID.begin(), tempNborGID.end(),hg_.nborGID);
  hg_.nborWEIGHT = (float *)malloc(sizeof(float) * tempNborGID.size());
  std::copy(tempNborWeight.begin(), tempNborWeight.end(),hg_.nborWEIGHT);
  hg_.nborPROC = (int *)malloc(sizeof(int) * tempNborGID.size());
  std::copy(tempNborProc.begin(), tempNborProc.end(),hg_.nborPROC);

  ///////// WRITE THE FIXED ELEMENTS INTO THE PROVIDED STRUCTURE
  if (fix_bnd_) // fixing element requires using hypergraph partitioning (which is perhaps better anyway?)
  {
    hg_.fixed_elmts.fixed_GID.resize(fixedElmtVector.size());
    std::copy(fixedElmtVector.begin(), fixedElmtVector.end(), hg_.fixed_elmts.fixed_GID.begin());
    hg_.fixed_elmts.fixed_Process.resize(fixedProcVector.size());
    std::copy(fixedProcVector.begin(), fixedProcVector.end(), hg_.fixed_elmts.fixed_Process.begin());
    hg_.fixed_elmts.fixed_entities = hg_.fixed_elmts.fixed_Process.size();
  }
}


// implement the required zoltan callback functions
template< class Grid >
int ZoltanLoadBalanceHandle<Grid>::
get_num_fixed_obj(void *data, int *ierr)
{
  HGraphData *graph = (HGraphData *)data;
  return graph->fixed_elmts.fixed_entities;
}
template< class Grid >
void ZoltanLoadBalanceHandle<Grid>::
get_fixed_obj_list(void *data, int num_fixed_obj,
                   int num_gid_entries, ZOLTAN_ID_PTR fixed_gids, int *fixed_part, int *ierr)
{
  HGraphData *graph = (HGraphData *)data;
  *ierr = ZOLTAN_OK;

  for (int i=0; i<num_fixed_obj*num_gid_entries; i++)
  {
    fixed_gids[i] = graph->fixed_elmts.fixed_GID[i];
  }

  for (int i=0; i<num_fixed_obj; i++)
  {
    fixed_part[i] = graph->fixed_elmts.fixed_Process[i];
  }
}


template< class Grid >
int ZoltanLoadBalanceHandle<Grid>::
get_number_of_vertices(void *data, int *ierr)
{
  HGraphData *temphg = (HGraphData *)data;
  *ierr = ZOLTAN_OK;
  return temphg->numMyVertices;
}

template< class Grid >
void ZoltanLoadBalanceHandle<Grid>::
get_vertex_list(void *data, int sizeGID, int sizeLID,
                ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                int wgt_dim, float *obj_wgts, int *ierr)
{
  int i;

  HGraphData *temphg= (HGraphData *)data;
  *ierr = ZOLTAN_OK;

  for (i=0; i<temphg->numMyVertices*sizeGID; i++)
  {
    globalID[i] = temphg->vtxGID[i];
  }

  for (i=0; i<temphg->numMyVertices; i++)
  {
    localID[i] = i;
    if (wgt_dim == 1)
      obj_wgts[i] = temphg->vtxWEIGHT[i];
  }
}

template< class Grid >
void ZoltanLoadBalanceHandle<Grid>::
get_num_edges_list(void *data, int sizeGID, int sizeLID,
                   int num_obj,
                   ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                   int *numEdges, int *ierr)
{
  HGraphData *temphg = (HGraphData *)data;
  *ierr = ZOLTAN_OK;
  for (int i=0;i<num_obj;++i)
    numEdges[i] = temphg->nborIndex[i+1]-temphg->nborIndex[i];
}
template< class Grid >
void ZoltanLoadBalanceHandle<Grid>::
get_edge_list(void *data, int sizeGID, int sizeLID,
              int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
              int *num_edges,
              ZOLTAN_ID_PTR nborGID, int *nborProc,
              int wgt_dim, float *ewgts, int *ierr)
{
  HGraphData *temphg = (HGraphData *)data;
  *ierr = ZOLTAN_OK;
  int k=0;
  for (int i=0;i<num_obj;++i)
  {
    int l = temphg->nborIndex[i];
    for (int j=0;j<num_edges[i];++j)
    {
      nborGID[k]  = temphg->nborGID[l];
      nborProc[k] = temphg->nborPROC[l];
      if (wgt_dim==1)
        ewgts[k] = temphg->nborWEIGHT[l];
      ++l;
      ++k;
    }
  }
}

#endif // HAVE_ZOLTAN
#endif // #ifndef LOADBALNCE_ZOLTAN_HH
