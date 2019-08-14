// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

// Needed for timing
#include <sys/time.h>

// Grid part declarations
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Needed for the femscheme outputs
#include <dune/fem/io/file/dataoutput.hh>

// Needed for the electrodes helper
#include <dune/fem/io/parameter.hh>

// Dynamic vector declaration for the sigmavector
#include <dune/common/dynvector.hh>
// For the function space
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/space/finitevolume.hh>

// Adaptation scheme
// #include "adaptation.hh"
#include "adaptation.hh"
// Type declarations for the electrodes
#include "electrode_helpers.hh"
// For the mathematical model used and elliptic solver
#include "femscheme.hh"
// Matlab helper
// TODO Uncomment when using Matlab
// #include "matlabhelpers.hh"
#include "poisson.hh"
// Piecewise class declaration
#include "piecewisefunction.hh"
// Declarations for Zoltan partitioning
#include "zoltaninterface.hh"

using namespace Dune;

template <class HGridType>
double algorithm (HGridType &grid,
                  Dune::DynamicVector<double> &sigmavector,
                  Dune::DynamicVector<double> &elementIDvector,
                  int step,
                  ElectrodePositions &electrodes,
                  CurrentProtocol &current_protocol)
{
  // TODO: Uncomment when using Matlab
  // Dune::Fem::MatlabHelper mhelp;

  
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 1 > FunctionSpaceType;

  // write sigma vector to a discrete function
  typedef Dune::Fem::FiniteVolumeSpace<FunctionSpaceType,GridPartType,0,Dune::Fem::CachingStorage> FVSpaceType;
  typedef Dune::DynamicVector<double> SigmaVectorType;
  typedef Dune::Fem::VectorDiscreteFunction<FVSpaceType,SigmaVectorType> SigmaFunctionType;
  FVSpaceType fvspace(gridPart);

  // get the time to put into name of save files later on
  char time_buf[21];
  time_t now;
  time(&now);
  strftime(time_buf, 21, "%Y-%m-%d_%H:%M:%S", gmtime(&now));

  const SigmaFunctionType sigma("sigma", fvspace, sigmavector);
  const SigmaFunctionType elementID("elementID", fvspace, elementIDvector);

  // write sigma vector to matlab readable file
  std::string sigmaname("sigmavector.bin");
  sigmaname.insert(11,time_buf);
  const string outputpath = Dune::Fem::Parameter::getValue< string >( "fem.io.outputpath" );
  sigmaname.insert(0,outputpath);
  // TODO Uncomment when using Matlab
  // mhelp.saveSigmaVectorBinary(sigmaname.c_str(), sigma, elementID, fvspace);  

  // if the sigma values should be loaded from a separate file do it here
  if (Dune::Fem::Parameter::getValue< bool >("fem.io.load_sigma_separately"))
  {
    const string sigmaUpdate = Dune::Fem::Parameter::getValue< string >( "fem.io.separate_sigma_file" );
    // TODO: Uncomment when using Matlab
    // int problemwithsigma = mhelp.loadSigmaVectorBinary(sigmaUpdate.c_str(), sigma, elementID, fvspace);
    // if (problemwithsigma)
      // return 1;
  }

  // type of the mathematical model used
  typedef DiffusionModel< FunctionSpaceType, GridPartType, SigmaFunctionType > ModelType;

  // define problem
  typedef typename ModelType :: ProblemType ProblemType ;
  ProblemType* problemPtr = 0 ;
  problemPtr = new EIT_solver< FunctionSpaceType > ();
  assert( problemPtr );
  ProblemType& problem = *problemPtr ;
  
  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart, sigma, elementID );
  if (implicitModel.error_)
    return 1;

  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
  typedef Dune::Fem::SparseRowMatrixObject<DiscreteFunctionSpaceType, DiscreteFunctionSpaceType> MatrixType;
  typedef Dune::Fem::SparseRowMatrix<double> FinalMatrixType;

#if HAVE_DUNE_ISTL && WANT_ISTL
#warning using ISTL
  typedef FemScheme< ModelType,istl > SchemeTypeOne;
#elif HAVE_PETSC && WANT_PETSC
#warning using PETSC
  typedef FemScheme< ModelType,petsc > SchemeTypeOne;
#else
  typedef FemScheme< ModelType,femoem > SchemeTypeOne;
#endif

  SchemeTypeOne scheme( gridPart, implicitModel, electrodes, current_protocol );
  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );

  //! input/output tuple and setup datawriter
  bool write_vtk = Dune::Fem::Parameter::getValue< bool >( "write.vtk" );
  typedef std::tuple< const typename SchemeTypeOne::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  timeval atime,astop;
  if( Dune::Fem::MPIManager::rank() == 0)
    {
      cout << "Assembling the System Matrix... ";
      gettimeofday(&atime, NULL);
    }
  scheme.assemble();
  if( Dune::Fem::MPIManager::rank() == 0)
    {
      gettimeofday(&astop, NULL);
      double assemblytimediff = astop.tv_sec-atime.tv_sec + (astop.tv_usec-atime.tv_usec)/1000000.0;
      cout << "done!        It took: " << assemblytimediff << " seconds" << endl << endl;
    }

  int noUnique = current_protocol.unique_patterns_.size();

  // this is the container for the unique current injections
  boost::ptr_vector< typename SchemeTypeOne::DiscreteFunctionType > injectionresults;
  boost::ptr_vector< vector<double> > elecVolts;

  if (Dune::Fem::MPIManager::rank() == 0)
    cout << "Now Solving the Forward for " << noUnique << " Unique Injections:" << endl;

  for (int protocol_step = 0; protocol_step < noUnique; ++protocol_step)
    {
      // if only electrode voltages are required, do NOT compute adjacent fields
      if (!(Dune::Fem::Parameter::getValue< bool >("fem.io.do_jacobian") || Dune::Fem::Parameter::getValue< bool >("fem.io.do_electrode_jacobian")))
        {
          bool indicator = false;
          for (int iterator = 0; iterator < current_protocol.current_protocol_.size(); iterator++)
            {
              if (current_protocol.retrace_current_pattern_[iterator][0] == protocol_step)
                { indicator = true; break; }
            }
          if (!indicator)
            {
              if( Dune::Fem::MPIManager::rank() == 0)
                cout << "skipping a forward simulation " << protocol_step << endl;
              
              // save previous injection result (is not going to be used - just there to fill space)
              injectionresults.push_back(new typename SchemeTypeOne::DiscreteFunctionType (scheme.solution()));

              // add empty vector to elecVolts, again just to fill the space
              elecVolts.push_back(new vector<double> (electrodes.electrode_positions_.size()));

              continue;
            }
        }

      // setup the right hand side
      timeval stime, updatetime, resolvetime;
      if( Dune::Fem::MPIManager::rank() == 0)
        gettimeofday(&stime, NULL);
      scheme.prepare(protocol_step);
      if( Dune::Fem::MPIManager::rank() == 0)
        gettimeofday(&updatetime, NULL);

      // solve
      scheme.solve(protocol_step);
      if( Dune::Fem::MPIManager::rank() == 0)
        {
          gettimeofday(&resolvetime, NULL);
          double updatetimediff = updatetime.tv_sec-stime.tv_sec + (updatetime.tv_usec-stime.tv_usec)/1000000.0;
          double solvetimediff = resolvetime.tv_sec-updatetime.tv_sec + (resolvetime.tv_usec-updatetime.tv_usec)/1000000.0;
          cout << "    " << protocol_step << endl << "    updatetime = " << updatetimediff << " ; solvetime = " << solvetimediff << endl;
        }

      // write initial solve
      if (write_vtk && !protocol_step)
        dataOutput.write();

      // save current injection results (for a unit current => scale with injection current if needed)
      injectionresults.push_back(new typename SchemeTypeOne::DiscreteFunctionType (scheme.solution()));

      // compute the resulting electrode potentials
      elecVolts.push_back(new vector<double> (electrodes.electrode_positions_.size()));
      scheme.calculateElecPot(injectionresults.back(), protocol_step, elecVolts.back());
    
    }

  // define entities
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::EntitySeed     EntitySeed;
  const vector<vector<EntitySeed>> electrodeElements = scheme.electrodeElements();
  const vector<double> electrodeAreas = scheme.electrodeAreas();

  // save local stiffness matrices (and electrode jacobian templates, if necessary)
  JacobianRowCalculator<typename SchemeTypeOne::DiscreteFunctionType, SigmaFunctionType> jacRowCalculator( injectionresults[0], electrodes, electrodeElements, electrodeAreas );
  SigmaVectorType jacobian_row_vector(sigmavector.size());
  SigmaFunctionType jacobianRow("jacrow",fvspace,jacobian_row_vector);

  int noProtocols = current_protocol.current_protocol_.size();
  if (Dune::Fem::MPIManager::rank() == 0)
    cout << endl << "Now Calculating Electrode Voltages and Jacobian Matrix for " << noProtocols << " Protocol Steps:" << endl;

  for (int protocol_step = 0; protocol_step < noProtocols; ++protocol_step)
    {
      if (Dune::Fem::MPIManager::rank() == 0)
        cout << "    " << protocol_step << endl;
      // retrace the correct current patterns and calculate the factor we have to multiply them with
      int indexdrive = current_protocol.retrace_current_pattern_[protocol_step][0];
      int indexmeasure = current_protocol.retrace_current_pattern_[protocol_step][1];
      double factor = current_protocol.input_current_, fwdfactor; // all injectionresults are for unit current, i.e. we need to scale the drive current now
      if (indexdrive<0)
        {indexdrive = -(indexdrive+1); factor = -factor;}
      fwdfactor = factor;
      if (indexmeasure<0)
        {indexmeasure = -(indexmeasure+1); factor = -factor;}

      // write electrode voltages
      if (Dune::Fem::Parameter::getValue< bool >("fem.io.do_elec_volts"))
        scheme.writeElecPot ( elecVolts[indexdrive], indexmeasure, fwdfactor, factor/fwdfactor, protocol_step );

      if (Dune::Fem::Parameter::getValue< bool >("fem.io.do_jacobian"))
        {
          timeval jacstart, jacstop;
          if( Dune::Fem::MPIManager::rank() == 0)
            {
              cout << "    Calculating Jacobian Row... ";
              gettimeofday(&jacstart, NULL);
            }

          // compute "normal" jacobian row
          jacRowCalculator.getJacobianRow( injectionresults[indexdrive], injectionresults[indexmeasure], jacobianRow, factor);
          
          // if wished for, compute the jacobian matrix with respect to the electrode boundaries for both mesh surface coordinates: x and y
          vector<double> elecJacRow;
          if (jacRowCalculator.do_electrode_jacobian_)
            jacRowCalculator.getElecJacRow( injectionresults[indexdrive], elecVolts[indexdrive], injectionresults[indexmeasure], elecVolts[indexmeasure], elecJacRow, factor);

          if( Dune::Fem::MPIManager::rank() == 0)
            {
              gettimeofday(&jacstop, NULL);
              double jactimediff = jacstop.tv_sec-jacstart.tv_sec + (jacstop.tv_usec-jacstart.tv_usec)/1000000.0;
              cout << "done!          It took: " << jactimediff << " seconds" << endl;
            }


          // SAVE EVERYTHING TO MATLAB-READABLE FILES
          std::string rhsfn("jacobian.bin");
          rhsfn.insert(8,time_buf);
          const string outputpath = Dune::Fem::Parameter::getValue< string >( "fem.io.outputpath" );
          rhsfn.insert(0,outputpath);
          int mode = 0;
          if (!protocol_step)
            mode = noProtocols;
          else if (protocol_step == noProtocols-1)
            mode = -1;

          timeval jacwritestart, jacwritestop;
          if( Dune::Fem::MPIManager::rank() == 0)
            {
              cout << "    Writing Jacobian Row... ";
              gettimeofday(&jacwritestart, NULL);
            }

          // TODO: Uncomment when Matlab use
          // mhelp.saveJacobianBinary(rhsfn.c_str(), jacobianRow, elecJacRow, mode, fvspace);

          if( Dune::Fem::MPIManager::rank() == 0)
            {
              gettimeofday(&jacwritestop,NULL);
              double jacwritetimediff = jacwritestop.tv_sec-jacwritestart.tv_sec + (jacwritestop.tv_usec-jacwritestart.tv_usec)/1000000.0;
              cout << "done!              It took: " << jacwritetimediff << " seconds" << endl;
            }
        }
    } // end of protocal iteration

  // calculate error
  double error = 0 ;
  
  //cout << "Iteration completed" << endl;
  return error ;
}

// main
// ----

int main(int argc, char** argv)
{
  try{
    // // Maybe initialize MPI
    // Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    // std::cout << "Hello World! This is dune-peits." << std::endl;
    // if(Dune::MPIHelper::isFake)
    //   std::cout<< "This is a sequential program." << std::endl;
    // else
    //   std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
    //     <<" processes!"<<std::endl;

    // Initialize MPI, if necessary
    Dune::Fem::MPIManager::initialize( argc, argv );

    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // Start timing
    timeval starttime, stoptime;
    if( myRank == 0 ) gettimeofday(&starttime, NULL);

    // Append overloaded parameters from the command line
    Dune::Fem::Parameter::append( argc, argv );

    // append default parameter file  
    Dune::Fem::Parameter::append( "../data/parameter" );
    
    // type of hierarchical grid 
    //typedef Dune :: AlbertaGrid< 2 , 2 > GridType;
    typedef Dune::GridSelector::GridType  HGridType ;

    /////////// first, read in the electrode helpers

    // read in electrode positions
    string electrodes_string;
    if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
      electrodes_string = Dune::Fem::Parameter::getValue< string >( "electrode.nodes" );
    else
      electrodes_string = Dune::Fem::Parameter::getValue< string >( "electrode.positions" );
    ElectrodePositions electrodes(electrodes_string.c_str());
    if (electrodes.error_)
      return 1;

    // read in current protocol
    const string current_protocol_string = Dune::Fem::Parameter::getValue< string >( "current.protocol" );
    CurrentProtocol current_protocol(current_protocol_string.c_str());
    if (current_protocol.error_)
      return 1;

    // load the mesh
    HGridType *hgrid;
    typedef Dune::DynamicVector<double> SigmaVectorType;
    SigmaVectorType eldat(0),elementID(0);

    if (Dune::Fem::Parameter::getValue< bool >("fem.io.loadPartitions"))
      {
        string fileName = Dune::Fem::Parameter::getValue< string >("fem.io.macroGrid");
        string partitionPath = Dune::Fem::Parameter::getValue< string >("fem.io.partitionPath");
        std::stringstream partitionName;
        partitionName << partitionPath << fileName << "." << numProcs << "_" << myRank;

        timeval loadpartstart, loadpartstop;
        if( myRank == 0 )
          {
            cout << endl << "Loading grid: " << fileName << " in from partitions... ";
            gettimeofday(&loadpartstart, NULL);
          }

        Dune::Fem::BinaryFileInStream instream(partitionName.str()+".grid");
        hgrid = Dune::BackupRestoreFacility<HGridType>::restore(instream.stream());
        eldat.resize(hgrid->size(0));
        elementID.resize(hgrid->size(0));

        ifstream indata(partitionName.str().c_str());
        indata.read( (char*)(&eldat[0]), sizeof(double)*eldat.size());
        indata.read( (char*)(&elementID[0]), sizeof(double)*eldat.size());

        if( myRank == 0 )
          {
            gettimeofday(&loadpartstop, NULL);
            double loadparttimediff = loadpartstop.tv_sec-loadpartstart.tv_sec + (loadpartstop.tv_usec-loadpartstart.tv_usec)/1000000.0;
            cout << "done!        It took: " << loadparttimediff << " seconds" << endl << endl;
          }
      }
    else
      {
        // create grid from DGF file
        string fileName = Dune::Fem::Parameter::getValue< string >("fem.io.macroGrid");
        string meshPath = Dune::Fem::Parameter::getValue< string >("fem.io.meshPath");
        std::stringstream meshName;
        meshName << meshPath << fileName;

        timeval dopartstart, dopartstop;
        if( myRank == 0 )
          {
            cout << endl << "Loading grid: " << meshName.str() << " - remember to load the created partitions for successive executions" << endl;
            gettimeofday(&dopartstart, NULL);
          }

        timeval firststart,firststop;
        if( myRank == 0 )
          {
            cout << "Now we have to create partitions: " << endl << "    Loading... ";
            gettimeofday(&firststart, NULL);
          }
        
        // construct macro using the DGF Parser 
        Dune::GridPtr< HGridType > gridPtr( meshName.str() );

        // do initial load balance
        gridPtr.loadBalance();

        ////////////////////////////////////////////////
        if( myRank == 0 )
          {
            gettimeofday(&firststop, NULL);
            double loadtimediff = firststop.tv_sec-firststart.tv_sec + (firststop.tv_usec-firststart.tv_usec)/1000000.0;
            cout << "            done!                  It took " << loadtimediff << " seconds" <<  endl;
          }
        ////////////////////////////////////////////////

        typedef HGridType::LeafGridView GridView;
        typedef GridView::IndexSet IndexSetType;

        // create Grid from DGF parser
        size_t nofElParams( 0 );
        GridView gridView = gridPtr->leafView();

        // write the conductivity and elementId into a vector that is supported by the partitioning
        typedef HGridType::Codim< 0 >::Entity Entity;
        typedef PiecewiseFunction< GridView, Dune::FieldVector< double, 2 > > TestVectorType;
        TestVectorType TestSolution( gridView );
        TestSolution.clear();

        const IndexSetType &indexSet = gridView.indexSet();
        nofElParams = gridPtr.nofParameters( 0 );

        if( nofElParams > 0 )
          {
            if( Dune::Fem::MPIManager::rank() == 0 )
              std::cout << "    Reading Element Parameters... ";
            const Dune::PartitionIteratorType partType = Dune::All_Partition;
            typedef GridView::Codim< 0 >::Partition< partType >::Iterator Iterator;
            const Iterator enditer = gridView.end< 0, partType >();

            for( Iterator iter = gridView.begin< 0, partType >(); iter != enditer; ++iter )
              {
                const std::vector< double > &param = gridPtr.parameters( *iter );
                assert( param.size() == nofElParams );
                
                const Entity &entity = *iter;
                TestVectorType::LocalDofVector dat(TestSolution[entity]);
                dat[0] = param[0];
                dat[1] = param[1];
                TestSolution.setLocalDofVector( entity, dat);
              }
            if( Dune::Fem::MPIManager::rank() == 0 )
              std::cout << "done!" << endl;
          }
        hgrid = gridPtr.release();
        HGridType &grid = *hgrid;

        ///////////////////////////////////////////////////
        timeval startzoltan,stopzoltan;
        if( myRank == 0 )
          {
            cout << "    Start Zoltan repartitioning... ";
            gettimeofday(&startzoltan, NULL);
          }
        ///////////////////////////////////////////////////

        int rc;
        int changes, numGidEntries, numLidEntries, numImport, numExport;
        ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
        int *importProcs, *importToPart, *exportProcs, *exportToPart;

        
        float ver;
        rc = Zoltan_Initialize((int)0,(char**)0, &ver);
        
        if (rc != ZOLTAN_OK){
          printf("sorry...\n");
          MPI_Finalize();
          exit(0);
        }

        // create Zoltan Hypergraph here
        HGRAPH_DATA hg;
        int NUM_GID_ENTRIES = 4;
        write_hypergraph(myRank, numProcs, grid, &hg, NUM_GID_ENTRIES, electrodes);
        //read_input_file(myRank, numProcs, "hypergraph.txt", &hg);
        struct Zoltan_Struct *zz;
        zz = Zoltan_Create(MPI_COMM_WORLD);

        // General parameters

        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
        Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   // partitioning method 
        Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); // version of method
        Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "4");// global IDs are integers
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");// local IDs are integers 
        Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); // export AND import lists 
        Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); // use Zoltan default vertex weights 
        Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");// use Zoltan default hyperedge weights 
        Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");

        // Application defined query functions
        Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &hg);
        Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &hg);
        Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, &hg);
        Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, &hg);

        // Register fixed object callback functions 
        if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_FIXED_OBJ_FN_TYPE,
                          (void (*)()) get_num_fixed_obj,
                          (void *) &hg) == ZOLTAN_FATAL)
          {
            return 1;
          }

        if (Zoltan_Set_Fn(zz, ZOLTAN_FIXED_OBJ_LIST_FN_TYPE,
                          (void (*)()) get_fixed_obj_list,
                          (void *) &hg) == ZOLTAN_FATAL)
          {
            return 1;
          }
        
        rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output)
                                 &changes,        // 1 if partitioning was changed, 0 otherwise 
                                 &numGidEntries,  // Number of integers used for a global ID 
                                 &numLidEntries,  // Number of integers used for a local ID 
                                 &numImport,      // Number of vertices to be sent to me 
                                 &importGlobalGids,  // Global IDs of vertices to be sent to me 
                                 &importLocalGids,   // Local IDs of vertices to be sent to me 
                                 &importProcs,    // Process rank for source of each incoming vertex 
                                 &importToPart,   // New partition for each incoming vertex 
                                 &numExport,      // Number of vertices I must send to other processes
                                 &exportGlobalGids,  // Global IDs of the vertices I must send 
                                 &exportLocalGids,   // Local IDs of the vertices I must send 
                                 &exportProcs,    // Process to which I send each of the vertices 
                                 &exportToPart);  // Partition to which each vertex will belong

        if( myRank == 0 )
          {
            gettimeofday(&stopzoltan, NULL);
            double zoltantimediff = stopzoltan.tv_sec-startzoltan.tv_sec + (stopzoltan.tv_usec-startzoltan.tv_usec)/1000000.0;
            cout << "done!          It took " << zoltantimediff << " seconds" << endl;
          }


        ZOLTAN_PARTITIONING new_partitioning;
        new_partitioning.changes = changes; // 1 if partitioning was changed, 0 otherwise 
        new_partitioning.numGidEntries = numGidEntries;  // Number of integers used for a global ID 
        new_partitioning.numExport = numExport;  // Number of vertices I must send to other processes
        new_partitioning.exportGlobalGids = exportGlobalGids;  // Global IDs of the vertices I must send 
        new_partitioning.exportProcs = exportProcs;    // Process to which I send each of the vertices

        typedef ZoltanLoadBalanceHandle<HGridType> LoadBalancer;
        LoadBalancer ldb(grid);


        Zoltan_Destroy(&zz);
        if (rc != ZOLTAN_OK){
          printf("sorry...\n");
          MPI_Finalize();
          Zoltan_Destroy(&zz);
          exit(0);
        }

        // create adaptation method
        typedef LeafAdaptation< HGridType, TestVectorType,  LoadBalancer> AdaptationType;
        AdaptationType adaptation( grid, ldb );

        // Apply the Zoltan partitioning
        adaptation( TestSolution );

        // Write back the new sigmavector
        const Dune::PartitionIteratorType partType = Dune::All_Partition;
        typedef GridView::Codim< 0 >::Partition< partType >::Iterator Iterator;
        eldat.resize( TestSolution.size() );
        elementID.resize( TestSolution.size() );
        const Iterator enditer = gridView.end< 0, partType >();

        for( Iterator iter = gridView.begin< 0, partType >(); iter != enditer; ++iter )
          {
            TestVectorType::LocalDofVector dat(TestSolution[*iter]);
            eldat[indexSet.index(*iter)] = dat[0];
            elementID[indexSet.index(*iter)] = dat[1];
          }
        timeval writestart;
        if( myRank == 0 )
          {
            cout << "    Writing the partitions to files... ";
            gettimeofday(&writestart, NULL);
          }

        string partitionPath = Dune::Fem::Parameter::getValue< string >("fem.io.partitionPath");
        std::stringstream partitionFile;
        partitionFile << partitionPath << fileName << "." << numProcs << "_" << myRank;
        Dune::Fem::BinaryFileOutStream outstream(partitionFile.str()+".grid");
        Dune::BackupRestoreFacility<HGridType>::backup(grid,outstream.stream());
        
        ofstream outdata(partitionFile.str().c_str());
        outdata.write( (char*)(&eldat[0]), sizeof(double)*eldat.size());
        outdata.write( (char*)(&elementID[0]), sizeof(double)*eldat.size());

        if( myRank == 0 )
          {
            gettimeofday(&dopartstop, NULL);
            double writetimediff = dopartstop.tv_sec-writestart.tv_sec + (dopartstop.tv_usec-writestart.tv_usec)/1000000.0;
            double totaltimediff = dopartstop.tv_sec-dopartstart.tv_sec + (dopartstop.tv_usec-dopartstart.tv_usec)/1000000.0;
            cout << "done!      It took " << writetimediff << " seconds" << endl;
            cout << "THE WHOLE PARTITIONING TOOK " << totaltimediff << " seconds" << endl << endl;
          }
      }
    HGridType &grid = *hgrid;

    timeval startsolving;
    if( myRank == 0 )
      gettimeofday(&startsolving, NULL);

    /////////  SECTION: RUN ALGORITHM AS MANY TIMES AS REQUESTED IN PARAMETER FILE

    // calculate first step
    algorithm ( grid, eldat, elementID, false, electrodes, current_protocol);

    if( myRank == 0 )
      {
        gettimeofday(&stoptime, NULL);
        double totalruntime = stoptime.tv_sec-starttime.tv_sec + (stoptime.tv_usec-starttime.tv_usec)/1000000.0;
        double totalsolvingtime = stoptime.tv_sec-startsolving.tv_sec + (stoptime.tv_usec-startsolving.tv_usec)/1000000.0;
        cout << endl << "TOTAL RUNTIME WAS " << totalruntime << " seconds" << endl;
        cout << "TOTAL TIME SPENT SETTING UP SYSTEM AND SOLVING WAS " << totalsolvingtime << " seconds" << endl;
      }

    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
