// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>

// Local includes
#include "piecewisefunction.hh"
#include "adaptation.hh"

using namespace Dune;

int main( int argc, char** argv )
try
{
  auto& mpiHelper = MPIHelper::instance( argc, argv );
  const int rank = mpiHelper.rank();

  std::string filename;
  if( argc > 1 )
  {
    filename = std::string(argv[1]);
  }
  else
  {
    std::cerr << "usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }

  using ALU3dSimplex = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
  Dune::GridPtr< ALU3dSimplex > gridPtr( filename );

  // gridPtr.loadBalance();

  std::cout << "P[ " << rank << " ] parameters = " << gridPtr.nofParameters( 0 ) << std::endl;
  auto lvlView = gridPtr->levelGridView( 0 );
  const auto end = lvlView.end< 0 > ();
  for( auto it = lvlView.begin< 0 > (); it != end; ++it )
  {
    const auto& entity = *it;
    std::cout << "P[ " << rank << " ], entity " << lvlView.indexSet().index( entity );
    const auto& param  = gridPtr.parameters( entity );
    for( const auto& p : param )
    {
      std::cout << " " << p;
    }
    std::cout << std::endl;
  }


  ALU3dSimplex &grid = *gridPtr;

  /* ... some global refinement steps */
  // if( verboseRank )
  // std::cout << "globalRefine: " << startLevel << std::endl;
  // grid.globalRefine( startLevel );

  /* get view to leaf grid */
  typedef ALU3dSimplex::Partition< Dune::Interior_Partition >::LeafGridView GridView;
  GridView gridView = grid.leafGridView< Dune::Interior_Partition >();

  /* construct data vector for solution */
  typedef PiecewiseFunction< GridView, Dune::FieldVector< double, 2 > > VectorType;
  VectorType TestSolution(gridView);
  TestSolution.clear();


   float version;
   int rc = Zoltan_Initialize((int)0, (char**)0, &version);
   if (rc != ZOLTAN_OK){
     printf("sorry zoltan did not initialize successfully...\n");
     MPI_Finalize();
     exit(1);
   }


  typedef ZoltanLoadBalanceHandle<ALU3dSimplex> LoadBalancer;
  LoadBalancer ldb(grid);
  ldb.repartition();

  typedef LeafAdaptation< ALU3dSimplex, VectorType, LoadBalancer > AdaptationType;
  AdaptationType adaptation( grid, ldb );

  adaptation(TestSolution);
  

  return 0;
}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e) {
  std::cerr << e.what() << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
