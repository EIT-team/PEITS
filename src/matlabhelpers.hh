/**************************************************************************
**       Title: matlabhelper.hh
**    $RCSfile: matlabhelper.hh,v $
**   $Revision: 1.4 $$Name:  $
**        Date: 28.11.2005
**   Copyright: GPL Bernard Haasdonk
** Description: collection of auxiliary functions for writing 
**              vector/matrix data for further external processing, e.g. in
**              MATLAB
**************************************************************************/

#ifndef __MATLABHELPERS_HH__
#define __MATLABHELPERS_HH__

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <config.h>


namespace Dune
{

  namespace Fem
  {

    /*======================================================================*/
    /*! @ingroup DiscFuncIO
     *  \class MatlabHelper
     *  \brief The MatlabHelper class provides functionality for exporting
     *         Dune Structures to Matlab.
     *
     *   These here are all written by MJE for the application in the 
     *   forward solver for Electrical Impedance Tomography dubbed
     *   "ultrasolver".
     */
    /*======================================================================*/

    class MatlabHelper
    {
      public:
  
	      template <class SigmaFunctionType, class FVSpace>
        int saveSigmaVectorBinary(const char* filename, SigmaFunctionType& sigma, SigmaFunctionType& indices, FVSpace& fvspace)
        /* Save the sigma vector and the corresponding element ID to a binary file, which can be loaded in Matlab */
        { 
		      typedef typename FVSpace::IteratorType IteratorType;
        	typedef typename IteratorType::Entity EntityType;

		      int myRank, numProcs;
      		MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

		      std::vector<double> sigmavector, IDvector;

        	IteratorType end = fvspace.end();
        	for( IteratorType it = fvspace.begin(); it != end; ++it )
        	{
	      		const EntityType &entity = *it;
          	typename SigmaFunctionType::LocalFunctionType sigmalocal = sigma.localFunction( entity );
          	typename SigmaFunctionType::LocalFunctionType indiceslocal = indices.localFunction( entity );
			      sigmavector.push_back(sigmalocal[0]);
			      IDvector.push_back(indiceslocal[0]);
        	}

		      // MPI COMMUNICATE THE VECTOR
		      std::vector<int> sizevector(numProcs);
		      int sigmasize = sigmavector.size();
		      MPI_Allgather(&sigmasize, 1, MPI_INT, &sizevector.front(), 1, MPI_INT, MPI_COMM_WORLD);
		      std::vector<int> displacements(numProcs);
        	displacements[0] = 0;
		      for (int process = 1; process < numProcs; ++process)
			      displacements[process] = std::accumulate(&sizevector[0],&sizevector[process],0);
		      int totalsize = std::accumulate(sizevector.begin(),sizevector.end(),0);
		      std::vector<double> totalsigmavector(totalsize), totalIDvector(totalsize);

		      MPI_Allgatherv(&sigmavector[0], sigmasize, MPI_DOUBLE, &totalsigmavector[0], &sizevector[0], &displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
		      MPI_Allgatherv(&IDvector[0], sigmasize, MPI_DOUBLE, &totalIDvector[0], &sizevector[0], &displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);

		      // WRITE THE VECTOR
          if (myRank == 0)
		      {    
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    

            // write magic number: type of binary file: 
            // DDV (Dune Dof Vector)            
            fid.write("DDV",3);
                
            // for debugging purposes: write an int and a double
            int magicint =    111;
            double magicdouble = 111.0;
            fid.write((char*)&magicint,sizeof(int));
            fid.write((char*)&magicdouble,sizeof(double));
                
            // write number of entries
            fid.write((char*)&totalsize,sizeof(int));
                
            // iterate over all elements defining the function
            for (int it = 0; it<totalsize; ++it)
            {
              double entry = totalIDvector[it];
              fid.write((char*)&entry,sizeof(double));
              entry = totalsigmavector[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 

            // write end-of-file marker            
            fid.write("EOF",3);

            int status = fid.good();
                
            fid.close();  

            return status;
          }
		      else
			      return 0;
        }

        template <class SigmaFunctionType, class FVSpace>
        int loadSigmaVectorBinary(const char* filename, SigmaFunctionType& sigma, SigmaFunctionType& indices, FVSpace& fvspace)
        /* Save the sigma vector and the corresponding element ID to a binary file, which can be loaded in Matlab */
        { 
		      typedef typename FVSpace::IteratorType IteratorType;
        	typedef typename IteratorType::Entity EntityType;

		      int myRank, numProcs;
      		MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

		      std::vector<double> sigmavector, IDvector;

          // WRITE OLD SIGMAVALUES AND INDICES INTO VECTORS
        	IteratorType end = fvspace.end();
        	for( IteratorType it = fvspace.begin(); it != end; ++it )
        	{
	      		const EntityType &entity = *it;
          	typename SigmaFunctionType::LocalFunctionType sigmalocal = sigma.localFunction( entity );
          	typename SigmaFunctionType::LocalFunctionType indiceslocal = indices.localFunction( entity );
			      sigmavector.push_back(sigmalocal[0]);
			      IDvector.push_back(indiceslocal[0]);
        	}

		      // MPI COMMUNICATE THE VECTOR SIZES AND RESULTING DISPLACEMENTS
		      std::vector<int> sizevector(numProcs);
		      int sigmasize = sigmavector.size();
		      MPI_Allgather(&sigmasize, 1, MPI_INT, &sizevector.front(), 1, MPI_INT, MPI_COMM_WORLD);
		      std::vector<int> displacements(numProcs);
        	displacements[0] = 0;
		      for (int process = 1; process < numProcs; ++process)
			      displacements[process] = std::accumulate(&sizevector[0],&sizevector[process],0);
		      int totalsize = std::accumulate(sizevector.begin(),sizevector.end(),0);
		      std::vector<double> totalsigmavector(totalsize), totalIDvector(totalsize);

		      // READ THE VECTOR
          if (myRank == 0)
		      {    
            // open file for writing  
            std::ifstream fid(filename, std::ios::binary | std::ios::in);
            int status = fid.good();
            if (!status)
            {
              std::cout << "\n The separate sigma file was not found! \n" << std::endl;
              return 1;
            }

            // DDV (Dune Dof Vector)
            static char buffer[3];
            fid.read(buffer,3);
            if (strcmp(buffer,"DDV") != 0)
              std::cout << "Problem with the file: DDV not written" << std::endl;
                
            // for debugging purposes: read an int and a double
            int magicint;
            double magicdouble;
            fid.read((char*)&magicint,sizeof(int));
            fid.read((char*)&magicdouble,sizeof(double));
            if ( magicint != 111 || magicdouble != 111.0 )
              std::cout << "Problem with the file: Magic Numbers wrong" << std::endl;
                
            // read number of entries
            int totalsizefile;
            fid.read((char*)&totalsizefile,sizeof(int));
            if (totalsizefile != totalsize)
              std::cout << "Problem with the file: The number of sigma values does not correspond" << std::endl;
                
            // iterate over all elements defining the function
            for (int it = 0; it<totalsizefile; ++it)
            {
              fid.read((char*)&totalIDvector[it],sizeof(double));
              fid.read((char*)&totalsigmavector[it],sizeof(double));
            } // end element iteration 

            // read end-of-file marker
            fid.read(buffer,3);
            if (strcmp(buffer,"EOF") != 0)
              std::cout << "Problem with the file: EOF not written" << std::endl;
                
            fid.close();
          }

          // SCATTER THE CORRESPONDING PARTS TO ALL PROCESSES
          std::vector<double> sigmavectorNEW(sizevector[myRank]), IDvectorNEW(sizevector[myRank]);
          MPI_Scatterv(&totalsigmavector[0], &sizevector[0], &displacements[0], MPI_DOUBLE, &sigmavectorNEW[0], sizevector[myRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		      MPI_Scatterv(&totalIDvector[0], &sizevector[0], &displacements[0], MPI_DOUBLE, &IDvectorNEW[0], sizevector[myRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

          // TEST IF ELEMENT ID'S CORRESPOND AND THEN UPDATE THE SIGMA VALUES// WRITE OLD SIGMAVALUES AND INDICES INTO VECTORS
          int index = 0;
        	for( IteratorType it = fvspace.begin(); it != end; ++it )
        	{
	      		const EntityType &entity = *it;
          	typename SigmaFunctionType::LocalFunctionType sigmalocal = sigma.localFunction( entity );
          	typename SigmaFunctionType::LocalFunctionType indiceslocal = indices.localFunction( entity );
            if (indiceslocal[0] == IDvectorNEW[index])
              sigmalocal[0] = sigmavectorNEW[index];
            else
              std::cout << "Indices do not add up!" << std::endl;
            index++;
        	}
          return 0;
        }


        void saveElectrodeVoltagesBinary(const char* filename, std::vector<double> evolts, int mode)
	      /*
	      Writes the first line of the Jacobian which contains "mode" lines.
	      If "mode" is zero, it just appends another line of the Jacobian.
	      If "mode" is negative, the Jacobian writing will be completed and the EOF set.
	      */
        { 
		      int myRank;
        	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
		      if (myRank>0)
			      return;

		      if (mode > 0)
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    

            // write magic number: type of binary file: 
            // DDV (Dune Electrode Voltages)            
            fid.write("DEV",3);
                  
            // for debugging purposes: write an int and a double
            int magicint =    111;
            double magicdouble = 111.0;
            fid.write((char*)&magicint,sizeof(int));
            fid.write((char*)&magicdouble,sizeof(double));
                  
            // write number of entries
            int ndofs = evolts.size();
            fid.write((char*)&ndofs,sizeof(int));

            // write number of rows
            int nrows = mode;
            fid.write((char*)&nrows,sizeof(int));
                  
            // iterate over all elements defining the function
			      double entry;
            for (int it = 0; it<ndofs; ++it)
            {
              entry = evolts[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 

			      if (mode == 1)
				      fid.write("EOF",3);
                  
            fid.close();
          }
		      else if (mode == 0)
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out | std::ios::app);    

            // iterate over all elements defining the function
            int ndofs = evolts.size();
			      double entry;
            for (int it = 0; it<ndofs; ++it)
            {
              entry = evolts[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 
                  
            fid.close();
		      }
		      else
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out | std::ios::app);    

            // iterate over all elements defining the function
            int ndofs = evolts.size();
			      double entry;
            for (int it = 0; it<ndofs; ++it)
            {
              entry = evolts[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 

            // write end-of-file marker            
            fid.write("EOF",3);
                  
            fid.close();
		      }
	      }



        template <class SigmaFunctionType, class FVSpace>
        void saveJacobianBinary(const char* filename, SigmaFunctionType& func, std::vector<double> elecJac, int mode, FVSpace& fvspace)
	      /*
	      Writes the first line of the Jacobian which contains "mode" lines.
	      If "mode" is zero, it just appends another line of the Jacobian.
	      If "mode" is negative, the Jacobian writing will be completed and the EOF set.
	      */
        { 
		      typedef typename FVSpace::IteratorType IteratorType;
        	typedef typename IteratorType::Entity EntityType;

		      int myRank, numProcs;
      		MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

		      std::vector<double> jacrow;

        	IteratorType end = fvspace.end();
        	for( IteratorType it = fvspace.begin(); it != end; ++it )
        	{
	      		const EntityType &entity = *it;
          	typename SigmaFunctionType::LocalFunctionType jaclocalrow = func.localFunction( entity );
			      jacrow.push_back(jaclocalrow[0]);
        	}

		      // MPI COMMUNICATE THE VECTOR
		      std::vector<int> sizevector(numProcs);
		      int jacrowsize = jacrow.size();
		      MPI_Allgather(&jacrowsize, 1, MPI_INT, &sizevector.front(), 1, MPI_INT, MPI_COMM_WORLD);
		      std::vector<int> displacements(numProcs);
        	displacements[0] = 0;
		      for (int process = 1; process < numProcs; ++process)
			      displacements[process] = std::accumulate(&sizevector[0],&sizevector[process],0);
		      int totalsize = std::accumulate(sizevector.begin(),sizevector.end(),0);
		      std::vector<double> totaljacrow(totalsize);

		      MPI_Allgatherv(&jacrow[0], jacrowsize, MPI_DOUBLE, &totaljacrow[0], &sizevector[0], &displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);

          if (myRank>0)
			      return;

		      if (mode > 0)
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    

            // write magic number: type of binary file: 
            // DDV (Dune Jacobian Matrix)            
            fid.write("DJM",3);
            
            // for debugging purposes: write an int and a double
            int magicint =    111;
            double magicdouble = 111.0;
            fid.write((char*)&magicint,sizeof(int));
            fid.write((char*)&magicdouble,sizeof(double));
            
            // write number of entries
            int ndofs = totaljacrow.size() + elecJac.size();
            fid.write((char*)&ndofs,sizeof(int));

            // write number of rows
            int nrows = mode;
            fid.write((char*)&nrows,sizeof(int));
            
            // iterate over all elements defining the function
			      double entry;
            for (unsigned int it = 0; it<totaljacrow.size(); ++it)
            {
              entry = totaljacrow[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 

            // iterate over all electrodes and write x and y Jacobian entry
            for (unsigned int it = 0; it<elecJac.size(); ++it)
            {
              entry = elecJac[it];
              fid.write((char*)&entry,sizeof(double));
            }

			      if (mode == 1)
				      fid.write("EOF",3);
            
            fid.close();
          }
		      else if (mode == 0)
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out | std::ios::app);    

            // iterate over all elements defining the function
            int ndofs = totaljacrow.size();
			      double entry;
            for (int it = 0; it<ndofs; ++it)
            {
              entry = totaljacrow[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration

            // iterate over all electrodes and write x and y Jacobian entry
            for (unsigned int it = 0; it<elecJac.size(); ++it)
            {
              entry = elecJac[it];
              fid.write((char*)&entry,sizeof(double));
            }
                  
            fid.close();
		      }
		      else
		      {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out | std::ios::app);    

            // iterate over all elements defining the function
            int ndofs = totaljacrow.size();
			      double entry;
            for (int it = 0; it<ndofs; ++it)
            {
              entry = totaljacrow[it];
              fid.write((char*)&entry,sizeof(double));
            } // end element iteration 

            // iterate over all electrodes and write x and y Jacobian entry
            for (unsigned int it = 0; it<elecJac.size(); ++it)
            {
              entry = elecJac[it];
              fid.write((char*)&entry,sizeof(double));
            }

            // write end-of-file marker            
            fid.write("EOF",3);
            
            fid.close();
		      }
	      }

      private:

    }; // end class MatlabHelper

  } // namespace Fem
  
} // namespace Dune

#endif
