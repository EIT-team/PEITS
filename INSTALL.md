# Downloading and installing the required modules (Tested on Ubuntu 16.04)

1. Install required libraries.
   
   ```bash
   sudo apt-get update
   sudo apt-get install git gcc build-essential clang make autotools-dev \
                        autoconf automake libtool cmake scons pkg-config \
                        libblas-dev liblapack-dev gfortran libboost-dev \
                        mpi-default-bin mpi-default-dev libgmp-dev libopenmpi-dev \
                        gfortran autotools-dev automake libug-dev \
                        libmetis-dev libparmetis-dev
   ```

1. Create a new folder — e.g. PEITS_root (Parallel EIT Solver) where the modules
   will be installed.

   ```bash
   PEITS_root=/path/where/you/want/to/install/PEITS_root
   mkdir -p ${PEITS_root}
   cd ${PEITS_root}
   ```
   **Important** Change the `PEITS_root=/path/where/you/want/to/install/PEITS_root`
   variable to reflect the actual location of the PEITS installation. /e.g./,
   `PEITS_root=${HOME}/PEITS_root`.

1. Download the PETSc library:

   ```bash
   git clone -b maint https://bitbucket.org/petsc/petsc
   ```
   This will create a `petsc` subfolder in the PEITS_root directory.

1. The PETSc libary needs to be configured, from within the `pestc` folder.
   Installation has been tested using a previous version of PETSc (commit
   `8695de0` - the nearest release is v3.6.3), newer versions cause an error
   during configuration:

   ```bash
   cd petsc
   git checkout 8695de0
   ./configure --prefix=${PEITS_root}/petscBUILD \
               --with-x=0 --with-debugging=0 \
               -CFLAGS="-O3 -DNDEBUG -ffast-math" \
               --with-parmetis=1 --download-parmetis=yes \
               --with-hypre=1 --download-hypre=yes \
               --with-superlu_dist=1 --download-superlu_dist=yes \
               --with-mumps=1 --download-mumps=yes \
               --with-ml=1 --download-ml=yes \
               --with-metis=1 --download-metis=yes \
               --download-scalapack=yes --download-blacs=yes
   ```

   The `petscBUILD` folder will be created during the configuration process.


1. Build the PETSc source code.

   ```bash
   make all test
   make install
   ```

1. Download Zoltan library (v3.83 required) and extract to the PEITS_root folder.

   ```bash
   cd ${PEITS_root}
   wget http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz
   tar xf zoltan_distrib_v3.83.tar.gz --warning=no-unknown-keyword
   ```
   
   There may be some warning messages `'Ignoring unknown extended header keyword'`, these can be ignored.
   The Zoltan library can also be downloaded through [their website](http://www.cs.sandia.gov/~web1400/1400_download.html).

1. Configure and build the Zoltan library, must be run in a build subdirectory:

   ```bash
   mkdir Zoltan_v3.83/BUILD_DIR
   cd Zoltan_v3.83/BUILD_DIR
   ../configure --prefix=${PEITS_root}/Zoltan_v3.83/BUILD_DIR \
                --with-parmetis --with-parmetis-incdir=/usr/include/ \
                --with-parmetis-libdir=/usr/lib
   make everything
   make install
   ```
   
1. Download Parallel EIT Solver (PEITS) code into PEITS_root folder:

   ```bash
   cd ${PEITS_root}
   git clone https://github.com/EIT-team/PEITS
   ```

1. Edit the `PEITS/config.opts_example` file and change `PETSCPATH` and `ZOLTANPATH`
   to the appropriate directories by replacing `/home/username/PEITS_root/` with
   the directory where `PEITS_root` is located, save the edited file as
   `config.opts`.

   If you have installed metis and parmetis using the package manager as
   described above, then you can remove the lines that refer to them: /i.e./,
   `--with-metis=...` and `--with-parmetis=...`

1. Run install script. Located in PEITS_root/PEITS/

   ```bash
   ./INSTALL.sh
   ```

   The most likely cause of failure at this step is errors in the `config.opts`
   file, double check `PETSCPATH` and `ZOLTANPATH` if this step fails.

1. Test installation.

   A simple regression test is includes in the PEITS/tests folder. This will check
   that the default configuration is producing the expected output (by comparing
   file sizes).
   
   ```bash
   cd ${PEITS_root}/PEITS/tests
   sh initial_test.sh
   ```
   
   If the test is successful, `file sizes match - test OK` will be displayed,
   meaning the installation process _should_ have been carried out correctly.
   
   If there is an error regarding loading `libparmetis.so` then we need to add the
   petsc library directory to the environmental path:
   
   ```bash
   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PEITS_root}/petscBUILD/lib/
   export LD_LIBRARY_PATH
   ```
   
1. Running PEITS. If installation has been completed successfully, the solver
   can be called from the `${PEITS_root}/PEITS/src` folder:

   ```bash
   cd ${PEITS_root}/PEITS/src
   mpirun -np 2 ./dune_peits
   ```
   
   where `-np` specifies the number of parallel processors the solver should be run
   on.
   
   If the solver needs an unreasonably long time for the assembly of the system
   matrix, then the pre-allocation of memory in PETSc needs to be adjusted. In
   _file `${PEITS_root}/PEITS/dune-fem-1.4.0/dune/fem/misc/petsc/petsccommon.hh` the
   number of allocated non-zeros can be changed in the command
   `MatMPIAIJSetPreallocation(mat,100,PETSC_NULL,40,PETSC_NULL)`. A safe way of
   adjusting this is to use very high numbers (e.g. 1000 and 150) and then running
   the solver with the option `-info`, which outputs the precise number of non-zeros
   required on the used mesh. Also, on some meshes ML preconditioning fails on some
   numbers of parallel processes. If this happens, either the number of processes
   can be changed or hypre preconditioning can be used.
