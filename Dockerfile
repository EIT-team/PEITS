FROM ubuntu:18.10
LABEL maintainer="r.guichard@ucl.ac.uk"

# Install basic dependencies
RUN apt-get -qq update
RUN apt-get install -y \
    git \
    wget \
    gcc-8 \
    g++-8 \
    gfortran-8 \
    cmake \
    libblas-dev \
    liblapack-dev \
    libboost-dev \
    mpi-default-bin \
    mpi-default-dev \
    libgmp-dev \
    libopenmpi-dev \
    pkg-config \
    m4 \
    libmetis-dev \
    libparmetis-dev

# Install PETSc from packet manager
RUN apt-get install -y petsc-dev

# Basic Zoltan install with CMake
RUN git clone https://github.com/trilinos/Trilinos.git
RUN mkdir Trilinos/BUILD_DIR
WORKDIR "/Trilinos/BUILD_DIR"
RUN cmake -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF -DTrilinos_ENABLE_Zoltan:BOOL=ON ..
RUN make && make install

# Get and untar all the necessary dune tarballs
WORKDIR "/"
RUN wget https://www.dune-project.org/download/2.6.0/dune-common-2.6.0.tar.gz
RUN wget https://www.dune-project.org/download/2.6.0/dune-geometry-2.6.0.tar.gz
RUN wget https://www.dune-project.org/download/2.6.0/dune-grid-2.6.0.tar.gz
RUN wget https://gitlab.dune-project.org/dune-fem/dune-fem/-/archive/releases/2.6/dune-fem-releases-2.6.tar.gz
RUN wget https://gitlab.dune-project.org/extensions/dune-alugrid/-/archive/releases/2.6/dune-alugrid-releases-2.6.tar.gz
RUN tar xvzf dune-common-2.6.0.tar.gz
RUN tar xvzf dune-geometry-2.6.0.tar.gz
RUN tar xvzf dune-grid-2.6.0.tar.gz
RUN tar xvzf dune-fem-releases-2.6.tar.gz
RUN tar xvzf dune-alugrid-releases-2.6.tar.gz
