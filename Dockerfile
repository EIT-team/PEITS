FROM spack/ubuntu-bionic
LABEL maintainer="r.guichard@ucl.ac.uk"

# Minimal dependencies
RUN apt-get update && \
    apt-get install -y build-essential software-properties-common

# Spack install PETSc and Zoltan
RUN spack install petsc ^parmetis ^metis ^hypre ^superlu-dist
RUN spack install zoltan
RUN spack install boost

# Get and untar the Dune modules
WORKDIR /home
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

# Append bashrc to load the necessary modules for PEITS build
RUN echo 'module load openmpi-3.1.4-gcc-7.4.0-asctus2' >> ~/.bashrc
RUN echo 'module load petsc-3.11.3-gcc-7.4.0-ctafkbm' >> ~/.bashrc
RUN echo 'module load zoltan-3.83-gcc-7.4.0-cmexe3l' >> ~/.bashrc
RUN echo 'module load pkgconf-1.6.1-gcc-7.4.0-j566ycz' >> ~/.bashrc
RUN echo 'module load cmake-3.15.1-gcc-7.4.0-4dvxaym' >> ~/.bashrc
RUN echo 'module load boost-1.70.0-gcc-7.4.0-d42gtzk' >> ~/.bashrc