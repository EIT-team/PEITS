FROM spack/ubuntu-bionic
LABEL maintainer="r.guichard@ucl.ac.uk"

# Minimal dependencies
RUN apt-get update && \
    apt-get install -y build-essential software-properties-common

# Spack install PETSc and Zoltan
RUN spack install petsc ^parmetis ^metis ^hypre ^superlu-dist
RUN spack install zoltan

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

