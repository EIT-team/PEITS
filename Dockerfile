FROM uclrits/cppdev:16.04

USER root
WORKDIR /build
ENV LD_LIBRARY_PATH /usr/local/lib

RUN docker-apt-install libblas-dev liblapack-dev gfortran libboost-dev \
    mpi-default-bin mpi-default-dev libgmp-dev libopenmpi-dev gfortran \
    libmetis-dev libparmetis-dev autotools-dev automake libug-dev

########## Build PETSC ##########
RUN git clone -b maint https://gitlab.com/petsc/petsc.git && \
    cd petsc && git checkout 8695de0 && \
    ./configure --prefix=/usr \
                --with-x=0 --with-debugging=0 \
                CFLAGS="-O3 -DNDEBUG -ffast-math" \
                --with-parmetis=1 \
                --with-metis=1 \
                --with-hypre=1 --download-hypre=yes \
                --with-superlu_dist=1 --download-superlu_dist=yes \
                --with-mumps=1 --download-mumps=yes \
                --with-ml=1 --download-ml=yes \
                --download-scalapack=yes --download-blacs=yes && \
    make all && make install

########## Build Zoltan ##########
# NOTE - Because Zoltan doesn't allow to build on the source directory, it needs to go into zoltanBUILD
#      - why does it need the path to incdir and libdir? it doesn't know how to look in the right place?
RUN wget http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz && \
    tar xf zoltan_distrib_v3.83.tar.gz && \
    mkdir Zoltan_v3.83/zoltanBUILD && cd Zoltan_v3.83/zoltanBUILD && \
    ../configure --prefix=/usr \
                 --with-parmetis -with-parmetis-incdir=/usr/include/ --with-parmetis-libdir=/usr/lib/ && \
    make everything && \
    make install

RUN rm -rf {petsc,?oltan}

# Copy the local files
COPY . /build/PEITS/
WORKDIR PEITS

# Configure PEITS build
RUN cp config.opts_example config.opts && \
    sed -i -e '/^PETSCPATH/d' -e '/ZOLTANPATH/d' \
        -e '/with-metis/d' -e '/with-parmetis/d' \
        config.opts && \
    sed -i -e '1 a PETSCPATH="/usr/lib"' \
        config.opts

# Build PEITS
RUN bash ./INSTALL.sh

# Create User and set environment
RUN useradd -ms /bin/bash peitsier
RUN chown -R peitsier:peitsier /build/PEITS
USER peitsier
ENV HOME /home/peitsier
RUN cp /etc/zsh/newuser.zshrc.recommended /home/peitsier/.zshrc
WORKDIR /mydata
CMD ["bash"]
