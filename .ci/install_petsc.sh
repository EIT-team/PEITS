#!/bin/bash -xe

if [ ! -f ${ci_path}/petsc_cached ]; then
    git clone -b maint https://gitlab.com/petsc/petsc.git
    pushd petsc
    git checkout 8695de0
    ./configure --prefix=${ci_path}/usr \
                --with-x=0 --with-debugging=0 \
                CFLAGS="-O3 -DNDEBUG -ffast-math" \
                --with-parmetis=1 \
                --with-metis=1 \
                --with-hypre=1 --download-hypre=yes \
                --with-superlu_dist=1 --download-superlu_dist=yes \
                --with-mumps=1 --download-mumps=yes \
                --with-ml=1 --download-ml=yes \
                --download-scalapack=yes --download-blacs=yes
    make
    make install
    popd
    touch petsc_cached
fi
