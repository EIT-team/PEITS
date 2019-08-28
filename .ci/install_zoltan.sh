#!/bin/bash -xe

if [ ! -f ${ci_path}/zoltan_cached ]; then
    wget http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz
    tar xf zoltan_distrib_v3.83.tar.gz --warning=no-unknown-keyword
    mkdir Zoltan_v3.83/zoltanBUILD
    pushd Zoltan_v3.83/zoltanBUILD
    ../configure --prefix=${ci_path}/usr \
                 --with-parmetis --with-parmetis-incdir=/usr/include/ --with-parmetis-libdir=/usr/lib/
    make everything
    make install
    popd
    touch zoltan_cached
fi
