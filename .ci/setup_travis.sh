#!/bin/bash -xe


ci_path=/build

# Check whether cache exists
if [ ! -d ${ci_path} ]; then
    mkdir -p ${ci_path}
fi

pushd ${ci_path}

# Set up petsc
source $TRAVIS_BUILD_DIR/.ci/install_petsc.sh
# Set up zoltan
source $TRAVIS_BUILD_DIR/.ci/install_zoltan.sh

popd
