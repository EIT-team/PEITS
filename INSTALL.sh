#!/bin/bash

function dl_dune () {
    # usage: dl_dune <url>

    filename=$(basename -- $1)
    pattern='dune-[[:alpha:]]+'
    [[ "${filename}" =~ $pattern ]]
    project=${BASH_REMATCH[0]}
    curl -O $1
    mkdir $project
    tar xzf ${filename} -C $project --strip-components=1
}

#########   Dune packages ##########
# We need to use the local copy of dune. It's unknown which exact version they
# are. Below we point to the git commit that they are closest to (but that breaks
# everything!)
####################################

######## Dune common #######
# This package is from before the 2.3 release, closest commit is: a29299b4
# That has 8 differences between the equal files, a total of 50 files that
# exists on one but not on the other.
# Also tried with the 2.3.0 released and failed
# dl_dune https://dune-project.org/download/2.3.0/dune-common-2.3.0.tar.gz
tar xvzf dune-common-2.3-svn.tar.gz

######## Dune alugrid  #######
# Attempted with release provided by Dune but failed
#  dl_dune https://gitlab.dune-project.org/extensions/dune-alugrid/-/archive/releases/2.3/dune-alugrid-releases-2.3.tar.gz
tar xvzf dune-alugrid-2.3.tar.gz

######## Dune fem  #######
# Attempted with closer commit: ff22d61a6037947c997c0b0bc034922e529063b4
# but failed.
# dl_dune https://gitlab.dune-project.org/dune-fem/dune-fem/-/archive/ff22d61a6037947c997c0b0bc034922e529063b4/dune-fem-ff22d61a6037947c997c0b0bc034922e529063b4.tar.gz
tar xvzf dune-fem-1.4.0.tar.gz

######### Dune geometry ########
# Attempted with release but failed
# dl_dune https://dune-project.org/download/2.3.0/dune-geometry-2.3.0.tar.gz
tar xvzf dune-geometry-2.3-svn.tar.gz

######### Dune grid ##########
# Attempted with release but failed
# dl_dune https://dune-project.org/download/2.3.0/dune-grid-2.3.0.tar.gz
tar xvzf dune-grid-2.3-svn.tar.gz

./dune-common-2.3-svn/bin/dunecontrol --opts=config.opts all
