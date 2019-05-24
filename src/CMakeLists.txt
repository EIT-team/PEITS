# Some variables that need to be passed for compilation
# and grid definition
if( NOT POLORDER )
  set( POLORDER 1)
endif()

if( NOT GRIDTYPE )
  set( GRIDTYPE ALUGRID_CONFORM )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 3 )
endif()

if( NOT WORLDDIM )
  set( WORLDDIM 3 )
endif()

# the default compiler flags
set( DEFAULTFLAGS "GRIDDIM=${GRIDDIM};${GRIDTYPE};WORLDDIM=${WORLDDIM};POLORDER=${POLORDER}" )

add_executable("dune-peits" dune-peits.cc)
target_compile_definitions("dune-peits" PUBLIC "${DEFAULTFLAGS}")
# target_compile_definitions("dune-peits" PUBLIC ALUGRID_CONFORM "GRIDDIM=3" "WORLDDIM=3" -DPOLORDER=${POLORDER})
# add_definitions(-DGRIDTYPE=${GRIDTYPE} -DGRIDDIM=${GRIDDIM} -DWORLDDIM=${WORLDDIM} -DPOLORDER=${POLORDER})
target_link_dune_default_libraries("dune-peits")
