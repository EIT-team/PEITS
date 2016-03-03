# $Id: duneproject 7100 2013-02-02 03:11:29Z gruenich $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = src m4 dune data doc

if BUILD_DOCS
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.9 or newer
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-alugrid=$(DUNE_ALUGRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

