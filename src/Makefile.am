GRIDTYPE = ALUGRID_CONFORM
GRIDDIM = 3
WORLDDIM = 3
POLORDER = 1
USE_ZOLTANLB = 1
WANT_ISTL = 0
WANT_PETSC = 1

PETSCADD = -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
 -lparmetis -lmetis -lpthread -lHYPRE \
 -lmpi_cxx -lstdc++ -lscalapack -lml -lmpi_cxx \
 -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl \
 -lstdc++ -llapack -lblas -lgfortran -lm -lm -lm -lm -lmpi_cxx \
 -lstdc++ -lmpi_cxx -lstdc++ -ldl -lmpi 

SUBDIRS =

noinst_PROGRAMS = dune_peits

dune_peits_SOURCES = dune_peits.cc

dune_peits_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) \
  $(ZOLTAN_CPPFLAGS) \
  $(GRIDDIM_CPPFLAGS) \
	$(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) \
  -DWANT_ISTL=$(WANT_ISTL) -DWANT_PETSC=$(WANT_PETSC) 
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.

dune_peits_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD) \
	$(ALL_PKG_LIBS) $(PETSCADD)

dune_peits_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS) \
	$(ALL_PKG_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-alugrid=$(DUNE_ALUGRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

