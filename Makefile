CPP          = g++
LA_SOFTWARE  = /home/dahlemd/software

ARPACKPP_DIR = $(LA_SOFTWARE)/arpack++
ARPACKPP_INC = $(ARPACKPP_DIR)/include
SUPERLU_DIR  = $(LA_SOFTWARE)/SuperLU-4.3
BOOST_DIR    = $(LA_SOFTWARE)/boost

LIBDIR = $(LA_SOFTWARE)/lib

ARPACK_LIB   = $(LIBDIR)/libarpack.a
LAPACK_LIB   = $(LIBDIR)/liblapack.a
UMFPACK_LIB  = $(LIBDIR)/libumfpack.5.5.2.a
SUPERLU_LIB  = $(LIBDIR)/libsuperlu_4.3.a
BLAS_LIB     = $(LIBDIR)/libf77blas.a $(LIBDIR)/libatlas.a
FORTRAN_LIBS = -lgfortran

CPP_WARNINGS = -Wall -ansi -pedantic-errors
CPP_DEBUG    = -g
CPP_OPTIM    = -O
CPP_LIBS     = 
CPP_INC      = $(ARPACKPP_DIR)/examples/matrices/sym
BOOST_INC    = ${BOOST_DIR}/include

CPP_FLAGS    = $(CPP_DEBUG) -I$(ARPACKPP_INC) -I$(CPP_INC) -I$(BOOST_INC) $(CPP_WARNINGS)
ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(FORTRAN_LIBS)

all: spectralClustering

spectralClustering: spectralClustering.o
	$(CPP) $(CPP_FLAGS) -o spectralClustering spectralClustering.o $(SUPERLU_LIB) $(ALL_LIBS)

.PHONY: clean
clean:
	rm -f *~ *.o core.* spectralClustering

%.o:    %.cc
	$(CPP) $(CPP_FLAGS) -I$(SUPERLU_DIR) -c $<
