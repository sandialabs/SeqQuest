#
##
## Makefile for the SeqQuest code, MPI-PARALLEL
##
## Last revised: 27Jun18-PAS for version 2.67c (git-conversion)
## Configured for parallel: openmpi, mkl with integrated scalapack/blacs
##
#

# directory with changes
qdev = d268
# version of changes (main & subs library)
qvrs = d268b

SRCROOT=../src/
VPATH = $(SRCROOT)/src_main

######
###### Choose compiler (pick one and only one):
######
## general:
# FC = f77
## g77:
# FC = g77
## gfortran (linux/O2, MacOS only tested to O0):
# FC = gfortran
## pgi compiler:
# FC = /opt/pgi/linux86-64/7.2-1/bin/pgf77
## f77+mpi:
FC = mpif77

######
###### Set compiler flags (pick one and only one):
######
## general/preferred:
# FFLAGS = -O -fast
## general/when -fast is not available, or with pgi compiler:
# FFLAGS = -O
## for some g77 type-change (type-casting?) causes compile failure:
# FFLAGS = -O -fno-globals
## the linux gfortran requires only tested to O2
# FFLGAS = O2
## the non-optimized version is safest (and slowest), test for MacOS/gfortran
# FFLGAS = O0
## on SGI machines:
# FFLAGS = -O3 -mips4 -OPT:IEEE_arithmetic=3:alias=restrict

## ICC/liberty:
# FFLAGS = -O2 -tpp7 -xW
## NWCC/spirit-old and tbird-old:
# FFLAGS = -O2 -tpp7 -xP
## NWCC/spirit-openmpi
# FFLAGS = -O2 -tpp7 -xP -assume 2underscores
## NWCC/tbird-mva-0.9.8
# FFLAGS = -O2 -xP
## NWCC/tbird-openmpi-1.1.2
# FFLAGS = -O2 -tpp7 -xP -assume 2underscores
## NWCC/tbird-openmpi-1.2.4 (-tpp7 & -assume 2underscores retired)
# FFLAGS = -O2 -xP
## redsky
FFLAGS = -O2
## chama
FFLAGS = -O2 -xSSE3 -mcmodel medium -shared-intel

## everyone ...
# LFLAGS =
## ... except chama:
LFLAGS = -mcmodel medium -shared-intel
LIBS =

######
###### Math libraries (to find basics for math and links, lapack/blas, mpi ...)
######
#
# Select appropriate math library to link in.
# Internally, code default is to make use of blas/3 routines,
# an option that is flag-controlled (doblas3).  With an
# optimized blas library, this gives the best performance,
# but without an optimized library, using blas/3 routines
# gives much slower execution.  If your platform does not
# have an optimized library and you must use the explicit
# source provided below, edit subroutine "dat0opt" in the
# file lcaosubs.f to set default "doblas3" to .false. for
# most efficient execution.
#
## dec math libraries:
# MATH = -ldxml
## compaq renamed math libraries:
# MATH = -lcxml
## SGI math libraries:
# MATH = -lcomplib.sgimath
## linux:
# MATH = -llapack -lblas
## linux/ibm:
# MATH = -static -llapack -lblas
## explicit source as provided:
# MATH = lapack.o blas.o
## Rick's atlas on edge:
# MATH = -L/home/rmuller/Programs/LAPACK/lib/ -llapack \
#        -L/home/rmuller/Programs/ATLAS/lib/ -lf77blas -lcblas -latlas -ltstatlas \
## pgi compiler: likely overcomplete, but it works
# MATH = -L/opt/pgi/linux86-64/7.2-1/lib/ -llapack -lblas -lpgftnrtl -lnspgc -lpgc -lm -lstd -lacml

###### Math libraries for parallel

## MPI link (simple):
# MATH = lapack.o blas.o -lmpi
## ICC/liberty math libraries:
## MATH = -L/apps/intel/mkl70/lib/32 -lmkl_ia32 -lguide -lmkl_lapack
# MATH = -L/apps/intel/mkl70/lib/32 -lmkl -lmkl_ia32 -lguide -lmkl_lapack \
#       -L/apps/intel/fce/lib -lsvml
## NWCC/spirit-old math libraries:
# MATH = -L/apps/intel/mkl701/lib/em64t -lmkl -lguide -lmkl_lapack64
## NWCC/spirit-openmpi math libraries:
# MATH = -L/apps/x86_64/mpi/openmpi/intel-9.1-f040-c045/openmpi-1.2.2_mx/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_9.1.018/lib/em64t \
#        -lmkl -lguide -lmkl_lapack

## tbird parallel math libraries
## tbird-old
# MATH = -L/apps/libraries/l_mkl_p_8.0.1.006/lib/em64t -lmkl -lguide -lmkl_lapack64
## tbird-mva-0.9.8
# MATH = -L/apps/x86_64/mpi/mvapich/intel-9.1/mvapich-0.9.8-ofed/lib -lmpich \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
## tbird-openmpi-1.1.2
# MATH = -L/apps/x86_64/mpi/openmpi/intel-9.1/openmpi-1.1.2-ofed/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
##
## >>>>> openmpi-1.2.4/mkl8 (had problems with mkl9), with separate scalapack
##
#MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#       -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#       -lmkl -lguide -lmkl_lapack64
##
## >>>>> OPENMPI with MKL10 with integrated scalapack 
##
#MKL_LIB = /projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t
#MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#       -L/projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t \
#       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
##
## >>>>> MKL libraries: redsky
##
### old 11.1 version, retired 4Jun14 (after big redsky TOSS2/12.1-compiler upgrade)
#MKL_LIB = /projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t
#MATH = -L/apps/x86_64/mpi/openmpi/intel-11.1-f064-c064/openmpi-1.4.2_static_nomemmgr_oobpr/lib/ -lmpi \
#       -L/projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t \
#       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
### new as of 6/4/2014 (needed to add -liomp5 to resolve openmp errors).
#OMPI_LIB = /opt/openmpi-1.6-intel/lib/
#MATH = -L$(OMPI_LIB) -lmpi \
#       -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
##
## >>>>> MKL libraries: chama
##
MKL_LIB = /opt/intel/16.0/mkl/lib/intel64/
MATH = -L${MKL_LIB} -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5


######
###### scalapack parallel math libraries
######
##
####  First, define patch to scalapack and blacs libraries
##
## tbird-old
# PMATHHOME = /apps/libraries/scalapack
# PBLACSHOME = /apps/libraries/blacs
## NWCC/spirit-openmpi
# PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.2_mx_intel-9.1-f040-c045/scalapack-1.7.4/lib
# PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.2_mx_intel-9.1-f040-c045/BLACS-1.1/lib
## tbird: mva-0.9.8
#PMATHHOME = /apps/x86_64/libraries/scalapack/mvapich-0.9.8-ofed-intel-9.1/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/mvapich-0.9.8-ofed-intel-9.1/BLACS-1.1/lib
## tbird: openmpi-1.1.2
#PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.1.2-ofed-intel-9.1/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.1.2-ofed-intel-9.1/BLACS-1.1/lib
## tbird: openmpi-1.2.4/mkl8 (had problems with mkl9)
#PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.4_ofed_intel-10.0-f025-c025/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.4_ofed_intel-10.0-f025-c025/BLACS-1.1/lib
##
#### >>>>> build pmath to libraries
##
#PMATH =  \
#${PMATHHOME}/libscalapack.a \
#${PBLACSHOME}/blacsCinit_MPI-LINUX-0.a \
#${PBLACSHOME}/blacsF77init_MPI-LINUX-0.a \
#${PBLACSHOME}/blacs_MPI-LINUX-0.a
##
## >>>>> MKL with integrated scalapack libraries, redsky and chama
##
PMATH = -L${MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

#########################################################################
#########################################################################
#########################################################################

######
###### Platform-dependent configuration file (serial or mpi-parallel)
######

## general serial code
# SUBSU = utl_gen.o
## general MPI code
SUBSU = utl_mpi.o

######
###### Source files
######

## main program
#MVRS=lcaomain.o
#MVRSTP=lcaomain_tp.o
MDEVTP=lcao_${qvrs}_tp.o
MAINTP=lcaomain_tp.o

## program subroutines
#SUBS=lcaosubs.o
#TPSUBS=lcaosubs_tp.o
DEVSUBS=dev_${qdev}/subs_${qvrs}.a
SUBS=subs/subs.a
TPSUBS=subs_tp/subs_tp.a

######
###### Where it all actually happens
######

dev: quest_${qvrs}.x
quest: quest_tp.x


help:
	@echo "Make SeqQuest program - openmpi parallel"
	@echo "Usage:"
	@echo "  make quest - to get production git version"
	@echo "  make dev   - to build local dev version"
	@echo "  make clean - to clear out .o, .x files"
	@echo "For parallel builds, you will have to specify library locations in makefile"
	@echo ""

clean:
	rm -f *.o quest*.x

params: quest_${qvrs}.x
	rm -f ${MVRSTP} quest_${qvrs}.x
	make quest

######
###### The actual build
######

quest_${qvrs}.x : ${MDEVTP} ${SUBS} ${TPSUBS} ${SUBSU} ${DEVSUBS}
	$(FC) $(LFLAGS) -o $@  ${MDEVTP} ${DEVSUBS} ${SUBSU} ${TPSUBS} ${SUBS} ${PMATH} ${MATH}

quest_tp.x : ${MAINTP} ${SUBS} ${TPSUBS} ${SUBSU}
	$(FC) $(LFLAGS) -o $@  ${MAINTP} ${SUBSU} ${TPSUBS} ${SUBS} ${PMATH} ${MATH}
#	$(FC) $(LFLAGS) -o $@  ${MVRSTP} ${SUBSU} ${TPSUBS} ${SUBS} ${PMATH} ${MATH}
