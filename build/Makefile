#
# Makefile (serial) for SeqQuest code, git-based
# Last updated  4May20-PAS(268e) (modified from JEM/HSE make) for general dev use
#

# The current release version of the code (This must have NO spaces after it)
VERS = d269a

# Dev stuff
qdev = d269
dvrs = d269

###################################
# INSTRUCTIONS FOR BUILDING QUEST #
###################################

#########
#   1   #  Create a new directory (e.g. "obj" or "build") in which to build Quest
#########

#########
#   2   #  Copy this Makefile into this build directory
#########

#########
#   3   #  Modify the following makefile variables:
#########

#########
#  3.0  #  Point QROOT to the Quest root source (git) directory
#########
#### general relative location when in git standard build directory:
# QROOT = ../../PGQuest/
#### personalized absolute locations ...
# QROOT = /home/paschul/lcao/quest/pgquest/
QROOT = /ascldap/users/paschul/lcao/quest/PGQuest/
# QROOT = /ascldap/users/jemouss/quest.HSE

#########
#  3.1  #  COMPILER  - Choose compiler (pick one and only one):
#########
#### f77+mpi:
FC = mpif77
#### pgi compiler:
# FC = /opt/pgi/linux86-64/7.2-1/bin/pgf77
#### hse build
#HSE# FC = mpif77 -g -traceback
#
#### general:
# FC = f77
#### g77:
# FC = g77
#### gfortran (linux/O2, MacOS only tested to O0):
# FC = gfortran

#########
#  3.2  #  FFLAGS  -  Set compiler flags (pick one and only one):
#########
####   chama/intel-mpi
FFLAGS = -O2 -xSSE3 -mcmodel medium -shared-intel
#### general/preferred:
# FFLAGS = -O -fast
#### general/when -fast is not available, or with pgi compiler:
# FFLAGS = -O
#### the linux gfortran requires only tested to O2 (MacOS tested to O0)
# FFLAGS = O2
#### for some g77 type-change (type-casting?) causes compile failure:
# FFLAGS = -O -fno-globals
#### on SGI machines:
# FFLAGS = -O3 -mips4 -OPT:IEEE_arithmetic=3:alias=restrict
#### various outdated parallel
# FFLAGS = -O2 -tpp7 -xW
# FFLAGS = -O2 -tpp7 -xP
# FFLAGS = -O2 -tpp7 -xP -assume 2underscores
# FFLAGS = -O2 -xP
#### hse build
#HSE# FFLAGS = -O

#########
#  3.3  #  LFLAGS  - Linker flags
#########
####  Usually blank for everyone ...
# LFLAGS = 
#### ... except chama/intel
LFLAGS = -mcmodel medium -shared-intel

#########
#  3.4  #  MPI-or-NOT - Platform configuration file
#########
#### general serial code
SUBSU = utl_gen.o
#### general MPI code
SUBSUTP = utl_mpi.o

#########
#  3.5  #  MATH  libraries (to find basics for math and links, lapack/blas, mpi ...)
#########
#
# Select appropriate math library to link in.
# Internally, code default is to use blas/3 routines (flag=doblas3/noblas3)
# WIth an optimized blas library, this is best performance.
# Without an optimized blas library, using blas/3 routines is much slower.
# If your platform does not have an optimized library:
#  1. You myst use the explicit source listed below
#  2. You should edit subroutine "DAT0OPT: to set default "doblas3: to .false.
#
#### pgi compiler: list is likely overcomplete, but it works
# MATH = -L/opt/pgi/linux86-64/7.2-1/lib/ -llapack -lblas -lpgftnrtl -lnspgc -lpgc -lm -lstd -lacml
####  mkl libraries
MKL_LIB = /opt/intel/16.0/mkl/lib/intel64/
MATH = -L${MKL_LIB} -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#### Peter's atlas build on edge
# MATH = -L/home/paschul/local/lapack/lib/ -llapack \
#        -L/home/paschul/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas # -ltstatlas
#### explicit source as provided:
# MATH = lapack.o blas.o
#### dec math libraries:
# MATH = -ldxml
#### compaq renamed math libraries:
# MATH = -lcxml
#### SGI math libraries:
# MATH = -lcomplib.sgimath
#### linux:
# MATH = -llapack -lblas
#### linux/ibm:
# MATH = -static -llapack -lblas
#### Rick's atlas on edge:
# MATH = -L/home/rmuller/Programs/LAPACK/lib/ -llapack \
#        -L/home/rmuller/Programs/ATLAS/lib/ -lf77blas -lcblas -latlas -ltstatlas \
#### various old complex math libraries
# MATH = -L/apps/libraries/l_mkl_p_8.0.1.006/lib/em64t -lmkl -lguide -lmkl_lapack64
# MATH = -L/apps/x86_64/mpi/mvapich/intel-9.1/mvapich-0.9.8-ofed/lib -lmpich \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
# MATH = -L/apps/x86_64/mpi/openmpi/intel-9.1/openmpi-1.1.2-ofed/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
# MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
# MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#        -lmkl -lguide -lmkl_lapack64
## OPENMPI/MKLIB withintegrated scalapack 
# MKL_LIB = /projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t
# MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#       -L/projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t \
#       #       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
# MKL_LIB = /projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t
# MATH = -L/apps/x86_64/mpi/openmpi/intel-11.1-f064-c064/openmpi-1.4.2_static_nomemmgr_oobpr/lib/ -lmpi \
#       -L/projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t \
## note: this needed -liomp5  esolve openmp errors:.
# OMPI_LIB = /opt/openmpi-1.6-intel/lib/
# MKL_LIB = /opt/intel-12.1/mkl/lib/intel64/
# MATH = -L$(OMPI_LIB) -lmpi \
#       -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

#########
#  3.6  #  PMATH  - scalapack/blacs parallel solver libraries
#########
####  serial leave blank
# PMATH = 
####  mkl libraries
PMATH = -L${MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64
####  previous/other scalapack/blacs definitions ...
# PMATHHOME = /apps/libraries/scalapack
#  PBLACSHOME = /apps/libraries/blacs
# PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.2_mx_intel-9.1-f040-c045/scalapack-1.7.4/lib
#  PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.2_mx_intel-9.1-f040-c045/BLACS-1.1/lib
# PMATHHOME = /apps/x86_64/libraries/scalapack/mvapich-0.9.8-ofed-intel-9.1/scalapack-1.7.4/lib
#  PBLACSHOME = /apps/x86_64/libraries/BLACS/mvapich-0.9.8-ofed-intel-9.1/BLACS-1.1/lib
# PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.1.2-ofed-intel-9.1/scalapack-1.7.4/lib
#  PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.1.2-ofed-intel-9.1/BLACS-1.1/lib
# PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.4_ofed_intel-10.0-f025-c025/scalapack-1.7.4/lib
#  PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.4_ofed_intel-10.0-f025-c025/BLACS-1.1/lib
# ####  ... to build a path to the scalapack libraries:
# PMATH =  \
# ${PMATHHOME}/libscalapack.a \
# ${PBLACSHOME}/blacsCinit_MPI-LINUX-0.a \
# ${PBLACSHOME}/blacsF77init_MPI-LINUX-0.a \
# ${PBLACSHOME}/blacs_MPI-LINUX-0.a

#########
#  3.7  #  LAPACK/BLAS - when using explicit source, or to reassign BLAS location
#########
#### - usually this will be blank
# LIBS = 
####  $BLASLIB is part of the MKL module in standard compute cluster environment at Sandia
LIBS = $(BLASLIB)
#### HSE/HEM builds
#HSE# LIBS = $(BLASLIB) -L${QROOT}/libint-quest/lib/.libs/ -lint2 -lstdc++
# LIBS = $(BLASLIB) -L/ascldap/users/jemouss/quest.HSE/libint-quest/lib/.libs/ -lint2 -lstdc++
#### - Compiled blas/lapack (must still list other math libs above)
# MATH = $(LAPO) $(BLASO)
#HSE# MATH = $(LINO)


#########
#   4   #  Type "make" in the obj directory. The program should build
#########  Email any problems to paschul@sandia.gov
## Edit makefile below to select BUILD OPTION

#########################################################################
######   Users should rarely/never change anything below here      ######
#########################################################################

#########################################################################
######                     SOURCE FILES
#########################################################################

## main program
MVRS=lcaomain.o
MVRSTP=lcaomain_tp.o

## program subroutines
SUBS=lcaosubs.o
TPSUBS=lcaosubs_tp.o

######################## Detailed files and breakdown ##################

############ Dev stuff
MDEVTP=lcao_${dvrs}_tp.o
SUBSA=subs/subs.a
TPSUBSA=subs_tp/subs_tp.a
DEVSUBS=subs_${dvrs}.a
DEVSUBA=dev_${qdev}/$(DEVSUBS)
DEVOBJ= setdata.o lmcc.o moment.o ionpot.o iongpot.o zpotout.o
############

SROOT = $(QROOT)/src
MROOT = $(SROOT)/src_math
VPATH = $(SROOT)/src_main:$(SROOT)/src_lib:$(SROOT)/src_subs\
	:$(SROOT)/src_tp\
	:$(MROOT)/linpack:$(MROOT)/lapack:$(MROOT)/blas:$(SROOT)/make

####
####  Math libraries
####

#### LINPACK ROUTINES - main code no longer uses these, but perhaps hse build does?
DIR_LIN = ${QROOT}/src/src_math/linpack
LINF = dsifa.f dsisl.f
LINO = $(LINF:.f=.o)

#### BLAS ROUTINES 
DIR_BLAS = ${QROOT}/src/src_math/blas
BLASF = caxpy.f cdotc.f cgemm.f dasum.f daxpy.f dcabs1.f dcopy.f	\
	ddot.f dgemm.f dgemv.f dger.f dnrm2.f dscal.f dswap.f dsymm.f	\
	dsymv.f dsyr2.f dsyr2k.f dsyrk.f dtrmm.f dtrmv.f dtrsm.f	\
	dtrsv.f dzasum.f dznrm2.f idamax.f izamax.f saxpy.f sdot.f	\
	sgemm.f zaxpy.f zcopy.f zdotc.f zdscal.f zgemm.f zgemv.f	\
	zgerc.f zhemm.f zhemv.f zher2.f zher2k.f zherk.f zscal.f	\
	zswap.f ztrmm.f ztrmv.f ztrsm.f ztrsv.f \
        scopy.f
BLASO = $(BLASF:.f=.o)

#### LAPACK ROUTINES
DIR_LAP = ${QROOT}/src/src_math/lapack
## lapack expert routines (sometimes missing from old std lapack)
LAPXF = dsygvx.f zhegvx.f 
## lapack mainstream routines
LAP0F = dgesv.f dgetf2.f dgetrf.f dgetrs.f dlacpy.f dladiv.f dlae2.f	\
	dlaebz.f dlaev2.f dlagtf.f dlagts.f dlamc1.f dlamc2.f dlamc3.f	\
	dlamc4.f dlamc5.f dlamch.f dlanst.f dlansy.f dlapy2.f dlapy3.f	\
	dlarfb.f dlarf.f dlarfg.f dlarft.f dlarnv.f dlartg.f dlaruv.f	\
	dlascl.f dlaset.f dlasr.f dlasrt.f dlassq.f dlaswp.f dlatrd.f	\
	dorg2l.f dorg2r.f dorgql.f dorgqr.f dorgtr.f dorm2l.f dorm2r.f	\
	dormql.f dormqr.f dormtr.f dpotf2.f dpotrf.f dstebz.f dstein.f	\
	dsteqr.f dsterf.f dsyev.f dsyevx.f dsygs2.f dsygst.f dsygv.f	\
	         dsytd2.f dsytrd.f ieeeck.f ilaenv.f lsame.f xerbla.f	\
	zheev.f zheevx.f zhegs2.f zhegst.f zhegv.f          zhetd2.f	\
	zhetrd.f zlacgv.f zlacpy.f zladiv.f zlanhe.f zlarfb.f zlarf.f	\
	zlarfg.f zlarft.f zlascl.f zlaset.f zlasr.f zlassq.f zlatrd.f	\
	zpotf2.f zpotrf.f zstein.f zsteqr.f zung2l.f zung2r.f zungql.f	\
	zungqr.f zungtr.f zunm2l.f zunm2r.f zunmql.f zunmqr.f zunmtr.f
LAPF = $(LAP0F) $(LAPXF)
LAPO = $(LAPF:.f=.o)

####
#### MAIN QUEST CODE
####

# Main programs
DIR_MAIN = ${QROOT}/src/src_main/
MAINS = COPYRIGHT.NOTICE \
	lcaomain.f lcaomain_tp.f \
	lcao.params lcao_4gb.params lcao_2gb.params lcao_1gb.params \
	utl_gen.f utl_mpi.f

# Library/utility subroutines that might be useful for other things
DIR_LIB = ${QROOT}/src/src_lib/
SUBS_LIB = expund.f mkszero.f mkzero.f \
	spline.f  splint.f  splind.f  splind2.f \
	timstuff.f fftfns.f  strpars.f flstuff.f mpstuff.f strnum.f \
	icopyoff.f bigio.f

# The core code subroutines
DIR_SUB = ${QROOT}/src/src_subs
SUBS_MAIN = ainit.f ainitf.f aofr.f aofrf.f assgn.f atdefct.f atmdat1.f	\
	atmrho.f atomnm.f bfsetup.f bfsplit.f         \
	blkvofn.f boxr.f boxrng.f boxwfn.f broyddj.f bsiall.f		\
	bsicbsv.f bsixgt4.f bulkpol.f c16low.f cee1.f cee.f celldat.f	\
	cellopt.f cellrupd.f cfft3db.f cfft3d.f cinit.f cinitf.f	\
	cmptorb.f cntrrho.f cofr.f cofrf.f config.f convchk.f		\
	corforc.f cross.f dat0ang.f dat0cell.f dat0geom.f	\
	dat0neb.f dat0opt.f dat0run.f datsave.f defbox.f df_gga0.f	\
	df_gga.f dftcgga.f dftget.f dftgga.f dftlda.f dftset.f		\
	dftslda.f dftstr2.f dftstr2sp.f dftstr.f dftxgga.f dippot.f	\
	dipslab.f dmd.f doebind.f doeng1c.f doengy.f dosdump.f		\
	dotrace.f dupemat.f dvstres.f dylm.f ef2dm.f eigsolv.f		\
	eqtrip.f erlsolv.f eslofrc.f evprnt.f          fdefect.f	\
	fermipop.f fftdrho.f          fieldpot.f fintsum.f fixframe.f	\
	          force_elastic.f force_neb.f frc2ctr.f ftorb.f		\
	geomdat.f gesser.f getmijk.f ggastress.f gnorm.f govlp.f	\
	gradrho.f gradwf.f grcofr.f grdovlp.f grelax.f grhodm1b.f	\
	grhodm1.f grhodm2b.f grhodm2.f grhoef1b.f grhoef1.f grhoef2b.f	\
	grhoef2.f grid0set.f grid0xc.f gridave.f gridcgga.f gridchg.f	\
	gridcstr.f griddft.f gridgga.f gridlda.f gridrho.f gridslda.f	\
	gridstr.f gridwf.f gridxgga.f gridxstr.f gselect.f gstamp.f	\
	gtypset.f gupdate.f gvecmag.f h0writ.f            indbox.f	\
	insetup.f intobox.f intsum.f invert.f iongpot.f ionpot.f	\
	kcart.f kcell.f kgrid.f kwrite.f latvec.f lmcc.f	        \
	lmccint.f lvproj.f mat0frc.f mat0set.f mesh.f mkcuts.f		\
	mkeikl.f mkgaunt.f mkgrho.f mkgvec.f mknrwt.f mkrbox.f		\
	mks1atm.f mkshl2d.f                    moment.f          	\
	nearby.f nebassgn.f nebdat.f nebfill.f nebforce.f nebgio.f	\
	nebntrnl.f nebout.f nebreap.f nebrun.f nebsetup.f nlang.f	\
	nlangf.f nlinteq.f nlint.f nlintf.f nlocfrc.f nlocmat.f		\
	nlrad.f nrscut.f onectr.f pop.f pop1k.f popclose.f potnear.f	\
	potnr0.f prodba.f prodhab.f projlv.f projuc.f qhist.f qmd.f	\
	qmddat.f rdstat.f reap.f rh1orb.f rhdmat.f rheigvec.f		\
	rho1ctr.f rhoij.f rintrnl.f rotgrp.f rotmat.f rundata.f		\
	runopt.f schroed.f setdata.f setgr1c.f setorig.f setxc1c.f	\
	sij.f                                      srint.f      	\
	storylm.f str2gpa.f                    symages.f symbuild.f	\
	symfrc.f symgrp.f symibz.f symlat.f symmap2.f symmap3.f		\
	symmat.f symordr.f symprod.f symrho2.f symrho3.f symrho.f	\
	symrmat.f symstr.f symvec.f tangent.f tij.f                     \
	tintsum.f treigvec.f ucatchk.f 	 				\
	v1ctrij.f vesslo.f vlang.f vljisrt.f vlocf0.f vlocf1.f		\
	vlocf2.f vlocfrc.f vlocmat.f vlocmii.f vlocpot.f vlrad.f	\
	vlradf.f vlradz.f vnorms.f vslof1b.f vslof1.f vslof2b.f		\
	vslof2.f vslofix.f vslofrc.f vslom1b.f vslom1.f vslom2b.f	\
	vslom2.f vslomat.f wf1cofk.f wf1cofr.f wkmem.f wrsetd.f		\
	xconmsh.f xcorb.f ylmreal.f ylp2m.f ymm.f yrecur.f              \
        broydata.f broydatf.f ebroyddj.f  rschroed.f rrhoij.f           \
        spinopt.f   matstr.f h0store.f dm0fetch.f                       \
        slabchrg.f zpotout.f   kparallel.f                              \
        banddata.f bandkpts.f bands.f                                   \
        elements.f vdwff.f vdwload.f vdwpair.f vdwstuff.f               \
        eblendasd.f nebstat.f runsteer.f    slics.f slicstuff.f         \
        occdft.f popoccdft.f wreigocc.f
#        focklatt.f fockcell.f fockshl.f fockint.f                       \
#        fockset.f fockdm0.f fockmat.f fockfrc.f fockband.f fockmain.f
SUBSF = $(SUBS_LIB) $(SUBS_MAIN)
SUBS = $(SUBSF:.f=.o)

DIR_TP = ${QROOT}/src/src_tp
TPSUBSF=mptools.f eigstuff.f \
       peigsolv.f perlsolv.f prhoij.f pschroed.f pschroed2.f \
       peigsolvkp.f pschroedkp.f
TPSUBS = $(TPSUBSF:.f=.o)

####################
# precision.f \

#HSE# OBJ = lcaomain.o $(SUBS) $(SUBSU) $(MATH)
OBJ   = lcaomain.o    $(SUBS) $(SUBSU)
OBJTP = lcaomain_tp.o $(TPSUBS) $(SUBS) $(SUBSUTP)

#########################################################################
######                     BUILD OPTIONS
#########################################################################
######   Standard Build Options: Move desired option to front

#### Build Option 1: parallel build
quest: quest_${VERS}.x

#### Build Option 2: serial build
serial: seqquest_$(VERS).x

#### Build Option 3: Dev/dvrs build: move these before quest:/serial: to activate dvrs builds
dev: quest_${dvrs}.x

#### Building free-standing dev subroutine directories, for developing convenience
#### Build Option 3/a: Dev/dvrs build: move these before quest:/serial: to activate dvrs builds
devs: $(DEVSUBS)
#### Build Option 3/b: to build free-standing subroutine object library:
subs: subs.a
#### Build Option 3/c: to build free-standing tp-subroutine object library:
tpsubs: subs_tp.a

seqq-math: seqq-math.x

clean:
	rm -f *.o *.x
#	rm -f *.o *.x *.tgz lcaosubs.f lcaosubs_tp.f Makefile_quest

params:
	rm -f lcao*.o quest_{VERS}.x seqquest_{VERS}.x

help:
	@echo "Make SeqQuest program"
	@echo "Usage: make <platform>"
	@echo "Standard builds: quest (mpi), serial, dev"
	@echo "Serial shortcuts: gnu, pgi, atlas, gfor(+atlas),ibm, sgi, dec, win-x86, linux-x86"
	@echo "You may have to specify location of compiler and blas/lapack libraries"
	@echo "Dev stuff options:"
	@echo "  shipper - make a compact distribution package"
	@echo "  devver  - make a development distribution package"
	@echo "  devs    - build local object library of development subroutines"
	@echo "  subs    - build local object library of serial subroutines"
	@echo "  tpsubs  - build local object library of parallel subroutines"

#########################################################################
######               SERIAL BUILD SHORTCUTS
#########################################################################

gnu:
	$(MAKE) serial   \
		"FC = g77" \
		"FFLAGS = -O -fno-globals" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -llapack_g77 -lf77blas_g77 -lcblas_g77 -latlas_g77 -ltstatlas_g77 -llapack"

pgi:
	$(MAKE) serial   \
		"FC = /opt/pgi/linux86-64/7.2-1/bin/pgf77" \
		"FFLAGS = -O" \
		"LFLAGS =  " \
		"LIBS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -L/opt/pgi/linux86-64/7.2-1/lib/ -llapack -lblas -lpgftnrtl -lnspgc -lpgc -lm -lstd -lacml"

gfor:
	$(MAKE) serial   \
		"FC = gfortran" \
		"FFLAGS = -O" \
		"LFLAGS =  " \
		"LIBS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -L/home/paschul/local/lapack/lib/ -llapack \
		-L/home/paschul/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas"

atlas:
	$(MAKE) serial   \
		"FC = g77" \
		"FFLAGS = -O -fno-globals" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -L/home/rmuller/Programs/LAPACK/lib/ -llapack \
		-L/home/rmuller/Programs/ATLAS/lib/ -lf77blas -lcblas -latlas -ltstatlas"

ibm:
	$(MAKE) serial   \
		"FC = f77" \
		"FFLAGS = -O -fno-globals" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -static -llapack -lblas"

sgi:
	$(MAKE) serial   \
		"FC = f77" \
		"FFLAGS = -O3 -mips4 -OPT:IEEE_arithmetic=3:alias=restrict" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -lcomplib.sgimath"

dec:
	$(MAKE) serial \
		"FC = f77" \
		"FFLAGS = -O -fast" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -ldxml"

win-x86:
	$(MAKE) seqq-math \
		"FC = g77" \
		"FFLAGS = -O -fast" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o"

linux-x86:
	$(MAKE) serial   \
		"FC = f77" \
		"FFLAGS = -O -Wall" \
		"LFLAGS =  " \
		"SUBSU = utl_gen.o" \
		"MATH = -llapack -lblas"

#########################################################################
######                       BUILDS
#########################################################################

quest_$(VERS).x : source_main $(MVRSTP) $(SUBS) $(TPSUBS) $(SUBSUTP)
	$(FC) $(LFLAGS)  -o $@  $(MVRSTP) $(SUBS) $(TPSUBS) $(SUBSUTP)  $(LIBS) $(PMATH) $(MATH)

seqquest_$(VERS).x : source_main ${MVRS} ${SUBS} ${SUBSU}
	$(FC) $(LFLAGS)  -o $@  ${MVRS} ${SUBS} ${SUBSU}  ${MATH} ${LIBS}

seqq-math.x : ${MVRS} ${SUBS} ${SUBSU} lapack.o blas.o
	$(FC) $(LFLAGS)  -o $@  $^


#########################################################################
######                Dev  OPTIONS 
#########################################################################

#subs: subs.a

#tpsubs: subs_tp.a

#dev: quest_${dvrs}.x

quest_$(dvrs).x : $(MDEVTP) $(DEVSUBA) $(TPSUBSA) $(SUBSA) $(SUBSUTP)
	$(FC) $(LFLAGS)  -o $@  $^  $(LIBS) $(PMATH) $(MATH)


#########################################################################
######                SHIP PACKAGES
#########################################################################

#DIST = lcaomain.f lcaomain_tp.f lcaosubs.f lcaosubs_tp.f blas.f lapack.f lcao.params	\
#	Makefile_quest utl_gen.f utl_mpi.f READ.ME
DIST = $(MAINS) lcaosubs.f lcaosubs_tp.f blas.f lapack.f 	\
	Makefile_quest

#DEV_DIST = lcaomain.f lcaomain_tp.f $(SUBS) $(TPSUBS) blas.f lapack.f lcao.params	\
#	Makefile Makefile_quest utl_gen.f utl_mpi.f READ.ME
DEV_DIST = $(MAINS) $(SUBSF) $(TPSUBSF) blas.f lapack.f 	\
	Makefile Makefile_quest

HSE_DIST = lcaomain.f $(SUBSF) $(TPSUBSF) blas.f lapack.f linpack.o lcao.params	\
	Makefile Makefile_quest utl_gen.f utl_mpi.f READ.ME

lcaosubs.f: $(SUBSF)
	cat $^ >> $@

lcaosubs_tp.f: $(TPSUBSF)
	cat $^ >> $@

linpack.f: $(LINF)
	cat $^ >> $@

lapack.f: $(LAPF)
	cat $^ >> $@

blas.f: $(BLASF)
	cat $^ >> $@

#### Build Option 4: create compact user ship package
#shipper: source_main source_subs source_tp source_math $(DIST)
shipper: source_main $(DIST)
	tar czfh quest_$(VERS).tgz $(DIST)
#
#	rm -f  lcaosubs.f lcaosubs_tp.f

#### Build Option 5: create expanded developer ship package
#devver : source_main source_subs source_tp source_math $(DEV_DIST)
devver : source_main source_subs source_tp $(DEV_DIST)
	tar czfh quest_dev_$(VERS).tgz $(DEV_DIST)

#########################################################################
######                DEV-friendly
#########################################################################

## Create object libraries to facilitate separate dev builds
#subs.a: $(SUBS)
$(DEVSUBS): $(DEVOBJ)
	ar ruv  $@ $^
	ranlib  $@

subs.a: source_subs $(SUBS)
	ar ruv  subs.a $(SUBS)
	ranlib  subs.a

#subs_tp.a: $(TPSUBS)
subs_tp.a: source_tp $(TPSUBS)
	ar ruv  subs_tp.a $(TPSUBS)
	ranlib  subs_tp.a

#### to create links to repo source files
#	-for file in ${OBJ_SUB:.o=.f} ; do 
source_main:
	-for file in ${MAINS} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_MAIN}/$$file ;     \
	   fi ;                             \
	 done
source_subs:
	-for file in ${SUBS_MAIN} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_SUB}/$$file ;     \
	   fi ;                             \
	 done
	-for file in ${SUBS_LIB} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_LIB}/$$file ;     \
	   fi ;                             \
	 done
source_tp:
	-for file in ${TPSUBSF} ; do \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_TP}/$$file ;    \
	   fi ;                             \
	 done

source_math: source_lin source_lap source_bla
source_lin:
	-for file in ${LINF} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_LIN}/$$file ;     \
	   fi ;                             \
	 done
source_lap:
	-for file in ${LAPF} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_LAP}/$$file ;     \
	   fi ;                             \
	 done
source_bla:
	-for file in ${BLASF} ; do  \
	   if test ! -f $$file ; then       \
	      ln -s ${DIR_BLAS}/$$file ;     \
	   fi ;                             \
	 done

