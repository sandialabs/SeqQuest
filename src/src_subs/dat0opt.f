c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0OPT
c
c
      subroutine DAT0OPT( lvlout, mkrhopt, ion_opt, ikoopman, kparopt,
     $ dotests, dosetup,doiters,doforce,dorelax,do_cell, do_pvcalc,
     $ opt_spin, do_md, do_neb, do_post,
     $ do_bands,bvl, do_optic, do_nonscf,
     $ dosplit,doblas3, redusym, dolocxc,
     $ dospars,doeigen,doprjct,dorspsy,dodefct,
     $ do_kppsolve, do_psolv )
c---------------------------------------------------------------
c Purpose: initialize options for code
c
c Written: Peter A. Schultz, 23-May-2001, for v2.47
c
c Revision history:
c   1Dec22-PAS/2.69: koopmans/occdft
c   5Aug15-PAS/2.67: initialize "bvl" along with do_bands
c  24Sep13-PAS/2.65a: install flag to compute ehthalpy (E+PV)
c  10Jan13-ACP/2.64: kppsolve flag
c  06Dec10-ACP+AHE/2.63: band structure; 05Dec12-PAS/2.63: psolv flag
c  22Nov08-PAS/2.62: k-parallel flag
c  10Jan08-PAS/2.61: spin optimization
c   7May04-RPM/2.58: add MD flag
c   5Jan02-PAS/2.52: do_cell flag installed
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL  dotests
      LOGICAL  dosetup,doiters,doforce,dorelax,do_cell, do_pvcalc
      LOGICAL  opt_spin, do_md, do_neb, do_post
      LOGICAL  do_bands,do_optic,do_nonscf
      LOGICAL  dosplit,doblas3, redusym, dolocxc
      LOGICAL  dospars,doeigen,doprjct,dorspsy,dodefct
      LOGICAL  do_kppsolve, do_psolv
      CHARACTER*3  bvl
c
c >>>> EXECUTABLE CODE:
c
c Set modest level of output (lvlout: 0=least ... 5=most):
      lvlout = 2
c Grid density scheme (if not 0, 1 compels DM scheme; 2,EF scheme)
      mkrhopt = 0
c Choice for treating charge/moments in supercells:
      ion_opt = 2
c        ABS(ion_opt) = 1: flat jellium countercharge
c        ABS(ion_opt) = 2: local gaussian countermoment
c        ion_opt<0: exclude dipole correction
c        ion_opt>0: include dipole correction
c Koopmans: 0=Off; 1=retrieve saved eigenvectors from reference calculation
      ikoopman = 0
c Use whatever parallelism code decides later (-1=compel std, -2=compel kpar)
      kparopt = 0
C        call SET_KPFLAGS( kparopt )
c By default, attempt  k-parallel-parallel solver (no CTP in subroutines)
      do_kppsolve = .true.
C      do_kppsolve = .false.
c Compel a sequence of all-parallel solves: NO
      do_psolv = .false.
c
c Do not merely test problem setup:
      dotests = .false.
c By default, do just thru scf/iters; no frc, relax, cells, or NEB
      dosetup = .true.
      doiters = .true.
      doforce = .false.
      dorelax = .false.
      do_cell = .false.
      do_pvcalc = .false.
      do_neb  = .false.
      do_md   = .false.
      opt_spin = .false.
c
c output of data for post-processing: OFF
      do_post = .false.
c do non-selfconsistent calculations: OFF
      do_nonscf = .false.
c do bands: OFF, preset Bravais lattice to "nil"
      do_bands=.false.
      bvl = 'nil'
c do optical properties: OFF
      do_optic=.false.
c split basis calculation is: ON
      dosplit = .true.
c blas level 3 use is: ON
      doblas3 = .true.
c allowing symmetry to be reduced: OFF
      redusym = .false.
c
c And some stuff which is not standard usage yet, or not
c yet integrated into production code, and these flags
c should not be altered ...
c
c dolocxc - treat reference atom xc potential in forces is OFF
c   (locxc gives bad forces in current code)
      dolocxc = .false.
c dospars - sparse matrix storage is OFF
      dospars = .false.
c eigensolver is: ON
      doeigen = .true.
c doprjct - project is OFF
      doprjct = .false.
c dorspsy - real space symmetrization is ON
      dorspsy = .true.
c dodefct - defect calculation for linear scaling is OFF
      dodefct = .false.
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
