c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0RUN
c
c
      subroutine DAT0RUN( doforce, do_post, dokdefct,icloseocc,
     $ ikoopman,ioccdft,
     $ itstart,itstop,meth_bl,nhiste,igesopt,  nocc0k,norb,nstate,
     $ convgr,convs,convsl,convii, conv1c, alfast,
     $ opt_spin, meth_sp,conv_sp,nstep_sp,blend_sp,istart_sp,
     $ edegen,etemp, scfblnd,scfbl2,scfconv, idosolv, scut )
c---------------------------------------------------------------
c Purpose: initialize parameters for run (scf) phase
c
c Written: Peter A. Schultz, 23-May-2001, for v2.47
c
c Revision history:
c   1Dev22-PAS/2.69: Koopmans/occdft capability
c  25Mar16-PAS/2.67c: nhist default to 4 (from 6)
c  14Dec14-PAS/2.66: install ASD e-blend option
c  11Mar13-PAS/2.65: change nhiste from 10 to 6, itstop fr. 30 to 50
c  10Jan08-PAS/2.61: spin optimization
c  14Jun05-RPM/2.59: added scut variable to set cutoff for
c                    small eigenvalues of the overlap matrix
c  14Feb05-PAS/2.59: defect k-sampling
c  25Nov03-PAS/2.58: remove scf blend reset feature
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c   9Jan03-PAS/2.55: blending changes, and cosmetics
c   9Aug02-PAS/2.53(f): fix nocc0k for spin-pol calculation
c  16May02-PAS/2.53: cutoff params moved to "rundata"
c   8Mar02-PAS/2.52: literal constants to d0
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      LOGICAL    doforce, do_post, dokdefct, opt_spin
c Output:
      DIMENSION  convs(8)
c
c Local:
      DIMENSION  convs0(8)
      DATA  convs0 / 3.0d-9 , 3.0d-9 ,1.5d-7 , 6.0d-8,
     $               3.0d-7 , 3.0d-7 ,0.5d-5 , 6.0d-6 /
c Defines blending type
      DATA    methbl_broy,methbl_asd / 1,2 /
c
c >>>> EXECUTABLE CODE:
c
c itstart,itstop = initial/final iteration to be max:
      itstart = 0
      itstop = 50
c     If koopmans set, default to stop after koopmans, before scf
      if( ikoopman .gt. 0 ) itstop = 1
c ioccdft = occ-constrained dft, defaults to OFF=0
c   ioccdft = 0 OFF
c   ioccdft = 1 excite-generated with nk=1, to be expanded to full nk
c   ioccdft - 2 from file with presumably full nk
      ioccdft = 0
c meth_bl = blending type: 1=broyden (default), 2=asd (alternate)
      meth_bl = methbl_broy
c nhiste = max number of steps in scf blend
      nhiste =  4
c igesopt = scf ges: 0=raw,1=full d-Ham,2=true d-Ham,3=d-atHam
      igesopt = -2
c
c Default number of occupied states ...
      nstate = MIN( nocc0k+10, norb )
c      ... if post-process, and not a force calc, do all states
c       (at AFRL request)
      if( do_post .and. .not. doforce ) nstate = norb
c
c icloseocc = 0/normal fermi occs; 1/binned by k; 2/0K (binned)
      icloseocc = 0
      if( dokdefct ) icloseocc = 1
c
c convgr = cutoff value for range of grid orbitals/shells:
      convgr = 2.d-8
       call FLGETIWR( IWR )
      write(IWR,*) '************** VERSION 2.65 CHANGE: ***************'
      write(IWR,*) '**** grid range default changed from 2.d-7 to 2.d-8'
      write(IWR,*) '**** This will cause numerical changes in results'
      write(IWR,*) '**** To recover pre-2.65, use input cutgrd/2.d-7'
      write(IWR,*) '***************************************************'
c convs() = cutoffs for the various (semi-)analytic routines:
c         setup: (1) ovlp (2) kin (3) nloc (4) local
c        &force: (5) ovlp (6) kin (7) nloc (8) local
      do  ic=1,8
        convs(ic) = convs0(ic)
      enddo
c convsl = cutoff for pairs of slow gaussians
      convsl = 1.d-1
c convii = cutoff Vnuc:II (bf on same atom
      convii = 1.d-2
c conv1c = cutoff for scf fast potential integrals (non-std)
      conv1c = 1.d-4
c
c alfast = Gaussian exponent defining "slow/fast"
      alfast = 0.6d0
c
c etemp = Temperature for occupation weighting
      etemp = 0.003d0
c edegen = criterion for degeneracy of two states for pop.
      edegen = 0.0002
c scfblnd = Initial blend factor for Hamiltonian
      scfblnd = 0.30d0
c scfbl2 = Initial blend factor for Ham after first geometry
      scfbl2 = 0.30d0
c scfconv =  SCF convergence criterion:
      scfconv = 0.00030d0
c
c  idosolv = 1 full eigensolve, using dsygv/zhegv lapack routines
c  idosolv = 2 uses faster dsygvx/zhegv selected eigenvec routines
c  idosolv = 3 uses dsyev with a routine to remove basis linear dependence
      idosolv = 2
      if (idosolv .eq. 3 .or. idosolv .eq. 1) nstate = norb
c
c  scut = Cutoff below which eigenpairs of S are ignored if idosolv=3
      scut = 1d-6
c
c Spin optimization:
c  meth_sp = method for spin optimization (0=simple,1=linear)
      meth_sp = 1
c  conv_sp = spin convergence criterions (polarization, in e)
      conv_sp = 0.004d0
c  blend_sp = slope for updating spin polarization (should be 1.0-2.0)
      blend_sp = 1.5d0
c  nstep_sp = maximum number of steps for spin convergence
      nstep_sp = 10
c  istart_sp = initial start (do not change from 0)
      istart_sp = 0
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
c End
