c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRADRHO
c
c
      subroutine GRADRHO( do_spin,
     $ idrhfl0,idrhofl,idrhsfl,igrhofl,
     $ n1r,n2r,n3r, nptr, gprim, ws1,ws2,ws3, ftarray,
     $ rho, rhob, wkr,wkrd )
c---------------------------------------------------------------
c Purpose: compute on grid gradients needed for gga calculation,
c          via analytic reference plus fft residual density
c
c Written: Peter A. Schultz, 13-December-2000, for v2.46 (gga)
c
c Revision history:
c   2Feb08-PAS/2.60e: bug master grho/file
c  27Jun07-PAS/2.60: task-parallel; dist-fftdrho, master-grho
c   1Dec03-PAS/2.58: purge unused argument
c  20Jul01-PAS/2.49: spin-polarized dft
c---------------------------------------------------------------
c
c Notes:
c   I combine all the calls for constructing the (scf) gradients
c   under this rtn for clarity, to make memory mgmt easier, and
c   to make file manipulations easier.
c
c   For spin calculation, it looks like an opportunity to skip
c   gradients if up/dn density empty.  However, need to alter
c   reference calculation to make separate up/dn densities/pots.
c   If empty, why not set to grads zero and ignore other spin?
c   Answer: what about when you have a core density?  Ouch.
c   Hence, a simple skip of gradients is only possible for a
c   very special case: H-atom (or similar core-less single spin).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( one = 1.d0 )
c Input declarations:
      LOGICAL    do_spin
      DIMENSION  gprim(3,3)
c Input array: (delta-density[up/dn] on grid)
      DIMENSION  rho(nptr)
      DIMENSION  rhob(nptr)
c Scratch array:
      DIMENSION  wkr(nptr,7), wkrd(*)
c fft arrays:
      DIMENSION  ws1(*),ws2(*),ws3(*), ftarray(2,*)
c
c Local declarations:
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Get MP info:
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
c
c  Assign memory slots for grid fields (need seven slots total):
c
c gvec(3r),irhog(2r),idrhog(2r)
      igvec  = 1
      irhog  = 4
      idrhog = 6
c
c grho(1r),g2rho(1r),gpggp(1r), drho(3r), ddrho(1r)
      igrho  = 1
      ig2rho = 2
      igpggp = 3
      idrho  = 4
      iddrho = 7
c
c  Make reciprocal space vectors:
c
      call MKGVEC( n1r,n2r,n3r,nptr, gprim, wkr(1,igvec) )
c -->                                       gvec-3o
c
      if( do_spin )then
c       Combine up and dn spin delta-densities ...
        call DAXPY( nptr, one, rhob,1, rho,1 )
      endif
c
c  Compute 1st and 2nd derivatives of delta-density on grid via fft
c  and save to idrhofl
c
      if( iproc.eq.master ) REWIND( unit=idrhofl )
c
      call FFTDRHO( idrhofl,  n1r,n2r,n3r,nptr,
     $ ws1,ws2,ws3,ftarray,
     $ rho,   wkr(1,igvec), wkr(1,irhog), wkr(1,idrhog) )
c -->  rho-io gvec-3i       rhog-2s       drhog-2s
c
c  And, if we have spin-polarized calculation, we need
c  to construct gradients for up/dn separately ...
c
      if( do_spin )then
c
c       Return to separate up/dn delta-densities ...
        call DAXPY( nptr, -one, rhob,1, rho,1 )
c
        if( iproc.eq.master ) REWIND( unit=idrhsfl )
c
c        ... and compute derivatives for up delta-density:
        call FFTDRHO( idrhsfl,  n1r,n2r,n3r,nptr,
     $   ws1,ws2,ws3,ftarray,
     $   rho,   wkr(1,igvec), wkr(1,irhog), wkr(1,idrhog) )
c -->    rho-io gvec-3i       rhog-2s       drhog-2s
c
      endif
c
      if( lstout.gt.2 ) call TIMER( 'After fftdrho' )
c
c     And now we can go and destroy igvec ...
c
c  Combine derivatives from reference density (saved on idrhfl0),
c  and from total delta-density (saved on idrhofl), to construct
c  total density gradient terms needed for gga, and save to igrhofl:
c
c
c Only the master will generate the gradients
      if( iproc .ne. master ) goto 999
c
      REWIND( unit=idrhfl0 )
      REWIND( unit=idrhofl )
c
      call MKGRHO( idrhfl0,idrhofl,   nptr,
     $ wkr(1,igrho),wkr(1,ig2rho),wkr(1,igpggp),
c -->  grho-1o      g2rho-1o      gpggp-1o
     $ wkr(1,idrho),wkr(1,iddrho) )
c -->  drho-3o      ddrho-1s
c
c Shift to beginning of igrhofl:
      REWIND( unit=igrhofl )
c Save total grho, total g2rho, and total gpggp to igrhofl ...
      call WRITBIG( igrhofl, nptr, wkr(1,igrho) )
      call WRITBIG( igrhofl, nptr, wkr(1,ig2rho) )
      call WRITBIG( igrhofl, nptr, wkr(1,igpggp) )
c
      if( do_spin )then
c
c       Need to get gradient terms for up-spin density ...
c
        REWIND( unit=idrhfl0 )
        REWIND( unit=idrhofl )
        REWIND( unit=idrhsfl )
c
        call MKGRHOA( idrhfl0,idrhofl,idrhsfl,   nptr,
     $   wkr(1,igrho),wkr(1,ig2rho),wkr(1,igpggp),
c -->    grho-1o      g2rho-1o      gpggp-1o
     $   wkr(1,idrho),wkr(1,iddrho), wkrd )
c -->    drho-3o      ddrho-1s       grhob-o
c
c       Save grho-up, grho-dn, g2rho-up, and gpggp-up to igrhofl:
        call WRITBIG( igrhofl, nptr, wkr(1,igrho) )
        call WRITBIG( igrhofl, nptr, wkrd )
        call WRITBIG( igrhofl, nptr, wkr(1,ig2rho) )
        call WRITBIG( igrhofl, nptr, wkr(1,igpggp) )
c
c        ... and gradient terms for dn-spin density ...
c
        REWIND( unit=idrhfl0 )
        REWIND( unit=idrhofl )
        REWIND( unit=idrhsfl )
c
        call MKGRHOB( idrhfl0,idrhofl,idrhsfl,   nptr,
     $   wkr(1,igrho),wkr(1,ig2rho),wkr(1,igpggp),
c -->    grho-1o      g2rho-1o      gpggp-1o
     $   wkr(1,idrho),wkr(1,iddrho), wkrd )
c -->    drho-3o      ddrho-1s       wkrd-1s
c
c       Save grho-dn, g2rho-dn, and gpggp-dn to igrhofl:
        call WRITBIG( igrhofl, nptr, wkr(1,igrho) )
        call WRITBIG( igrhofl, nptr, wkr(1,ig2rho) )
        call WRITBIG( igrhofl, nptr, wkr(1,igpggp) )
c
      endif
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
