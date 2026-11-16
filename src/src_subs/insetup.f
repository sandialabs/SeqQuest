c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INSETUP
c
c
      subroutine INSETUP( IWR, dokdefct,
     $ natm,ntyp, norb,norbd, nstbulk,nedefct,
     $ elchrg, nocc0k, do_spin,opt_spin,spinpol, nspin,elecno,spinfac,
     $ itypa, ioffa, norba, znuc )
c---------------------------------------------------------------
c Purpose: process some system information
c
c Written: Peter A. Schultz,  9-May-2000, for v2.45
c
c Revision history:
c  22May13-JEM/2.65: fix check of # electrons > zero
c  10Jan08-PAS/2.61: do not pass unused argument (elnumbr); spin-opt
c  15Feb05-PAS/2.59: defect k-sampling
c   8Aug02-PAS/2.53(f): elnum(spin)
c  28Jun02-PAS/2.53: write out number of basis functions
c  21Jun01-PAS/2.48: replace STOPs, clean a little
c   7Jun01-PAS/2.47: strip out nstd dimension check
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER( zero=0.d0, half=0.5d0, one=1.0d0, two=2.d0 )
c
c input/output arrays:
      LOGICAL    dokdefct
      DIMENSION  nstbulk(2), nedefct(2)
      DIMENSION  itypa(natm), ioffa(natm)
      DIMENSION  norba(ntyp)
      DIMENSION  znuc(ntyp)
      LOGICAL    do_spin, opt_spin
      DIMENSION  elecno(2)
c
c >>>> EXECUTABLE CODE:
c
 9020 format(1x,a,3f20.10)
C      any1ctr = .true.
C      anycore = .true.
C      anyxc1c = .true.
      norb = 0
      elnumbr = -elchrg
      do  iatm=1,natm
        ityp = itypa(iatm)
        ioffa(iatm) = norb
        norb = norb + norba(ityp)
        elnumbr = elnumbr + znuc(ityp)
C        if( lmx1ctr(ityp).ge.0 ) any1ctr = .true.
C        if( nrcor(ityp).ne.0 ) anycore = .true.
C        if( nfityp(ityp).gt.0 ) anyxc1c = .true.
      enddo
c
c Check, and set up electron counts
      write(IWR,9020)   '>>>> Total number of electrons=       ',elnumbr
C      if( ABS( elchrg ) .gt. elnumbr ) call STOPXERR( 'charge > # elec')
      if( elnumbr .lt. zero ) call STOPXERR( '# elec < 0')
      if( do_spin )then
        if( spinpol .gt. elnumbr ) call STOPXERR( 'spinpol > # elec' )
        nspin = 2
        spinfac = half
        elecno(1) = ( elnumbr + spinpol ) / two
        elecno(2) = elnumbr - elecno(1)
        elstate = (elecno(1) + 0.9999999d0)
        write(IWR,9020) '>>>> Spin!  Number of up/dn electrons=',elecno
        if( dokdefct )then
          nupstate = elstate
          eldown = (elecno(2) + 0.9999999d0)
          ndnstate = eldown
          nstbulk(1) = nupstate - nedefct(1)
          nstbulk(2) = ndnstate - nedefct(2)
          write(IWR,9020) 'Defect sampling, bulk states, up/dn=',nedefct
        endif
      else
        opt_spin = .false.
        nspin = 1
        spinfac = one
        elecno(1) = elnumbr
        elecno(2) = zero
        elstate = (elnumbr   + 1.9999999d0) / two
        if( dokdefct )then
          nupstate = elstate
          nstbulk(1) = nupstate - nedefct(1)
          write(IWR,9020) 'Defect sampling, bulk states=',nedefct(1)
        endif
      endif
c
c Minimum number of occupied states
      nocc0k = elstate
c
      write(IWR,*) '>>>> Minimum number of occupied states=',nocc0k
      write(IWR,*) '>>>> Number of basis functions, norb=',norb
      if( norb.gt.norbd )then
        write(IWR,*) '***** max norb dimension norbd =',norbd
        call STOPXERR( 'norb2big/norb is larger than dimension' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
