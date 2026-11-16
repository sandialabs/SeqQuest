c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POPOCCDFT
c
c
      subroutine POPOCCDFT( nstate,nk,nspin,
     $ elnumbr,npop, occfull,
     $ efermi,egap, esum1e,
C     $ etemp,edegen,
     $ wtk, eigval, eigpop )
c---------------------------------------------------------------
c Purpose: populated 1e-states using constrained-occ-dft occupations
c
c Written: Peter A. Schultz, 10-December-2022, for v2.69 (koopmans/occdft)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c   This sets (in principle, non-Aufbau) occupations of 1e-states
c   We obtain the occ-dft state populations from a file (unit=iflcdft)
c     ioccdft = 1 ... only have single k-pt, expand out to all nk
c     ioccdft = 2 ... full nstate*nk from saved file
c
c Notes:
c   Need to define efermi,egap (ill-defined in general occdft):
c    efermi = lowest (across k-pts) of highest occ'd state
c    egap   = range bottom to top (across k-pts) of this highest state
c   We ignore "dokdecft" G-sampled defect state stuff (~stale anyway)
c   For now, ignore etemp/edegen, install as strictly closed shell
c    - likely need this later, if near-degens, or to deal with cbe/vbe things
c    
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.d0 )
      PARAMETER ( idebug=3 )
c
c Input:
C      DIMENSION  elecno(*)
      DIMENSION  wtk(nk)
      DIMENSION  eigval(nstate,nk)
      DIMENSION  eigpop(nstate,nk)
c Output:
C      DIMENSION  efermi(*), egap(*)
C      DIMENSION  npops(*)
c Local:
      DATA   delocc / 1.d-8 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
      if( idebug .gt. 0 )then
        write(IWR,*) 'POPOCCDFT: set up occDFT occupations'
      endif
c
c We have the occupations upon input already
c If we have no electrons, skip everything
      if( elnumbr .le. zero )then
        npop = 0
        esum1k = zero
        gap = -one
        efermi = eigval(1,1)
        goto 900
      endif
c
c Count occupied states
c
      npop0gs = ( elnumbr+occfull-delocc ) / occfull
      if( npop0gs .gt. nstate )then
        write(IWR,*) 'popoccdft: npop0gs,nstate=',npop0gs,nstate
        call STOPXERR( 'POPOCCDFT: fewer states than occups' )
      endif
c
c     Assume pop always non-negative, assume pop's make sense
      esum1e = zero
      do 200 ik=1,nk
        pop1k  = zero
        esum1k = zero
        do 110 n=1,nstate
          popnk = eigpop(n,ik)
          eignk = eigval(n,ik)
          if( popnk.gt.zero )then
            npop1k = n
            pop1k  = pop1k  + popnk
            esum1k = esum1k + popnk*eignk
          endif
  110   continue
        esum1e = esum1e + esum1k*wtk(ik)
c       Do an idiot-check
        if( ABS( elnumbr-pop1k ) .gt. delocc )then
          write(IWR,*) 'popoccdft: pops .ne. electron #, k=',ik
          write(IWR,*) ' elnumbr(target), pops(occs)=',elnumbr, pop1k
          write(IWR,'(20f5.2)') (eigpop(n,ik),n=1,nstate)
          call STOPXERR( 'POPOCCDFT: populations wrong' )
        endif
c
        ehomo1k = eigval(npop1k,ik)
        elumo1k = eigval(npop1k+1,ik)
        if( ik.eq.1 ) then
          npop = npop1k
          efermi = ehomo1k
          ehomo  = ehomo1k
          elumo  = elumo1k
        else
c         ***** NB: redefine efermi,egap here
          npop = MAX( npop, npop1k )
          efermi = MIN( efermi, ehomo1k )
          ehomo  = MAX( ehomo, ehomo1k )
          elumo  = MIN( elumo, elumo1k )
        endif
c
        if( idebug.gt.2 .and. IWR.gt.0 )then
c         Dump out diagnostics
          write(IWR,9013) 'POPOCCDFT: nk,ik,wtk=',nk,ik,wtk(ik)
          write(IWR,9013)
     $      'POPOCCDFT: ik,npopk, esum1k, ehomo.elumo=',
     $      ik,npop1k, esum1k, ehomo,elumo, (elumo-ehomo)
        endif
c
  200 continue
      egap = elumo - ehomo
c
  900 continue
c
      if( idebug.gt.1 .and. IWR.gt.0 )then
c       Dump out diagnostics
 9013   format(1x,a,2i6,4f16.6)
        write(IWR,9013)
     $   'POPOCCDFT: sum nk,npop, esum1e, efermi,egap=',
     $    nk,npop, esum1e, efermi,egap
      endif
c
c    That's all Folks!
c
      RETURN
      END
