c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POP
c
c
      subroutine POP( nspin, spindif,
     $ fermilvl, egap,
     $ etemp, elnumbr, esumsp, occfull,
     $ nstate,nk, wtk, eigval,eigpop,numpop,
     $ elvl,wlvl,klvl,nlvl,flvl )
c---------------------------------------------------------------
c Purpose: compute common (across spin) fermi level, up/dn-electrons
c
c Written: Peter A. Schultz, 12-January-2008, for 2.61
c
c Revision history:
c  30Sep08-PAS/2.62i: fix spin/non-spin initialization
c  12Jan08-PAS:2.61: adapted from function of original pop
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.d0, two=2.d0 )
c
c Input:
      DIMENSION  wtk(nk), eigval(*)
c Output:
      DIMENSION  eigpop(nstate,*)
c Scratch:
      DIMENSION  elvl(nstate*nk*nspin),klvl(*),nlvl(*)
      DIMENSION  wlvl(*),flvl(*)
c
c Local:
      DIMENSION  elspin(2), bandsum(2), npop(2)
      LOGICAL    xch,pos,neg
      DATA  elcerr,elconv / 1.d-6,1.d-11 /
      DATA  occnil / 1.d-9 /
c
c >>>> EXECUTABLE CODE:
c
      nstk = nstate*nk
      nstall = nstk*nspin
c
c  Sort the eigenstates by eigenvalue across all k-points
c  Do this by treating dn-spin (if I have two spins) as
c  second set of nk points
c
c First, load k-pt weights for sort, doubling nk for spin
c Use flvl as scratch for this.
      iksp = 0
      do  ispin=1,nspin
        do  ik=1,nk
          iksp = iksp + 1
          flvl(iksp) = wtk(ik)
        enddo
      enddo
      nksp = iksp
c
c Sort eigenstates across all k-points, and spin
      call LVLSORT( nstate,nksp, eigval,flvl,
     $ elvl,klvl,nlvl,wlvl )
c
c Initialize occupation factors to zero
      call MKZERO( nstall, flvl )
      call MKZERO( nstall, eigpop )
c
c >>>>> Find the fermi level
c
c Find the highest occupied (0K) orbital
        elnumt = occnil
        do 140 ist=1,nstall
        elnumt = elnumt + wlvl(ist)*occfull
c If no weight, skip level:
        if( wlvl(ist) .gt. zero )then
          if( elnumt .lt. elnumbr ) goto 140
          ihomo = ist
          ehomo = elvl(ihomo)
          goto 141
        endif
  140 continue
 1323 format(1x,a,3f12.6)
      call STOPXERR( 'elecnum - POP electron number' )
  141 continue
c
c Find first safely empty orbital
      ilumo = ihomo + 1
      elumo = ehomo + 40.d0*etemp
      do 150 ist=ihomo+1,nstall
        if( wlvl(ist) .gt. zero )then
          ilumo = ist
          elumo = elvl(ilumo)
          goto 151
        endif
  150 continue
  151 continue
c
c Compute first estimate of fermi level, and 0K gap
      efguess = ( ehomo + elumo ) / two
      gap = elumo - ehomo
      if( (elnumt-occnil) .gt. (elnumbr+occnil) )then
c       partially occupied 0K orbital -> gap=0
        gap = zero
      endif
      egap = gap
c
c Find average level spacing around fermi level
      de = ( elumo - ehomo ) / 4.d0
      if( de .lt. 10.d0*etemp ) de = 10.d0*etemp
c
c Find correct fermi level by interval halving
c
      itmax = 50
      pos = .false.
      neg = .false.
c
      do 200 itfermi=1,itmax
c
c       Compute electron counts for this temp/fermi level
        call FERMIPOP( efguess,etemp, elnumt, nstall,nstpop,
     $   occfull, flvl,elvl,wlvl )
c
c       Check if we have converged a fermi level (or close enough)
c        Look at total electron count ...
        elmiss = elnumt - elnumbr
        if( ABS( elmiss ) .lt. elconv ) goto 300
        if( itfermi .eq. itmax .and.
     $      ABS( elmiss ) .lt. elcerr ) goto 300
c
c       Adjust fermi level by interval halving
        if( elmiss .lt. zero )then
          neg = .true.
          if( pos .and. neg ) de = de / two
          efguess = efguess + de
        else
          pos = .true.
          if( pos .and. neg ) de = de / two
          efguess = efguess - de
        endif
c
  200 continue
      call STOPXERR( 'SPINfermi/ efermi not found' )
c
c >>>>> Fermi energy converged
c
  300 continue
      fermilvl = efguess
c
      do  ispin=1,2
        npop(ispin) = 0
        bandsum(ispin) = zero
        elspin(ispin) = zero
      enddo
c
c Count up number of electrons in each spin, with this fermi level
      do 310 ist=1,nstpop
        n = nlvl(ist)
        k = klvl(ist)
        if( flvl(ist) .eq. zero ) goto 311
c
        eigpop(n,k) = flvl(ist)
        elst = wlvl(ist)*flvl(ist)
        enst = elst*elvl(ist)
        if( k .gt. nk )then
c         dn-spin electrons
          elspin(2) = elspin(2) + elst
          bandsum(2) = bandsum(2) + enst
          if( nlvl(ist) .gt. npop(2) ) npop(2) = nlvl(ist)
        else
c         up-spin electrons
          elspin(1) = elspin(1) + elst
          bandsum(1) = bandsum(1) + enst
          if( nlvl(ist) .gt. npop(1) ) npop(1) = nlvl(ist)
        endif
  310 continue
  311 continue
c
      numpop = MAX( npop(1), npop(2) )
      esumsp = bandsum(1) + bandsum(2)
      spindif = elspin(1) - elspin(2)
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVLSORT
c
c
      subroutine LVLSORT( nstate,nk,
     $ eigval, wtk,
     $ elvl,klvl,nlvl,wlvl )
c---------------------------------------------------------------
c Purpose: sort of full eigenspectrum over multiple k-points
c
c Written: Peter A. Schultz, 10-January-2008
c
c Revision history:
c  none
c---------------------------------------------------------------
c Notes on this eigenstate sorting scheme:
c   I know that states within each k-point are already sorted
c   I will use this knowledge to optimize this sort
c   (1a) Initialize global sort with the lowest state from each k-point
c   (1b) Sort those nk states within that buffer by energy, tracking n,k
c   (2) Accept the lowest state in that active buffer into next position
c       in global sorted list, bumping active sort buffer up one.
c   (3) Take the next-highest state (if there are any more) in
c       the k-set corresponding to the state just accepted,
c       insert it into the correct location in the active buffer.
c       Start from the top in this small sort, as the new state is
c       likely to be close to the end.
c   (4) Return to (2), until all states have been ranked
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  eigval(nstate,nk)
      DIMENSION  wtk(nk)
c Output:
      DIMENSION  elvl(nstate*nk),klvl(*),nlvl(*),wlvl(*)
c
c Local:
      LOGICAL    xch
c
c >>>> EXECUTABLE CODE:
c
      nlevels = nstate*nk
c
c Load bottom state in each k into sort buffer
      lvls = 0
      do  ik=1,nk
        lvls = lvls + 1
        elvl(lvls) = eigval(1,ik)
        klvl(lvls) = ik
        nlvl(lvls) = 1
      enddo
c
c Sort these with a simple pair exchange (k is never huge)
  100 continue
      xch = .false.
      do  jlvl=1,lvls-1
        if( elvl(jlvl).gt.elvl(jlvl+1) )then
c         Move this state upwards
          ekj = elvl(jlvl)
          elvl(jlvl) = elvl(jlvl+1)
          elvl(jlvl+1) = ekj
          kkj = klvl(jlvl)
          klvl(jlvl) = klvl(jlvl+1)
          klvl(jlvl+1) = kkj
c         These nlvl() are all the first state, exchange is redundant
          xch = .true.
        endif
      enddo
      if( xch ) goto 100
c
c Now sort through the remainder of the states
c
      ilvl = 1
c
  200 continue
c
c The state ilvl is lowest in active sort buffer (lowest-energy state
c from each k-point), and thereby guaranteed to be next lowest in
c global sort:
      kprev = klvl(ilvl)
      nprev = nlvl(ilvl)
      ilvl = ilvl + 1
c Insert next higher state from that k-pt into active set (if exists):
      if( nprev .lt. nstate )then
        lvls = lvls + 1
c       Rank eigval from next state in this k-pt in active sort buffer
        ekn = eigval(nprev+1,kprev)
        lvlnew = ilvl
        do 210 jlvl=lvls,ilvl+1,-1
          if( ekn .lt. elvl(jlvl-1) )then
            elvl(jlvl) = elvl(jlvl-1)
            klvl(jlvl) = klvl(jlvl-1)
            nlvl(jlvl) = nlvl(jlvl-1)
          else
            lvlnew = jlvl
            goto 211
          endif
  210   continue
  211   continue
        elvl(lvlnew) = ekn
        klvl(lvlnew) = kprev
        nlvl(lvlnew) = nprev + 1
      endif
c
c If not all levels ranked, go sort some more
      if( lvls .lt. nlevels ) goto 200
c
c Done sorting, associate states with correct k-sample weighting factors:
      do  ilvl=1,nlevels
        k = klvl(ilvl)
        wlvl(ilvl) = wtk(k)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
