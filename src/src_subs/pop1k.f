c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POP1K
c
c
      subroutine POP1K( closeshl,
     $ fermilvl,gap, ehomo,elumo,
     $ etemp,edegen, elnumbr,occfull,
     $ nstate,numpop, wtocc,
     $ eigval,eigpop, wlvl )
c---------------------------------------------------------------
c Purpose: construct level populations for eigenstates,
c          with option to impose 0K (closed shell) occupancies,
c          with proper treatment of degeneracies.
c
c Written: Peter A. Schultz, 20-February-2005, for v2.59
c          (extracted from real solver, and expanded)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0 )
c
c Input:
      LOGICAL    closeshl, doblas3
      DIMENSION  eigval(nstate)
c Output:
      DIMENSION  eigpop(nstate)
c Scratch arrays:
      DIMENSION  wlvl(nstate)
c
c Local:
      LOGICAL    pos,neg,closed
c  elconv=difference criterion for convergence of electron count
c  elcerr=worst tolerated error in electron count
      DATA  elconv,elcerr / 1.d-11,1.d-6 /
c  occnil=small number to make arithmetic safe for level occ manipulation
      DATA  occnil / 1.d-9 /
c  itmax=maximum number of iterations for fermi level finding:
      DATA  itmax / 50 /
c
c >>>> EXECUTABLE CODE:
c
c Clear out occupancy levels:
      numpop = 0
      call MKZERO( nstate, eigpop )
      if( elnumbr .le. zero )then
        elumo = eigval(1)
        ehomo = elumo
        fermilvl = ehomo
        gap = -one
        numpop = 0
        goto 900
      endif
c
c  Compute number of 0K occupied levels:
      nocc0k = ( elnumbr + occfull - occnil ) / occfull
      elfull = occfull*DBLE( nocc0k )
      elpart = ABS( elfull - elnumbr )
      nstfull = nocc0k
      if( elpart .gt. occnil ) nstfull = nocc0k - 1
c
      ihomo = nocc0k
      ehomo = eigval(ihomo)
c  Default fermi level to highest occupied state:
      elumo = ehomo
      fermilvl = ehomo
c
c Number of degenerate states at top, default to non-degen:
      ndegen = 1
c
      if( nstate .eq. nocc0k )then
c       Have only the occupied states, no unoccupied states
        closed = .true.
        if( elpart .gt. occnil )then
          gap = zero
        else
          gap = -one
        endif
c
      elseif( nstate .gt. nocc0k )then
c       Have empty states - can compute a gap
        closed = .false.
        ilumo = ihomo + 1
        elumo = eigval(ilumo)
        gap = elumo - ehomo
        fermilvl = (ehomo + elumo) / two
        if( closeshl )then
c         Impose closed shell (0K) occups, while checking for degens
          closed = .true.
c         I will allow for 5-fold degeneracy at most
c         Check 4 levels above homo for degeneracy with homo:
          itop = ihomo + 4
          itop = MIN( itop, nstate )
          do  ilvl=ihomo+1,itop
            eigi = eigval(ilvl)
            if( ABS( eigi-ehomo ) .lt. edegen )then
c             Add to degeneracy, and to 0K occup list
              ndegen = ndegen + 1
              if( nstfull.eq.ihomo ) nstfull = nstfull - 1
              nocc0k = nocc0k + 1
            endif
          enddo
c         If we have a degeneracy or partial occupancy, need
c         to check below to see if degeneracy extends deeper:
          if( ndegen .gt. 1 .or. elpart.gt.occnil )then
            nhomo = ihomo - (5-ndegen)
            do  ilvl=ihomo-1,nhomo,-1
              eigi = eigval(ilvl)
              if( ABS( eigi-ehomo ) .lt. edegen )then
c               Add to degeneracy, remove from full list
                ndegen = ndegen + 1
                nstfull = nstfull - 1
              endif
            enddo
          endif
        else
c         Determine if we have a closed shell system
          if( elpart .lt. occnil .and. gap .gt. 19*etemp )then
c           Closed shell anyway
            closed = .true.
          else
c           Will have partially occupied orbitals
            closed = .false.
            if( elpart .ge. occnil )then
              gap = zero
              fermilvl = ehomo
            endif
          endif
        endif
c
      endif
c
c Do 0K occupancies if we have closed shell system
c
      if( closed )then
        numpop = nocc0k
        do  ilvl=1,nstfull
          eigpop(ilvl) = occfull
        enddo
        if( ndegen .gt. 1 )then
          eldegen = elnumbr - occfull*DBLE( nstfull )
          occdegen = eldegen / DBLE( ndegen )
          do  ilvl=nstfull+1,nstfull+ndegen
            eigpop(ilvl) = occdegen
          enddo
        elseif( elpart .gt. occnil )then
c         Top level is only partially occupied
          eigpop(nocc0k) = elpart
        endif
c       Skip the fermi level determination
        goto 900
      endif
c
c Find average level spacing around fermi level
      de = ( elumo - ehomo ) / 4.d0
      if( de .lt.  5.d0*etemp ) de =  5.d0*etemp
c
c Find correct fermi level by interval halving
c
c     Set up weighting functions for electron counting:
      do  ilvl=1,nstate
        wlvl(ilvl) = one
      enddo
c
      pos = .false.
      neg = .false.
c
      do 200 itfermi=1,itmax
c       Compute electron counts for this temp and fermi level
        call FERMIPOP( fermilvl,etemp, elsum, nstate,numpop,
     $   occfull, eigpop, eigval, wlvl )
c
c       Check if we have good electron total (or close enough)
        elnumt = elsum*wtocc
        elmiss = elnumt - elnumbr
        if( ABS( elmiss ) .lt. elconv ) goto 900
        if( itfermi .eq. itmax .and.
     $      ABS( elmiss ) .lt. elcerr ) goto 900
c
c       Adjust fermi level by interval halving
        if( elmiss .lt. zero )then
          neg = .true.
          if( pos .and. neg ) de = de / two
          fermilvl = fermilvl + de
        else
          pos = .true.
          if( pos .and. neg ) de = de / two
          fermilvl = fermilvl - de
        endif
c       Go try this new fermi level
  200 continue
      call STOPXERR( 'efermirl/erlsolv efermi conv' )
c
c >>>>> Fermi energy converged
c
  900 continue
c
      RETURN
c
c    That's all Folks!
c
      END
