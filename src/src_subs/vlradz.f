c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLRADZ
c
c
      subroutine VLRADZ( venint,ltop,alsum,pmag,znuc,exparg,expfac )
c---------------------------------------------------------------
c Purpose: calculates radial integrals with znuc/r
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c   6Mar00-PAS/2.43: rename "radint" to "vlradz"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c expfac = EXP(-exparg) with underflow control
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c output array:
      DIMENSION  venint(*)
c local declarations:: ( ERF has machine/MACHINE dependence)
      DATA  tol / 1.d-14 /
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      rtpid2 = rtpi/two
      twoz = two*znuc
      ltopm1 = ltop - 1
c
      arg = pmag/SQRT(alsum)
      argsq = arg*arg
      erfarg = arg/two
      phiarg = erfarg*erfarg
c
      if( ltop.gt.1 ) goto 10
      if( arg.gt.zero )then
        venint(1) = EXP( phiarg-exparg )*DERF( erfarg )*rtpid2/arg
      else
        venint(1) = expfac/two
      endif
      goto 100
c
   10 if( phiarg.lt. DBLE( ltop )-1.5 ) goto 30
c
c      do forward recursion for phiarg gt lmax-0.5
c
      venint(1) = EXP( phiarg-exparg )*DERF( erfarg )*rtpid2/arg
      venint(2) = (half-one/argsq)*venint(1) + half*expfac/argsq
      if( ltop.gt.2 )then
        do  ll=3,ltop
          l = ll - 1
          venint(ll) = (half/phiarg)*
     $     ( venint(ll-2)* DBLE( l-1 ) / two
     $     + venint(ll-1)*(phiarg- DBLE( l )+half) )
        enddo
      endif
      goto 100
c
c ---- land here if phiarg is small
c
   30 continue
c
      if( pmag.eq.zero )then
c       >>>> phiarg equals zero
        fac = half
        do  l=1,ltop
          venint(l) = fac*EXP( -exparg )
          fac = fac* DBLE( l ) / DBLE( 2*l+1 )
        enddo
        goto 100
      endif
c
c      >>>> phiarg is small, but non-zero
c      taylor series for largest l-s and do backward recursion
c
      fac = one
      do  l=1,ltopm1
        fac = fac* DBLE( l ) / DBLE( 2*l+1 )
      enddo
      ta1 = fac*EXP(phiarg-exparg)/two
      ta2 = ta1* DBLE( 2*ltopm1+1 ) / DBLE( ltopm1 )
      ltopm2 = ltopm1 - 1
      t1 = zero
      t2 = zero
      do  k=1,100
        t1 = t1 + ta1
        t2 = t2 + ta2
        ta1 = -ta1* DBLE(2*k-1)*phiarg/ DBLE(k*(2*(ltopm1+k)+1))
        ta2 = -ta2* DBLE(2*k-1)*phiarg/ DBLE(k*(2*(ltopm2+k)+1))
        if( ABS(ta1/t1).le.tol .and. ABS(ta2/t2).lt.tol ) goto 70
      enddo
      call STOPXERR( 'radint' )
c
c      back recursion
c
   70 continue
      venint(ltopm1+1) = t1
      venint(ltopm2+1) = t2
      do  l=3,ltop
        ll = ltop - l + 1
        fll =  DBLE(ll)
        venint(ll) = (two*phiarg*venint(ll+2) -
     $                        (phiarg-fll-half)*venint(ll+1) )
     $                   / (fll/two)
      enddo
c
c >>>> scale results
c
  100 continue
      alscal = alsum
      do  l=1,ltop
        venint(l) = twoz*venint(l)/alscal
        alscal = alscal*alsum
      enddo
c
c    That's all Folks!
c
      RETURN
      END
