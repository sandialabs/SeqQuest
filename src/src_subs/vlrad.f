c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLRAD
c
c
      subroutine VLRAD( vint1,vint2,vint3,vint4,
     $ lmx,alsum,pmag,exparg,
     $ nrs, radmsh,rwght,  pot1,pot2,pot3,pot4,
     $ besstor,w1,w2,w3,w4,w5,w6 )
c---------------------------------------------------------------
c Purpose: does necessary radial integrals for local potentials
c
c Revision history:
c   6Mar00-PAS?2.43: rename "shtrng" to "vlrad"
c   7Dec99-PAS/2.42: clean-up loops
c  10Dec98-PAS/2.29: promote source to explicit double precision
c    May89-PAS: install opted bessel routines
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vint1(*),vint2(*),vint3(*),vint4(*)
c input arrays:
      DIMENSION  radmsh(nrs),rwght(nrs)
      DIMENSION  pot1(nrs),pot2(nrs),pot3(nrs),pot4(nrs)
c scratch arrays:
      DIMENSION  besstor(nrs,*)
     $   ,w1(nrs),w2(nrs),w3(nrs),w4(nrs),w5(nrs),w6(nrs)
c local declarations:
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      ltop = lmx
c
c ***** compute spherical bessel functions *****
c
      if( pmag.gt.zero )then
c
        call BSIALL( nrs, ltop, pmag, besstor, radmsh,
     $   w1,w2,w3,w4,w5,w6 )
c
        do  nr=1,nrs
          r = radmsh(nr)
          w1(nr) = rwght(nr)*EXP(-alsum*r**2+r*pmag-exparg)
        enddo
c
        factor = one/pmag
c
      else
c
        do  nr=1,nrs
          r = radmsh(nr)
          w1(nr) = rwght(nr)*EXP(-alsum*r**2-exparg)
          besstor(nr,1) = one
        enddo
c
        do 140 ll=2,ltop
          denom = one/ DBLE(2*ll-1)
          do  nr=1,nrs
            besstor(nr,ll) = besstor(nr,ll-1)*radmsh(nr)*denom
          enddo
  140   continue
c
        factor = one
c
      endif
c
c ***** bessel functions done, evaluate radial integrals *****
c
      do 240 ll=1,ltop
        do  nr=1,nrs
          rmul = besstor(nr,ll)*w1(nr)
          vint1(ll) = vint1(ll) + rmul*pot1(nr)
          vint2(ll) = vint2(ll) + rmul*pot2(nr)
          vint3(ll) = vint3(ll) + rmul*pot3(nr)
          vint4(ll) = vint4(ll) + rmul*pot4(nr)
          w1(nr) = w1(nr)*radmsh(nr)*factor
        enddo
  240 continue
c
c    That's all Folks!
c
      RETURN
      END
