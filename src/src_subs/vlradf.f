c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLRADF
c
c
      subroutine VLRADF( vint1,vint2,vint3, calc1,calc2,calc3,
     $ lmx,alsum,pmag,exparg,
     $ nrs, radmsh,rwght,  pot1,pot2,pot3,
     $ besstor, w1,w2,w3,w4,w5,w6 )
c---------------------------------------------------------------
c Purpose: integrals for n=0,ltop-1 off radial potentials of:
c       exp(-alsum*r**2) * r**(2+n) * i-sub-n(pmag*r) / pmag**n
c      "i-sub-n" is the spherical bessel function, imag argument
c
c Revision history:
c   6Mar00-PAS/2.43: rename "shtrngf" to "vlradf"; simple opt.
c  19Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c   r**2 factor needed for radial integral included in rwght
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vint1(*),vint2(*),vint3(*)
c input arrays:
      LOGICAL    calc1,calc2,calc3
      DIMENSION  radmsh(nrs),rwght(nrs)
      DIMENSION  pot1(nrs),pot2(nrs),pot3(nrs)
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
          w1(nr) = rwght(nr)*EXP(-alsum*r*r-exparg)
          besstor(nr,1) = one
        enddo
c
        do  ll=2,ltop
          denom = one/ DBLE(2*ll-1)
          do  nr=1,nrs
            besstor(nr,ll) = besstor(nr,ll-1)*radmsh(nr)*denom
          enddo
        enddo
c
        factor = one
c
      endif
c
c ***** bessel functions done, evaluate radial integrals *****
c
      do 200 ll=1,ltop
c
        do  nr=1,nrs
          w2(nr) = besstor(nr,ll)*w1(nr)
        enddo
c
        if( calc1 .and. calc2 )then
c         This is most common combination(es and pspot):
          do  nr=1,nrs
            vint1(ll) = vint1(ll) + w2(nr)*pot1(nr)
            vint2(ll) = vint2(ll) + w2(nr)*pot2(nr)
          enddo
        elseif( calc1 )then
c         This is second most common path (es only):
          do  nr=1,nrs
            vint1(ll) = vint1(ll) + w2(nr)*pot1(nr)
          enddo
        elseif( calc2 )then
          do  nr=1,nrs
            vint2(ll) = vint2(ll) + w2(nr)*pot2(nr)
          enddo
        endif
c
        if( calc3 )then
c         This is xc-local, usually inactive:
          do  nr=1,nrs
            vint3(ll) = vint3(ll) + w2(nr)*pot3(nr)
          enddo
        endif
c
        if( ll.lt.ltop )then
          do  nr=1,nrs
            w1(nr) = w1(nr)*radmsh(nr)*factor
          enddo
        endif
c
  200 continue
c
c    That's all Folks!
c
      RETURN
      END
