c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> YLMREAL
c
c
      subroutine YLMREAL( ylm,lmax,r )
c---------------------------------------------------------------
c Purpose: compute real spherical harmonics for m=0,l=0,1,...lmax
c          m=1,2,...,lmax;m=lmax,l=lmax
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c  The real ylm's are defined by:
c      ylm-plus  = ylm, if m=0
c      ylm-minus = 0, if m=0
c      ylm-plus  = sqrt(2.)*real(ylm), if m.ne.0
c      ylm-minus = sqrt(2.)*imag(ylm), if m.ne.0
c
c      the ylm-pluses are returned in the real part of array "ylm"
c      the ylm-minuses are returned in the imaginary part
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DOUBLE COMPLEX  ylm(*),  ylp2m
c
      DIMENSION  r(3)
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
      rt2 = SQRT(two)
      rpsq = r(1)**2 + r(2)**2
      rpar = SQRT(rpsq)
      rmag = SQRT(r(3)**2+rpsq)
      if( rmag.eq.zero ) call STOPXERR( 'ylmrlr=0' )
      costh = r(3)/rmag
      sinth = rpar/rmag
      if( rpar.eq.zero )then
        sinph = zero
        cosph = one
      else
        cosph = r(1)/rpar
        sinph = r(2)/rpar
      endif
c
c      compute  ylm's:
c
      lmaxp1 = lmax + 1
      do 20 mm=1,lmaxp1
        m = mm - 1
        mstore = 1 + (lmax+2)*m - (m*(m+1))/2
        call YMM( ylm(mstore),sinth,cosph,sinph,m )
        if( m.ne.0 ) ylm(mstore) = ylm(mstore)*rt2
        if( m.eq.lmax ) goto 20
        mstore = mstore + 1
        ylm(mstore) = SQRT( DBLE(2*m+3))*costh*ylm(mstore-1)
        if( m.eq.lmax-1 ) goto 20
c generate remaining ylm's by recursion
        mp2 = m + 2
        do 10 l=mp2,lmax
          mstore = mstore + 1
          ylm(mstore) = YLP2M(ylm(mstore-2),ylm(mstore-1),costh,l-2,m)
   10   continue
   20 continue
c
      RETURN
      END
