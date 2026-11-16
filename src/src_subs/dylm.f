c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DYLM
c
c
      subroutine DYLM( ylm,gradylm,lmax,r )
c---------------------------------------------------------------
c Purpose: Calculate real spherical harmonics for
c      m=0,l=0,1,...,lmax;m=1,l=1,2,...,lmax;...;m=lmax,l=lmax
c
c Revision history:
c  11Mar02-PAS/2.52: last literal reals to d0
c  21Jun01-PAS/2.48: replace STOPs
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  2*(lmax+1)**2 real ylm's are stored in the complex array "ylm."
c  The real ylm's are defined by:
c                   m=0         m.ne.0
c      ylm-plus  =  ylm    sqrt(2.)*real(ylm)
c      ylm-minus =   0     sqrt(2.)*imag(ylm)
c  ylm-pluses returned real part of "ylm", ylm-minus in imaginary
c
c  Also calculate:   rmag**(-l) * grad (rmag**l * ylmreal)
c  These results are returned in the array "gradylm"
c  To generate gradylm's, make use of the inhomogeneous recursion
c  formula they satisfy, which can be applied once the ylm's known
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DOUBLE COMPLEX    ylm(*),gradylm(3,*)
c local declarations:
      DOUBLE COMPLEX    recur
      DIMENSION  r(3)
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      sq2 = SQRT( two )
      rpsq = r(1)**2 + r(2)**2
      rpar = SQRT( rpsq )
      rsq = r(3)**2 + rpsq
      rmag = SQRT( rsq )
      if( rmag.eq.zero ) call STOPXERR( 'dylm-r=0' )
c
c    Compute theta and phi
c
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
c    Compute ylm's
c
      lmaxp1 = lmax + 1
      rtom = one
      do 100 mm=1,lmaxp1
        m = mm - 1
        mstore = 1 + (lmax+2)*m - (m*(m+1))/2
        mstmm1 = 1 + (lmax+2)*(m-1) - (m*(m-1))/2
        call YMM( ylm(mstore),sinth,cosph,sinph,m )
        gradylm(3,mstore) = zero
        if( m.eq.0 )then
          gradylm(1,1) = zero
          gradylm(2,1) = zero
        else
          ylm(mstore) = ylm(mstore)*sq2
          recur = -SQRT( DBLE((2*m+1)*2*m) )*ylm(mstmm1)*rtom/two
          rtom = rtom*rmag
          if( mstmm1.eq.1 ) recur = recur*sq2
          gradylm(1,mstore) = recur
          gradylm(2,mstore) = recur*(0.d0,1.d0)
        endif
        if( m.eq.lmax ) goto 110
c
        mstore = mstore + 1
        sq2mp3 = SQRT( DBLE(2*m+3) )
        ylm(mstore) = sq2mp3*costh*ylm(mstore-1)
        recur = sq2mp3*r(3)
        do 20 n=1,3
          gradylm(n,mstore)=recur*gradylm(n,mstore-1)
   20   continue
        gradylm(3,mstore) = gradylm(3,mstore) + sq2mp3*rtom*
     $   ylm(mstore-1)
        if(m.eq.lmax-1) goto 100
c
c    Generate remaining ylm's by recursion
c
        mp2 = m + 2
        rtol = rmag**m
        do 90 l=mp2,lmax
          mstore = mstore + 1
          call YRECUR( clp1m,clm,l-2,m )
          ylm(mstore) = costh*clp1m*ylm(mstore-1) - clm*ylm(mstore-2)
          do 30 n=1,3
            gradylm(n,mstore) = clp1m*r(3)*gradylm(n,mstore-1) -
     $       clm*(rsq*gradylm(n,mstore-2)+two*r(n)*rtol*ylm(mstore-2))
   30     continue
          gradylm(3,mstore) = gradylm(3,mstore) +
     $     clp1m*rtol*rmag*ylm(mstore-1)
          rtol = rtol*rmag
   90   continue
c
  100 continue
c
c    Scale gradylm's so the output includes rmag**(-l)
c
  110 mstore = 0
      do 200 mm=1,lmaxp1
        m = mm - 1
        rtol = rmag**m
        do 190 ll=mm,lmaxp1
          mstore = mstore + 1
          do 180 n=1,3
            gradylm(n,mstore) = gradylm(n,mstore)/rtol
  180     continue
          rtol = rtol*rmag
  190   continue
  200 continue
c
c    That's all Folks!
c
      RETURN
      END
