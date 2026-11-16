c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GNORM
c
c
      subroutine GNORM( l, n, a, c )
c---------------------------------------------------------------
c Purpose: normalize contracted gaussians
c
c Written: Peter A. Schultz, 27-July-1998, for v2.23 (split basis)
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(n),c(n)
c local declarations:
      DATA  one / 1.d0 /
c
      call GOVLP( s, l, n,a,c, n,a,c )
c
      snorm = one/SQRT(s)
      do  i=1,n
        c(i) = snorm*c(i)
      enddo
c
      RETURN
      END
