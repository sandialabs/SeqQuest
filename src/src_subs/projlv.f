c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PROJLV
c
c
      subroutine PROJLV( r, a,b,c, v )
c---------------------------------------------------------------
c Purpose: project vector r() onto linearly independent vectors
c          a(), b(), c(), putting result (Cramer's rule) in v().
c
c Written: Peter A. Schultz, 8-May-1991
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  r(3),a(3),b(3),c(3)
c output arrays:
      DIMENSION  v(3)
c
c >>>> EXECUTABLE CODE:
c
      denom =  a(1)*(b(2)*c(3)-b(3)*c(2)) +
     $         b(1)*(c(2)*a(3)-c(3)*a(2)) +
     $         c(1)*(a(2)*b(3)-a(3)*b(2))
c
      v(1) = ( r(1)*(b(2)*c(3)-b(3)*c(2)) +
     $         b(1)*(c(2)*r(3)-c(3)*r(2)) +
     $         c(1)*(r(2)*b(3)-r(3)*b(2)) ) / denom
      v(2) = ( a(1)*(r(2)*c(3)-r(3)*c(2)) +
     $         r(1)*(c(2)*a(3)-c(3)*a(2)) +
     $         c(1)*(a(2)*r(3)-a(3)*r(2)) ) / denom
      v(3) = ( a(1)*(b(2)*r(3)-b(3)*r(2)) +
     $         b(1)*(r(2)*a(3)-r(3)*a(2)) +
     $         r(1)*(a(2)*b(3)-a(3)*b(2)) ) / denom
c
      RETURN
      END
