c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVPROJ
c
c
      subroutine LVPROJ( a, r, v )
c---------------------------------------------------------------
c Purpose: project vector r() onto primitive vectors; and v.v.
c
c Written: Peter A. Schultz, 7-Jul-2000, for v2.45
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Uses Cramer's Rule to do projections.
c   a() are the input primitive vectors
c   r() is a vector in Cartesian space
c   v() is projection onto the primitive vectors
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  a(3,3), r(3)
c output arrays:
      DIMENSION  v(3)
c
c >>>> EXECUTABLE CODE:
c
      denom =  a(1,1)*( a(2,2)*a(3,3) - a(3,2)*a(2,3) ) +
     $         a(1,2)*( a(2,3)*a(3,1) - a(3,3)*a(2,1) ) +
     $         a(1,3)*( a(2,1)*a(3,2) - a(3,1)*a(2,2) )
c
      v(1) = ( r(1)  *( a(2,2)*a(3,3) - a(3,2)*a(2,3) ) +
     $         a(1,2)*( a(2,3)*r(3)   - a(3,3)*r(2)   ) +
     $         a(1,3)*( r(2)  *a(3,2) - r(3)  *a(2,2) ) ) / denom
      v(2) = ( a(1,1)*( r(2)  *a(3,3) - r(3)  *a(2,3) ) +
     $         r(1)  *( a(2,3)*a(3,1) - a(3,3)*a(2,1) ) +
     $         a(1,3)*( a(2,1)*r(3)   - a(3,1)*r(2)   ) ) / denom
      v(3) = ( a(1,1)*( a(2,2)*r(3)   - a(3,2)*r(2)   ) +
     $         a(1,2)*( r(2)  *a(3,1) - r(3)  *a(2,1) ) +
     $         r(1)  *( a(2,1)*a(3,2) - a(3,1)*a(2,2) ) ) / denom
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVXPND
c
c
      entry      LVXPND( a, r, v )
c
c >>>> EXECUTABLE CODE:
c
      r(1) = v(1)*a(1,1) + v(2)*a(1,2) + v(3)*a(1,3)
      r(2) = v(1)*a(2,1) + v(2)*a(2,2) + v(3)*a(2,3)
      r(3) = v(1)*a(3,1) + v(2)*a(3,2) + v(3)*a(3,3)
c
      RETURN
c
c    That's all Folks!
c
      END
