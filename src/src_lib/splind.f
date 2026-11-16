c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPLIND
c
c
      subroutine SPLIND( na, xa,ya, y2,  nb, xb, dyb )
c---------------------------------------------------------------
c Purpose: compute derivative of cubic spline
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  xa(na),ya(na), y2(na)
      DIMENSION  xb(*),dyb(*)
c
      xmin = xa(1)
      xmax = xa(na)
      do 100 i=1,nb
        x = xb(i)
c
c out of range test
        if( x .ge. xmax )then
          dyb(i) = 0.d0
          goto 100
        elseif( x .lt. xmin )then
          dyb(i) = 0.d0
          goto 100
        endif
c
c in range do binary search
        klo = 1
        khi = na
    1   continue
        if( khi-klo .gt. 1 )then
          k = (khi+klo)/2
          if( xa(k) .gt. x )then
            khi = k
          else
            klo = k
          endif
          goto 1
        endif
c
        dx = xa(khi) - xa(klo)
        a = (xa(khi) - x)/dx
        b = 1.d0 - a
        h2 = (dx**2)/6.d0
        c = (3.d0*a**2 - 1.d0)*h2
        d = (3.d0*b**2 - 1.d0)*h2
        dyb(i) = (  -ya(klo) +   ya(khi)
     $           - c*y2(klo) + d*y2(khi) ) / dx
  100 continue
c
      RETURN
      END
