c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPLIND2
c
c
      subroutine SPLIND2( na, xa,ya, y2,  nb, xb, ddyb )
c---------------------------------------------------------------
c Purpose: compute 2nd derivative of cubic spline
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  xa(na),ya(na), y2(na)
      DIMENSION  xb(*),ddyb(*)
c
      xmin = xa(1)
      xmax = xa(na)
      do 100 i=1,nb
        x = xb(i)
c
c out of range test
        if( x .gt. xmax )then
          ddyb(i) = 0.d0
          goto 100
        elseif( x .lt. xmin )then
          ddyb(i) = 0.d0
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
        ddyb(i) = a*y2(klo) + b*y2(khi)
  100 continue
c
      RETURN
      END
