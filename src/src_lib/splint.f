c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPLINT
c
c
      subroutine SPLINT( na, xa,ya,y2,  nb, xb,yb )
c---------------------------------------------------------------
c Purpose: compute function defined by cubic spline
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  xa(na),ya(na), y2(na)
      DIMENSION  xb(*),yb(*)
c
      xmin = xa(1)
      xmax = xa(na)
      do 100 i=1,nb
        x = xb(i)
c
c out of range test
        if( x .ge. xmax )then
          yb(i) = ya(na)
          goto 100
        elseif( x .le. xmin )then
          yb(i) = ya(1)
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
        c = (a**3 - a)*h2
        d = (b**3 - b)*h2
        yb(i) = a*ya(klo) + b*ya(khi)
     $        + c*y2(klo) + d*y2(khi)
100   continue
c
      RETURN
      END
