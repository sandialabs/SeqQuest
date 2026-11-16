c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPLINE
c
c
      subroutine SPLINE( n, x, y, ddy0, ddy1, y2,u )
c---------------------------------------------------------------
c Purpose: compute cubic spline fit
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  x(n),y(n), y2(n)
      DIMENSION  u(n)
c
c compute the data for the spline
      if( ddy0 .gt. .99d30 ) then
        y2(1) = 0.d0
        u(1) = 0.d0
      else
        y2(1) = -.5d0
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        u(1) = 3.d0*(dy/dx - ddy0)/dx
      endif
c
      do 11 i=2,n-1
        dx = x(i) - x(i-1)
        dx1 = x(i+1) - x(i)
        d2x = x(i+1) - x(i-1)
        sig = dx/d2x
        p = sig*y2(i-1) + 2.d0
        y2(i) = (sig - 1.d0)/p
        dy = y(i) - y(i-1)
        dy1 = y(i+1) - y(i)
        u(i) = ( 6.d0*(dy1/dx1-dy/dx)/d2x - sig*u(i-1) )/p
   11 continue
c
      if( ddy1 .gt. .99d30 )then
        qn = 0.d0
        un = 0.d0
      else
        qn = .5d0
        dx = x(n) - x(n-1)
        dy = y(n) - y(n-1)
        un = (ddy1 - dy/dx)*3.d0/dx
      endif
c
      y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
c
c back substitute
      do 12 i=n-1,1,-1
        y2(i) = y2(i)*y2(i+1) + u(i)
   12 continue
c
      RETURN
      END
