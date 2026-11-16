c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COFRF
c
c
      subroutine COFRF( c,gradc, r )
c---------------------------------------------------------------
c Purpose: calculate r-dependent part of coefficient matrix
c          for evaluation of electron-nucleus force elements
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c assume c() initialized elsewhere to include all r-dependence
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  r(3), c(10,9),gradc(3,10,9)
      DATA  two,four / 2.d0,4.d0 /
c
      x = r(1)
      y = r(2)
      z = r(3)
c
      c(1,2) = -x
      c(1,3) = -y
      c(1,4) = -z
      c(1,5) = x*y
      c(1,6) = y*z
      c(1,7) = z*x
      c(1,8) = x**2-y**2
      c(1,9) = two*z**2-x**2-y**2
c
      c(2,5) = -y
      c(3,5) = -x
      c(3,6) = -z
      c(4,6) = -y
      c(2,7) = -z
      c(4,7) = -x
      c(2,8) = -two*x
      c(3,8) =  two*y
      c(2,9) =  two*x
      c(3,9) =  two*y
      c(4,9) = -four*z
c
      gradc(1,1,5) =  y
      gradc(2,1,5) =  x
      gradc(2,1,6) =  z
      gradc(3,1,6) =  y
      gradc(1,1,7) =  z
      gradc(3,1,7) =  x
      gradc(1,1,8) =  two*x
      gradc(2,1,8) = -two*y
      gradc(1,1,9) = -two*x
      gradc(2,1,9) = -two*y
      gradc(3,1,9) =  four*z
c
      RETURN
      END
