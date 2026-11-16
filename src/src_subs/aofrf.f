c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AOFRF
c
c
      subroutine AOFRF( a,grada, r )
c---------------------------------------------------------------
c Purpose: calculate r-dependent part of coefficient matrix
c          for evaluation of non-local force matrix elements
c
c Revision history:
c   7Mar02-PAS/2.52: explicitly promote constants to dp
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c assume a() initialized elsewhere to include all r-dependence
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  r(3), a(9,9),grada(3,9,9)
c
      DATA  rt4pi   / 3.544907701811d0 /
      DATA  rt4pid3 / 2.046653415893d0 /
      DATA  two,four / 2.d0,4.d0 /
c
c >>>> EXECUTABLE CODE:
c
      x = r(1)
      y = r(2)
      z = r(3)
c
      a(1,2) = -x * rt4pi
      a(1,3) = -y * rt4pi
      a(1,4) = -z * rt4pi
      a(1,5) = x*y * rt4pi
      a(1,6) = y*z * rt4pi
      a(1,7) = z*x * rt4pi
      a(1,8) = (x**2-y**2) * rt4pi
      a(1,9) = (two*z**2-x**2-y**2) * rt4pi
c
      a(2,5) =  y * rt4pid3
      a(3,5) =  x * rt4pid3
      a(3,6) =  z * rt4pid3
      a(4,6) = -y * rt4pid3
      a(2,7) =  z * rt4pid3
      a(4,7) = -x * rt4pid3
      a(2,8) =  two*x * rt4pid3
      a(3,8) = -two*y * rt4pid3
      a(2,9) = -two*x * rt4pid3
      a(3,9) = -two*y * rt4pid3
      a(4,9) = -four*z * rt4pid3
c
      grada(1,1,5) =  y * rt4pi
      grada(2,1,5) =  x * rt4pi
      grada(2,1,6) =  z * rt4pi
      grada(3,1,6) =  y * rt4pi
      grada(1,1,7) =  z * rt4pi
      grada(3,1,7) =  x * rt4pi
      grada(1,1,8) =  two*x * rt4pi
      grada(2,1,8) = -two*y * rt4pi
      grada(1,1,9) = -two*x * rt4pi
      grada(2,1,9) = -two*y * rt4pi
      grada(3,1,9) =  four*z * rt4pi
c
c    That's all Folks!
c
      RETURN
      END
