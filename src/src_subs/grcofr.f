c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRCOFR
c
c
      subroutine GRCOFR( c, r )
c---------------------------------------------------------------
c Purpose: calculates c() matrix for 2nd derivatives in v1ctr
c
c Revision history:
c  19Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION r(3), c(10,3,3)
      DATA  zero,one,two,four / 0.d0,1.d0,2.d0,4.d0 /
c
c >>>> EXECUTABLE CODE:
c
      x = r(1)
      y = r(2)
      z = r(3)
c
      do  l=1,3
        do  ll=1,3
          do  n=1,10
            c(n,l,ll) = zero
          enddo
        enddo
      enddo
c
      c(5,1,2) = one
      c(6,2,3) = one
      c(7,1,3) = one
      c(8,1,1) = one
      c(9,2,2) = one
      c(10,3,3)= one
      c(1,1,2) = x*y
      c(1,2,3) = y*z
      c(1,1,3) = x*z
      c(1,1,1) = x**2
      c(1,2,2) = y**2
      c(1,3,3) = z**2
      c(2,1,2) = -y
      c(3,1,2) = -x
      c(3,2,3) = -z
      c(4,2,3) = -y
      c(2,1,3) = -z
      c(4,1,3) = -x
      c(2,1,1) = -two*x
      c(3,2,2) = -two*y
      c(4,3,3) = -two*z
c
c    That's all Folks!
c
      RETURN
      END
