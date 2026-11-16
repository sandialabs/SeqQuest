c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKZERO
c
c
      subroutine MKZERO( n, a )
c---------------------------------------------------------------
c Purpose: zero out a double precision real array
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(n)
      DATA  zero / 0.d0 /
c
      do  i=1,n
        a(i) = zero
      enddo
c
      RETURN
      END
