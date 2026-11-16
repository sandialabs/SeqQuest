c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKSZERO
c
c
      subroutine MKSZERO( n, a )
c---------------------------------------------------------------
c Purpose: zero out a single precision real array
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  a,zero
c
      DIMENSION  a(n)
      DATA  zero / 0.e0 /
c
      do  i=1,n
        a(i) = zero
      enddo
c
      RETURN
      END
