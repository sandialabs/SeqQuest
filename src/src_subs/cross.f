c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CROSS
c
c
      subroutine CROSS( a,b, axb )
c---------------------------------------------------------------
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c axb() = a() X b()
      DIMENSION  a(3),b(3), axb(3)
c
      axb(1) = a(2)*b(3) - a(3)*b(2)
      axb(2) = a(3)*b(1) - a(1)*b(3)
      axb(3) = a(1)*b(2) - a(2)*b(1)
c
      RETURN
      END
