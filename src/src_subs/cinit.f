c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CINIT
c
c
      subroutine CINIT( c )
c---------------------------------------------------------------
c Purpose: initialize c() for electron-nucleus matrix element
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  c(10,9)
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
      do  j=1,9
        do  i=1,10
          c(i,j) = zero
        enddo
      enddo
c
      do  i=1,8
        c(i,i) = one
      enddo
      c(9,8)  = -one
      c(8,9)  = -one
      c(9,9)  = -one
      c(10,9) =  two
c
      RETURN
      END
