c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CINITF
c
c
      subroutine CINITF( c,gradc )
c---------------------------------------------------------------
c Purpose: initialize c() for electron-nucleus force element
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  c(10,9),gradc(3,10,9)
      DATA  zero,one,two,four / 0.d0,1.d0,2.d0,4.d0 /
c
      do  j=1,9
        do  i=1,10
          c(i,j) = zero
          gradc(1,i,j) = zero
          gradc(2,i,j) = zero
          gradc(3,i,j) = zero
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
      gradc(1,1,2) = -one
      gradc(2,1,3) = -one
      gradc(3,1,4) = -one
      gradc(2,2,5) = -one
      gradc(1,3,5) = -one
      gradc(3,3,6) = -one
      gradc(2,4,6) = -one
      gradc(3,2,7) = -one
      gradc(1,4,7) = -one
      gradc(1,2,8) = -two
      gradc(2,3,8) =  two
      gradc(1,2,9) =  two
      gradc(2,3,9) =  two
      gradc(3,4,9) = -four
c
      RETURN
      END
