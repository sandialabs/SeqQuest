c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> C16LOW
c
c
      subroutine C16LOW( n, rm, cm )
c--------------------------------------------------------------
c Purpose: take complex matrix loaded in r*8 packed array
c          and load into lower triangle of c*8 array
c          Also initialize other triangle to zero.
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c   3Oct04-PAS/2.59: zero the upper triangle, to fix parallel
c--------------------------------------------------------------
c
c Note: need to conjugate the existing matrix to get right phase
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0 )
      DIMENSION  rm(n,n), cm(2,n,n)
c
c >>>> EXECUTABLE CODE:
c
      do  j=1,n
        do  i=1,j-1
          cm(1,i,j) = zero
          cm(2,i,j) = zero
        enddo
        cm(1,j,j) = rm(j,j)
        cm(2,j,j) = zero
        do  i=j+1,n
          cm(1,i,j) =  rm(i,j)
          cm(2,i,j) = -rm(j,i)
        enddo
      enddo
c
      RETURN
      END
