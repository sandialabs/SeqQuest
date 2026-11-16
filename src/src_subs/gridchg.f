c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDCHG
c
c
      subroutine GRIDCHG( chg,dum,weight,nptr,arr )
c---------------------------------------------------------------
c Purpose: count up movement on grid of grid quantity
c
c Written: Peter A. Schultz, 16-June-1993
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  arr(nptr)
c
      chg = 0.d0
      do  nr=1,nptr
        chg = chg + ABS(arr(nr))
      enddo
      chg = weight*chg
c
      RETURN
      END
