c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> YRECUR
c
c
      subroutine YRECUR( clp1m,clm,l,m )
c---------------------------------------------------------------
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c      produces coefficients for ylp2m=clp1m*ylp1m+clm*ylm
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      clp1m = SQRT( DBLE((2*l+5)*(2*l+3))/ DBLE((l+2)**2-m**2))
      clm = clp1m* SQRT( DBLE((l+1)**2-m**2)/ DBLE((2*l+1)*(2*l+3)))
c
      RETURN
      END
