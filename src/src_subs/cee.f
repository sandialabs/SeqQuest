c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CEE
c
c
      FUNCTION CEE( l1,m1,l2,m2,k )
c---------------------------------------------------------------
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL  ODD
c statement function:
      ODD(int) = (AND(int,1).eq.1)
c
      pf = SQRT(  DBLE( (2*l1+1)*(2*l2+1) ) )
      if( ODD(m1) ) pf = -pf
c
      cee = CEE1( l1,l2,k, -m1,m2,m1-m2 )*CEE1( l1,l2,k, 0,0,0 )*pf
c
      RETURN
      END
