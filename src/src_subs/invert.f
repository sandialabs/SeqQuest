c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INVERT
c
c
      subroutine INVERT( a,b,m, td,ad,bd, itmaxd )
c---------------------------------------------------------------
c Purpose: gaussian elimination
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c  21Dec99-PAS/2.42: comsetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(itmaxd,itmaxd),b(itmaxd,itmaxd)
      DIMENSION  td(itmaxd),ad(itmaxd),bd(itmaxd)
c
      DATA  zero,one / 0.d0,1.d0 /
c
      n = m
      if( n.gt.itmaxd )then
        call STOPXERR( 'amat2big/matrix too large in INVERT' )
      endif
c
      do 10 i=1,n
        atmp = a(i,i)
        if( ABS( atmp ) .lt. 1.0d-08 )then
          call STOPXERR( 'inv0diag/zero diag matrix in INVERT' )
        endif
  10  continue
c
      if( n.eq.1 ) goto 300
c
      do 100 i=1,n
c
        do  j=1,n
          td(j) = a(j,i)/a(i,i)
        enddo
        td(i) = zero
c
        do  k=1,n
          bd(k) = b(i,k)
          ad(k) = a(i,k)
        enddo
c
        do  k=1,n
          do  j=1,n
            b(j,k) = b(j,k) - (td(j)*bd(k))
            a(j,k) = a(j,k) - (td(j)*ad(k))
          enddo
        enddo
c
  100 continue
c
      do  i=1,n
        do  j=1,n
          b(j,i) = b(j,i)/a(j,j)
        enddo
      enddo
c
      RETURN
c
  300 continue
      b(1,1) = one/a(1,1)
c
      RETURN
      END
