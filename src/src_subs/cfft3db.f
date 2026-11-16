c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CFFT3DB
c
c
      subroutine CFFT3DB( array, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c---------------------------------------------------------------
c Purpose: compute normalized complex backward fourier transform
c          of the 3-dimensionally periodic "array."
c
c Revision history:
c   5Dec07-PAS/2.60c: replace cmplx<->2dble with more generic
c  21Dec99-PAS/2.42: comsetic clean-up
c   9Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  For case n1r, n2r, or n3r odd, phase factor multiplication
c  translates in k-space by -(n1r-1)/2, -(n2r-1)/2, etc
c  For case n1r, n2r, or n3r even, phase factor multiplication
c  translates in k-space by -(n1r/2-1), -(n2r/2-1), etc
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,expamx
c
c input/output arrays:
      DOUBLE COMPLEX    array(n1r,n2r,n3r)
      DIMENSION  ws1(*),ws2(*),ws3(*), ftarray(*)
c
c local declarations:
      DIMENSION  phx(2),phy(2),phz(2)
C      EQUIVALENCE  (phx,ph1),(phy,ph2),(phz,ph3)
      DOUBLE COMPLEX    ph1,ph2,ph3, ph3j,ph32j,ph321j
      LOGICAL    EVEN
      DATA  zero,one / 0.d0,1.d0 /
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      do  j2=1,n2r
        do  j3=1,n3r
          call CFFTB( n1r, array(1,j2,j3),ws1 )
        enddo
      enddo
c
      do  j1=1,n1r
        do  j3=1,n3r
          call ZCOPY( n2r, array(j1,1,j3),n1r, ftarray,1 )
          call CFFTB( n2r, ftarray,ws2 )
          call ZCOPY( n2r, ftarray,1, array(j1,1,j3),n1r )
        enddo
      enddo
c
      do  j1=1,n1r
        do  j2=1,n2r
          call ZCOPY( n3r, array(j1,j2,1),n1r*n2r, ftarray,1 )
          call CFFTB( n3r, ftarray,ws3 )
          call ZCOPY( n3r, ftarray,1, array(j1,j2,1),n1r*n2r )
        enddo
      enddo
c
c      restore phase coming from translation in k-space (see cfft3d)
c      in reciprocal space
c
c 18Dec92-PAS: cexp()-> generics; 5Dec07-PAS: even more generic
      if( EVEN(n1r) )then
        ph1 = -one
      else
c       ph1  =-cexp(cmplx(0.,pi/float(n1r)))
        f1r = DBLE(n1r)
        phx(1) = -COS( pi/f1r )
        phx(2) = -SIN( pi/f1r )
        ph1 = DCMPLX( phx(1), phx(2) )
      endif
c
      if( EVEN(n2r) )then
        ph2 = -one
      else
c       ph2 = -cexp(cmplx(0.,pi/float(n2r)))
        f2r = DBLE(n2r)
        phy(1) = -COS( pi/f2r )
        phy(2) = -SIN( pi/f2r )
        ph2 = DCMPLX( phy(1), phy(2) )
      endif
c
      if( EVEN(n3r) )then
        ph3 = -one
      else
c       ph3 = -cexp(cmplx(0.,pi/float(n3r)))
        f3r = DBLE(n3r)
        phz(1) = -COS( pi/f3r )
        phz(2) = -SIN( pi/f3r )
        ph3 = DCMPLX( phz(1), phz(2) )
      endif
c
      ph3j = one
      do  j3=1,n3r
        ph32j = ph3j
        do  j2=1,n2r
          ph321j = ph32j
          do  j1=1,n1r
            array(j1,j2,j3) = array(j1,j2,j3)*ph321j
            ph321j = ph1*ph321j
          enddo
          ph32j = ph2*ph32j
        enddo
        ph3j = ph3*ph3j
      enddo
c
      RETURN
      END
