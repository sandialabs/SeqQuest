c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMRHO3
c
c
      subroutine SYMRHO3( n1r,n2r,n3r, isymv,nsymv, rho )
c---------------------------------------------------------------
c Purpose: 3D-symmetrize density, done in-place in real space
c
c Written: Peter A. Schultz, 3-November-1994
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c symmetry point mapping arrays:
      DIMENSION  isymv(*),nsymv(*)
c density on grid
      DIMENSION  rho(*)
c local declarations:
      DATA  one / 1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      do  ir=1,n1r*n2r*n3r
        ifst = isymv(ir)
        symnum = one/DBLE(nsymv(ifst))
        if( ifst.eq.ir )then
          rho(ifst) = rho(ifst)*symnum
        else
          rho(ifst) = rho(ifst) + rho(ir)*symnum
        endif
      enddo
c
c  Have average density in first of each symmetry set,
c  now broadcast to other members of symmetry sets:
c
      do  ir=1,n1r*n2r*n3r
        ifst = isymv(ir)
        rho(ir) = rho(ifst)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
