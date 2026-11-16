c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMRHO2
c
c
      subroutine SYMRHO2( n1r,n2r,n3r, isymv,nsymv, rho )
c---------------------------------------------------------------
c Purpose: 2D-symmetrize charge density (in real space)
c
c Written: Peter A. Schultz, 20-February-1990
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  Note: symmetrization of rho() is done in-place, the
c        density of the symmetry-related mesh points
c        being accumulated in the first mesh point of
c        a symmetry set and then being broadcast to
c        the other members of that set afterwards.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  isymv(*),nsymv(*)
c input/output arrays:
      DIMENSION  rho(n1r*n2r,n3r)
c local declarations:
      DATA  one / 1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      do 100 ir=1,n1r*n2r
        ifst = isymv(ir)
        symnum = one/ DBLE(nsymv(ifst))
c
        if( ifst.eq.ir )then
          do  i3=1,n3r
            rho(ifst,i3) = rho(ifst,i3)*symnum
          enddo
        else
          do  i3=1,n3r
            rho(ifst,i3) = rho(ifst,i3) + rho(ir,i3)*symnum
          enddo
        endif
c
  100 continue
c
c  Have average density in first of each symmetry set,
c  now broadcast to other members of symmetry sets:
c
      do  ir=1,n1r*n2r
        ifst = isymv(ir)
        do  i3=1,n3r
          rho(ir,i3) = rho(ifst,i3)
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
