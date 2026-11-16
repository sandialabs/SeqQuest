c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TREIGVEC
c
c
      subroutine TREIGVEC( norb,nstate, veceig,eigvec )
c---------------------------------------------------------------
c Purpose: transpose complex eigenvector array.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
c  23Jul01-PAS/2.49: push down from rhoij
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input array:
      DIMENSION  veceig(2,norb,nstate)
c output arrays:
      DIMENSION  eigvec(2,nstate,norb)
c
c >>>> EXECUTABLE CODE:
c
      do  n=1,nstate
        do  i=1,norb
          eigvec(1,n,i) = veceig(1,i,n)
          eigvec(2,n,i) = veceig(2,i,n)
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
