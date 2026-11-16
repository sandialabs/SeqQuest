c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETMIJK
c
c
      subroutine GETMIJK( imatfl, norb,nk, dmat )
c---------------------------------------------------------------
c Purpose: retrieve full density (or analogous) matrix.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output array:
      DIMENSION  dmat(norb*norb,nk)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
        call READBIG( imatfl, norb*norb, dmat(1,k) )
      enddo
c
c    That's all Folks!
c
      RETURN
      END
