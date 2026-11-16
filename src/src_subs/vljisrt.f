c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLJISRT
c
c
      subroutine VLJISRT( norb,nk, vto, vfrom )
c---------------------------------------------------------------
c Purpose: copy matrices from k-first to k-last form
c          Not very efficient, but I do not care.
c
c Written: Peter A. Schultz, 15-September-1999, for v2.39
c
c Revision History:
c  27Mar00-PAS/2.43: rename "enjisrt" to "vljisrt"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output array:
      DIMENSION  vto(norb,norb,nk)
c input array:
      DIMENSION  vfrom(nk,norb,norb)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
        do  j=1,norb
          do  i=1,norb
            vto(i,j,k) = vfrom(k,j,i)
          enddo
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
