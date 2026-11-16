c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKEIKL
c
c
      subroutine MKEIKL( nk,nlat, vk,rlat,  coskr,sinkr )
c---------------------------------------------------------------
c Purpose: computes Bloch phase factors for lattice vectors
c
c Written: Peter A. Schultz, 27-March-2000, for v2.43
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  vk(3,nk), rlat(3,nlat)
c output arrays:
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
c local declarations:
      DATA  zero,one / 0.d0,1.d0 /
      cutoff = 1d-12
c
c >>>> EXECUTABLE CODE:
c
      do  ilat=1,nlat
c
        do  k=1,nk
c
          dot = vk(1,k)*rlat(1,ilat) + vk(2,k)*rlat(2,ilat)
     $        + vk(3,k)*rlat(3,ilat)
c
          if( dot .eq. zero )then
            coskr(k,ilat) = one
            sinkr(k,ilat) = zero
          else
            coskr(k,ilat) = COS( dot )
            sinkr(k,ilat) = SIN( dot )
          endif
c     RPM put in cutoff to make sure small values don't mess up
c     the floating point math:
          if (abs(coskr(k,ilat)) .lt. cutoff)
     $         coskr(k,ilat) = zero
          if (abs(sinkr(k,ilat)) .lt. cutoff)
     $         sinkr(k,ilat) = zero
c     
c
        enddo
c
      enddo
c
c    That's all Folks!
c
      RETURN
      END
