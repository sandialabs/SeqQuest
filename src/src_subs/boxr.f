c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BOXR
c
c
      subroutine BOXR( hh, ira0,irb0,irc0, idra,idrb,idrc, r0, rmid )
c---------------------------------------------------------------
c Purpose: give cartesian coord of box corner, and midpoint
c
c Written:  Peter A. Schultz, 1-May-1997, for v2.20 (boxes)
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input array: (mesh definition vectors)
      DIMENSION  hh(3,3)
c output array: (r-offset, r-midpoint)
      DIMENSION  r0(3),rmid(3)
c local:
      DATA  half / 0.5d0 /
c
c >>>> EXECUTABLE CODE:
c
      do  i=1,3
        r0(i) = ira0*hh(i,1) + irb0*hh(i,2) + irc0*hh(i,3)
        rmid(i) = (idra-1)*hh(i,1) + (idrb-1)*hh(i,2) + (idrc-1)*hh(i,3)
      enddo
c
      do  i=1,3
        rmid(i) = r0(i) + half*rmid(i)
      enddo
c
      RETURN
      END
