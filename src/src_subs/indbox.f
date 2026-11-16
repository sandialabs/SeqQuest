c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INDBOX
c
c
      subroutine INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $ ira0,irb0,irc0, idra,idrb,idrc )
c---------------------------------------------------------------
c Purpose: take box index and break into sub-indices
c
c Written: Peter A. Schultz, 23-April-1997, for v2.20 (boxes)
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  ndbox(2,2,3)
c
c >>>> EXECUTABLE CODE:
c
      nboxa1 = ndbox(1,1,1)
      nboxb1 = ndbox(1,1,2)
      nboxc1 = ndbox(1,1,3)
c The number of boxes in the three dimensions:
      nboxa = ndbox(1,1,1) + ndbox(2,1,1)
      nboxb = ndbox(1,1,2) + ndbox(2,1,2)
      nboxc = ndbox(1,1,3) + ndbox(2,1,3)
c
      iboxbc = (ibox-1)/nboxa
c The box index in the three dimensions, in grid of boxes
      iboxc = iboxbc/nboxb + 1
      iboxb = iboxbc - (iboxc-1)*nboxb + 1
      iboxa = ibox - iboxbc*nboxa
c
c Get box grid sides and corner index on global grid of points
      idra1 = ndbox(1,2,1)
      if( iboxa.le.nboxa1 )then
        idra = idra1
        ira0 = (iboxa-1)*idra1
      else
        idra = idra1 - 1
        ira0 = nboxa1*idra1 + (iboxa-nboxa1-1)*idra
      endif
c
      idrb1 = ndbox(1,2,2)
      if( iboxb.le.nboxb1 )then
        idrb = idrb1
        irb0 = (iboxb-1)*idrb1
      else
        idrb = idrb1 - 1
        irb0 = nboxb1*idrb1 + (iboxb-nboxb1-1)*idrb
      endif
c
      idrc1 = ndbox(1,2,3)
      if( iboxc.le.nboxc1 )then
        idrc = idrc1
        irc0 = (iboxc-1)*idrc1
      else
        idrc = idrc1 - 1
        irc0 = nboxc1*idrc1 + (iboxc-nboxc1-1)*idrc
      endif
c
c    That's all Folks!
c
      RETURN
      END
