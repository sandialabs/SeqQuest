c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INTOBOX
c
c
      subroutine INTOBOX( n1r,n2r,n3r, vgrid, ndbox, vbox )
c---------------------------------------------------------------
c Purpose: copy grid field into boxes
c
c Written: Peter A. Schultz, 21-April-1997, for v2.20 (boxes)
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  vgrid(n1r,n2r,n3r), vbox(*)
      DIMENSION  ndbox(2,2,3)
c
c >>>> EXECUTABLE CODE:
c
      irbox = 0
c
      nboxa1 = ndbox(1,1,1)
      nboxb1 = ndbox(1,1,2)
      nboxc1 = ndbox(1,1,3)
      nboxa = ndbox(1,1,1) + ndbox(2,1,1)
      nboxb = ndbox(1,1,2) + ndbox(2,1,2)
      nboxc = ndbox(1,1,3) + ndbox(2,1,3)
      idra1 = ndbox(1,2,1)
      idrb1 = ndbox(1,2,2)
      idrc1 = ndbox(1,2,3)
c
      irc0 = 0
      idrc = idrc1
      do 103 iboxc=1,nboxc
        if( iboxc.gt.nboxc1 ) idrc = idrc1 - 1
c
        irb0 = 0
        idrb = idrb1
        do 102 iboxb=1,nboxb
          if( iboxb.gt.nboxb1 ) idrb = idrb1 - 1
c
          ira0 = 0
          idra = idra1
          do 101 iboxa=1,nboxa
            if( iboxa.gt.nboxa1 ) idra = idra1 - 1
c
            do  irc=irc0+1,irc0+idrc
              do  irb=irb0+1,irb0+idrb
                do  ira=ira0+1,ira0+idra
                  irbox = irbox + 1
                  vbox(irbox) = vgrid(ira,irb,irc)
                enddo
              enddo
            enddo
c
            ira0 = ira0 + idra
  101     continue
          irb0 = irb0 + idrb
  102   continue
        irc0 = irc0 + idrc
  103 continue
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FROMBOX
c
c
      ENTRY FROMBOX( n1r,n2r,n3r, vgrid, ndbox, vbox )
c---------------------------------------------------------------
c Purpose: get grid field from boxes
c
c Written: Peter A. Schultz, 21-April-1997, for v2.20
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c---------------------------------------------------------------
c
      irbox = 0
c
      nboxa1 = ndbox(1,1,1)
      nboxb1 = ndbox(1,1,2)
      nboxc1 = ndbox(1,1,3)
      nboxa = ndbox(1,1,1) + ndbox(2,1,1)
      nboxb = ndbox(1,1,2) + ndbox(2,1,2)
      nboxc = ndbox(1,1,3) + ndbox(2,1,3)
      idra1 = ndbox(1,2,1)
      idrb1 = ndbox(1,2,2)
      idrc1 = ndbox(1,2,3)
c
      irc0 = 0
      idrc = idrc1
      do 203 iboxc=1,nboxc
        if( iboxc.gt.nboxc1 ) idrc = idrc1 - 1
c
        irb0 = 0
        idrb = idrb1
        do 202 iboxb=1,nboxb
          if( iboxb.gt.nboxb1 ) idrb = idrb1 - 1
c
          ira0 = 0
          idra = idra1
          do 201 iboxa=1,nboxa
            if( iboxa.gt.nboxa1 ) idra = idra1 - 1
c
            do  irc=irc0+1,irc0+idrc
              do  irb=irb0+1,irb0+idrb
                do  ira=ira0+1,ira0+idra
                  irbox = irbox + 1
                  vgrid(ira,irb,irc) = vbox(irbox)
                enddo
              enddo
            enddo
c
            ira0 = ira0 + idra
  201     continue
          irb0 = irb0 + idrb
  202   continue
        irc0 = irc0 + idrc
  203 continue
c
c    That's all Folks!
c
      RETURN
      END
