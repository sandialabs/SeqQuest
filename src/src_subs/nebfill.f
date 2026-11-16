c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFILL
c
c
      subroutine NEBFILL( IWR, natm, nimg_neb,imgstat,
     $ ratma,ratmb, ratm_neb, ratmx )
c---------------------------------------------------------------
c Purpose: fill in missing NEB images
c
c Written: Peter Schultz, 14-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c The current code uses linear interpolation to fill in missing
c images in the NEB.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/output declarations ...
c  ... coordinates for initial, final, and intermediate NEB images
      DIMENSION  ratma(3,natm),ratmb(3,natm), ratm_neb(3,natm,*)
c  ... existence of coordinates for given image
      DIMENSION  imgstat(*)
c Scratch:
      DIMENSION  ratmx(3,natm,2)
c Local:
      DATA       zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Initialize lead point coordinates with initial start-point:
      do  iatm=1,natm
        ratmx(1,iatm,1) = ratma(1,iatm)
        ratmx(2,iatm,1) = ratma(2,iatm)
        ratmx(3,iatm,1) = ratma(3,iatm)
      enddo
      imglast = 0
c
  100 continue
c
c Find the first empty image:
c (and, thereby, the filled image "img0" at start of open interval)
c
      do 150 img=imglast+1,nimg_neb
        if( imgstat(img).eq.0 )then
          img0 = img - 1
          goto 200
        endif
  150 continue
c     All images already determined - done.
      goto 900
c
  200 continue
      if( img0.gt.0 )then
c       Transfer coordinates of filled image at start of open interval:
        do  iatm=1,natm
          ratmx(1,iatm,1) = ratm_neb(1,iatm,img0)
          ratmx(2,iatm,1) = ratm_neb(2,iatm,img0)
          ratmx(3,iatm,1) = ratm_neb(3,iatm,img0)
        enddo
      endif
c
c Find filled image "imgn" at end of open interval:
c
      do 250 img=img0+1,nimg_neb
        if( imgstat(img).gt.0 )then
c         Found filled image at end of open interval
          imgn = img
c         Transfer coords of filled image
          do  iatm=1,natm
            ratmx(1,iatm,2) = ratm_neb(1,iatm,imgn)
            ratmx(2,iatm,2) = ratm_neb(2,iatm,imgn)
            ratmx(3,iatm,2) = ratm_neb(3,iatm,imgn)
          enddo
          goto 300
        endif
  250 continue
c     Did not find filled final image in band, use final end-product:
      imgn = nimg_neb + 1
      do  iatm=1,natm
        ratmx(1,iatm,2) = ratmb(1,iatm)
        ratmx(2,iatm,2) = ratmb(2,iatm)
        ratmx(3,iatm,2) = ratmb(3,iatm)
      enddo
c
  300 continue
c
c Fill in the empty images in the band (linear interpolation):
c
      nebgap = imgn - img0
      gapint = one / DBLE( nebgap )
      gapneb = zero
c
      do  img=img0+1,imgn-1
        imgstat(img) = 1
        gapneb = gapneb + gapint
        do  iatm=1,natm
          ratm_neb(1,iatm,img) = ratmx(1,iatm,1) +
     $     gapneb*( ratmx(1,iatm,2) - ratmx(1,iatm,1) )
          ratm_neb(2,iatm,img) = ratmx(2,iatm,1) +
     $     gapneb*( ratmx(2,iatm,2) - ratmx(2,iatm,1) )
          ratm_neb(3,iatm,img) = ratmx(3,iatm,1) +
     $     gapneb*( ratmx(3,iatm,2) - ratmx(3,iatm,1) )
        enddo
      enddo
c
      imglast = imgn
c Go back and try to find other open intervals in the band:
      if( imglast.lt.nimg_neb ) goto 100
c
  900 continue
c
c    That's all Folks!
c
      RETURN
      END
