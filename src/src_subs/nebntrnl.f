c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBNTRNL
c
c
      subroutine NEBNTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $ nimg_neb, imgstat,  ratma,ratmb,
     $ ratm_neb,  ratm,ratm0 )
c---------------------------------------------------------------
c Purpose: scale neb images between internal/external coordinates
c
c Written: Peter A. Schultz, 6-February-2002, for v2.54
c
c Revision history:
c  20Jan03-PAS/2.54e: NEBXTRNL bug
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      LOGICAL  dolvatm
      DIMENSION  rprim(3,4),orig(3),rscale(3)
      DIMENSION  imgstat(*)
c Input/Output:
      DIMENSION  ratma(3,natm),ratmb(3,natm)
      DIMENSION  ratm_neb(3,natm,*)
c Scratch:
      DIMENSION  ratm(3,natm),ratm0(3,natm)
c
c >>>> EXECUTABLE CODE:
c
      nvec = 3*natm
c
c Convert start-point coordinates:
c      ... copy coords into temp work space
      call DCOPY( nvec, ratma,1, ratm,1 )
c      ... take external coords in ratm(), convert to internal ratm0()
      call RINTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $ ratm,ratm0 )
c      ... and move converted coordinates back into right place:
      call DCOPY( nvec, ratm0,1, ratma,1 )
c
c Convert end-point coordinates:
      call DCOPY( nvec, ratmb,1, ratm,1 )
      call RINTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $ ratm,ratm0 )
      call DCOPY( nvec, ratm0,1, ratmb,1 )
c
c Convert all initialized images:
      do  img=1,nimg_neb
        if( imgstat(img).gt.0 )then
          call DCOPY( nvec, ratm_neb(1,1,img),1, ratm,1 )
          call RINTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $     ratm,ratm0 )
          call DCOPY( nvec, ratm0,1, ratm_neb(1,1,img),1 )
        endif
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBXTRNL
c
c
      entry      NEBXTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $ nimg_neb, imgstat,  ratma,ratmb,
     $ ratm_neb,  ratm,ratm0 )
c---------------------------------------------------------------
c
      nvec = 3*natm
c
c Convert start-point coordinates:
c      ... copy coords into temp work space
      call DCOPY( nvec, ratma,1, ratm0,1 )
c      ... take internal coords in ratm(), convert to external ratm0()
      call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $ ratm,ratm0 )
c      ... and move converted coordinates back into right place:
      call DCOPY( nvec, ratm,1, ratma,1 )
c
c Convert end-point coordinates:
      call DCOPY( nvec, ratmb,1, ratm0,1 )
      call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $ ratm,ratm0 )
      call DCOPY( nvec, ratm,1, ratmb,1 )
c
c Convert all initialized images:
      do  img=1,nimg_neb
        if( imgstat(img).gt.0 )then
          call DCOPY( nvec, ratm_neb(1,1,img),1, ratm0,1 )
          call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $     ratm,ratm0 )
          call DCOPY( nvec, ratm,1, ratm_neb(1,1,img),1 )
        endif
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
