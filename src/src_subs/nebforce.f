c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFORCE
c
c
      subroutine NEBFORCE( IWRNEB,  Lantikink, Loptall,
     $ Lclimb, Lnoninteract,
     $ nimg_neb, spring_neb,
     $ natm, imgstat, engyneb, ea_neb, eb_neb,
     $ fatm_neb, ratm_neb, ratma, ratmb,
     $ rilo, rihi, tmpa, pathout ) 
c---------------------------------------------------------------
c Purpose: generate band forces for NEB calculation
c
c Written: Peter A. Schultz, 13-February-2002, for v2.54
c          Based on original code by Aidan Thompson (SNL)
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c  11May05-RMV/2.59: add climb, noninteracting, analysis
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
c   Option flags for NEB:
      LOGICAL    Lantikink, Loptall, Lclimb
      LOGICAL    Lnoninteract
c   Status and energies for images
      DIMENSION  imgstat(*), engyneb(*)
c   Coordinates for the initial-point, and end-product
      DIMENSION  ratma(3,natm),ratmb(3,natm)
c   Coordinates and forces for all images
      DIMENSION  ratm_neb(3,natm,nimg_neb)
c   Energy, distance and force projections for each image
      DIMENSION pathout(5,nimg_neb)
c
c Input/Output:
c   Forces over all images
      DIMENSION  fatm_neb(3,natm,nimg_neb)
c
c Scratch space:
      DIMENSION  rilo(3,natm),rihi(3,natm)
      DIMENSION  tmpa(natm,*)
c
c Local declarations:
      LOGICAL    Lnoenergy
      LOGICAL    Loptalltmp,Lantikinktmp,Lclimbtmp
      REAL*8 pathlength
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Set up the various NEB options flags:
c
      Lnoenergy = .false.
c
      if( .not.Loptall )then
        Loptalltmp = .false.
      elseif( Lnoenergy )then
        Loptalltmp = .false.
      else
        Loptalltmp = .true.
      endif
c
      if( .not.Lclimb )then
        Lclimbtmp = .false.
      elseif( Lnoenergy )then
        Lclimbtmp = .false.
      else
        Lclimbtmp = .true.
      endif
c
      if( .not.Lantikink )then
        Lantikinktmp = .false.
      elseif( Lnoenergy )then
        Lantikinktmp = .false.
      else
        Lantikinktmp = .true.
      endif
c
c Cycle over all the images:
c
c     Set initial total path distance to zero
      pathlength = zero
      do 1000 img=1,nimg_neb
c
        engyimg = engyneb(img)
c
c       Get front and back-end coordinates:
c
        if( img.eq.1 )then
c         Load initial-state into down-chain coords
          do  iatm=1,natm
            rilo(1,iatm) = ratma(1,iatm)
            rilo(2,iatm) = ratma(2,iatm)
            rilo(3,iatm) = ratma(3,iatm)
          enddo
          engylo = ea_neb
        else
c         Load next lower image into down-chain coords
          imglo = img - 1
          do  iatm=1,natm
            rilo(1,iatm) = ratm_neb(1,iatm,imglo)
            rilo(2,iatm) = ratm_neb(2,iatm,imglo)
            rilo(3,iatm) = ratm_neb(3,iatm,imglo)
          enddo
          engylo = engyneb(imglo)
        endif
c
        if( img.eq.nimg_neb )then
c         Load final-state into up-chain coords
          do  iatm=1,natm
            rihi(1,iatm) = ratmb(1,iatm)
            rihi(2,iatm) = ratmb(2,iatm)
            rihi(3,iatm) = ratmb(3,iatm)
          enddo
          engyhi = eb_neb
        else
c         Load next higher image into up-chain coords
          imghi = img + 1
          do  iatm=1,natm
            rihi(1,iatm) = ratm_neb(1,iatm,imghi)
            rihi(2,iatm) = ratm_neb(2,iatm,imghi)
            rihi(3,iatm) = ratm_neb(3,iatm,imghi)
          enddo
          engyhi = engyneb(imghi)
        endif
c     
c       Non-interacting image option can just return from here
        if( Lnoninteract ) goto 1000
c
c       Apply the band force for this image along the band:
c
        call FORCE_NEB( IWRNEB,
     $     natm, rilo, ratm_neb(1,1,img),
     $     rihi, fatm_neb(1,1,img),
     $     tmpa(1,1), tmpa(1,4),     tmpa(1,7),
c -->      ti(3a)-s   driminus(3a)-s driplus(3a)-s
     $     spring_neb,
     $     engylo,engyimg,engyhi,   drsqplus,drsqminus,
     $     fpt,fpn,fet,
     $     Lantikinktmp,Loptalltmp,Lclimbtmp,
     $      fdottan,dprev,dnext)
c
c       Load array of projected force,dprev,dnext,
c        and cumulative path distance (RMV)
        pathlength = pathlength + dprev
        pathout( 1,img ) = pathlength
        pathout( 2,img ) = engyimg
        pathout( 3,img ) = fdottan
        pathout( 4,img ) = dprev
        pathout( 5,img ) = dnext
c
c       End loop over images
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
