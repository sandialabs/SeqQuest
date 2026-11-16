c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMFRC
c
c
      subroutine SYMFRC( IWRX, rtnlbl, stress, frcsym, frctmp, fdefct,
     $ ndim, natm, nsym, rprim, rmatsym, naofsym )
c---------------------------------------------------------------
c Purpose: symmetrize forces
c
c Written: Peter A. Schultz, 4-November-1994
c
c Revision history:
c  14Mar13-PAS/2.64: write format for 1000+ atoms
c   4Oct04-PAS/2.59: add grep-able pressure output diagnostic
c  26Jun02-PAS/2.54: cosmetic - to emphasize IWRX can be zero
c   1May02-PAS/2.53: push stress units conversion to subroutine
c  21Jun01-PAS/2.48: replace STOPs
c   2May01-PAS/2.47: convert stress units into standard
c  12Apr00-PAS/2.43: specify units in stress output
c  28Sep99-PAS/2.40: full three-dimensional stress tensor;
c                    change output formats.
c   5Jan99-PAS/2.30: redo symmetry scheme
c  17Dec98-PAS/2.30: clean force/stress arrays
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  26Sep97-PAS/2.21: force defect calculation
c  28Apr97-PAS/2.17: naming conventions
c   3Sep96-PAS/2.15: fix rotmat indice bug
c---------------------------------------------------------------
c
c ***** Note: IWRX that is passed in can be passed in as zero!
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output arrays:
      CHARACTER  rtnlbl*(*)
      DIMENSION  frcsym(3,natm), frctmp(3,natm), fdefct(3)
      DIMENSION  stress(3,3)
c  ... primitive unit cell vectors
      DIMENSION  rprim(3,3)
c  ... cartesian representations of non-identity symmetry ops
      DIMENSION  rmatsym(3,4,*)
c  ... atom images under inverse symmetry ops
      DIMENSION  naofsym(natm,*)
c local declarations
      DIMENSION  fi(3),sfi(3), strgpa(3,3)
      CHARACTER*8  strlbl

crpm  Create a common block to save the pressure for reference later
c     We are probably going to want to tear this code out and replace
c     with a subroutine at some point.
      COMMON /savepress/ presssave,prgpasave
c
      DATA  zero,one / 0.d0,1.d0 /
      DATA  vsmall / 1.d-6 /
c
c >>>> EXECUTABLE CODE:
c
c copy identity first:
      do  iatm=1,natm
        frctmp(1,iatm) = frcsym(1,iatm)
        frctmp(2,iatm) = frcsym(2,iatm)
        frctmp(3,iatm) = frcsym(3,iatm)
      enddo
c
      if( nsym.le.0 ) goto 500
c
c  Loop over all atoms, accumulate symmetry force
c
      do  iatm=1,natm
        do  isym=1,nsym
          jatm = naofsym(iatm,isym+1)
          fi(1) = frcsym(1,jatm)
          fi(2) = frcsym(2,jatm)
          fi(3) = frcsym(3,jatm)
c compute symmetry force:
          sfi(1) = rmatsym(1,1,isym)*fi(1) + rmatsym(2,1,isym)*fi(2)
     $           + rmatsym(3,1,isym)*fi(3)
          sfi(2) = rmatsym(1,2,isym)*fi(1) + rmatsym(2,2,isym)*fi(2)
     $           + rmatsym(3,2,isym)*fi(3)
          sfi(3) = rmatsym(1,3,isym)*fi(1) + rmatsym(2,3,isym)*fi(2)
     $           + rmatsym(3,3,isym)*fi(3)
c have jatm = S(iatm) or iatm=S**-1(jatm)
          frctmp(1,iatm) = frctmp(1,iatm) + sfi(1)
          frctmp(2,iatm) = frctmp(2,iatm) + sfi(2)
          frctmp(3,iatm) = frctmp(3,iatm) + sfi(3)
        enddo
      enddo
c
      symscal = one/ DBLE( nsym+1 )
      do  iatm=1,natm
        frctmp(1,iatm) = symscal*frctmp(1,iatm)
        frctmp(2,iatm) = symscal*frctmp(2,iatm)
        frctmp(3,iatm) = symscal*frctmp(3,iatm)
      enddo
c
c Symmetrize stress tensor:
      call SYMSTR( ndim, stress, nsym, rmatsym )
c
c  Summarize results:
c
  500 continue
      fdefct(1) = zero
      fdefct(2) = zero
      fdefct(3) = zero
      do  iatm=1,natm
c if negligible, set to zero:
        if( ABS(frctmp(1,iatm)).lt.vsmall ) frctmp(1,iatm) = zero
        if( ABS(frctmp(2,iatm)).lt.vsmall ) frctmp(2,iatm) = zero
        if( ABS(frctmp(3,iatm)).lt.vsmall ) frctmp(3,iatm) = zero
c put into right place:
        frcsym(1,iatm) = frctmp(1,iatm)
        frcsym(2,iatm) = frctmp(2,iatm)
        frcsym(3,iatm) = frctmp(3,iatm)
c accumulate defect:
        fdefct(1) = fdefct(1) + frctmp(1,iatm)
        fdefct(2) = fdefct(2) + frctmp(2,iatm)
        fdefct(3) = fdefct(3) + frctmp(3,iatm)
      enddo
c
      if( IWRX.gt.0 )then
c
        write(IWRX,*) ' '
c
        if( ndim.gt.0 )then
c
c Output stress result ...
c
c  ... convert stress into the appropriate units for dimension:
          call STR2GPA( ndim, rprim, stress, strgpa, strlbl )
c
c  ... and output the stress, both raw, and converted:
          write(IWRX,9550)
     $     rtnlbl,' stress tensor, in Rydberg and in ',strlbl
          do  jd=1,ndim
            if( ndim.eq.1 )then
              write(IWRX,9551)  (stress(id,jd),id=1,ndim),
     $                          (strgpa(id,jd),id=1,ndim)
            elseif( ndim.eq.2 )then
              write(IWRX,9552)  (stress(id,jd),id=1,ndim),
     $                          (strgpa(id,jd),id=1,ndim)
            elseif( ndim.eq.3 )then
              write(IWRX,9553)  (stress(id,jd),id=1,ndim),
     $                          (strgpa(id,jd),id=1,ndim)
            endif
          enddo
          press = zero
          prgpa = zero
          do  jd=1,ndim
            press = press + stress(jd,jd)
            prgpa = prgpa + strgpa(jd,jd)
          enddo
          press = press / DBLE( ndim )
          prgpa = prgpa / DBLE( ndim )
          write(IWRX,9555) press, prgpa,  strlbl
crpm  Here are the explicit save statements for the Pressure
c     and the GPA pressure.
          presssave = press
          prgpasave = prgpa
c
 9550     format(1x,a,a,a)
 9551     format(1x,2(4x, f14.6))
 9552     format(1x,2(4x,2f14.6))
 9553     format(1x,2(4x,3f14.6))
 9555     format(1x,'pressure=Tr(stress)/ndim=',f14.6,' Ry=',f14.6,1x,a)
c
        endif
c
c Output force results:
c
        write(IWRX,*) rtnlbl,' force contributions(Ry/bohr):'
        write(IWRX,9560)
        do  iatm=1,natm
          write(IWRX,9561)  iatm, ( frcsym(i,iatm), i=1,3 )
        enddo
        write(IWRX,9562)  fdefct
c
 9560   format(8x,'atom',7x,'x force',10x,'y force',10x,'z force')
 9561   format(6x, i5,1x,        3f17.10 )
 9562   format(6x,'f-defect',3x, 3f17.10 )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
