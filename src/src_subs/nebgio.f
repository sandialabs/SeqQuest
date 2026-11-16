c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBRGET
c
c
      subroutine NEBRGET( IWR,IDAT,
     $ neb_step, nimg_neb, ea_neb,eb_neb,
     $ natm, imgstat, ratma,ratmb, ratm_neb )
c---------------------------------------------------------------
c Purpose: get geometries out of NEB configuration file
c
c Written: Aidan Thompson and Art Edwards, Dec-2001
c
c Revision history:
c  14Jun05-PAS/2.59: IERR changes
c  10Jul03-PAS/2.56: modify some output formats/labels
c  10Jan02-PAS/2.54: adapted for SeqQuest production code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c  coordinates for initial, final, and intermediate NEB images
      DIMENSION  ratma(3,natm),ratmb(3,natm), ratm_neb(3,natm,*)
      DIMENSION  imgstat(*)
c
c local declarations:
      CHARACTER    label*72, readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
      LOGICAL      printon
      DATA         lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a72)
 9000 format(a)
 9030 format(1x,3f20.8)
c
      printon = .false.
      if( IWR.ge.0 .and. lstout.gt.1 ) printon = .true.
c
      if( printon ) write(IWR,*) '>>>> READ NEB GEOMETRY FILE'
c Initialize to empty NEB images:
      do  img=1,nimg_neb
        imgstat(img) = 0
      enddo
c
      REWIND( IDAT )
c
c Header for NEB section
      wantlbl = 'neb da'
      read(IDAT,8000,err=1300,end=1300)  label
      if( printon ) write(IWR,9000)  label
c
      if( readlbl.ne.wantlbl )then
        goto 1301
      endif
c
c The number of this NEB step
      wantlbl = 'neb st'
      read(IDAT,8000,err=1300,end=1300)  label
      if( printon ) write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
c
        read(IDAT,*) neb_step
        if( printon ) write(IWR,'(2(1x,i6))') neb_step
        if( neb_step.lt.0 ) call STOPXERR( 'NEB restart: bad step' )
c
        wantlbl = 'images'
        read(IDAT,8000,err=1300,end=1300)  label
        if( printon ) write(IWR,9000)  label
      endif
c
c Number of images and atoms
      wantlbl = 'images'
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  nimgx, natmx
        if( printon ) write(IWR,'(2(1x,i6))') nimgx,natmx
c       Check these numbers against this calculation
        if( nimgx.ne.nimg_neb .or. natmx.ne.natm )
     $   call STOPXERR( 'NEB restart: # of images or atoms is wrong' )
      else
        goto 1301
      endif
c
c Total energy of initial start-point
      wantlbl = 'einiti'
      read(IDAT,8000,err=1300,end=1300)  label
      if( printon ) write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  ea_neb
        if( printon ) write(IWR,9030)  ea_neb
      else
        goto 1301
      endif
c
c Start-point coordinates
      wantlbl = 'coordi'
      read(IDAT,8000,err=1300,end=1300)  label
      if( printon ) write(IWR,9000)  label
c
      call ATFMTIN( jatmfmt, label )
      if( jatmfmt.ge.0 )then
        call COORDRD( IDAT, jatmfmt, natm, ratma )
        if( printon ) call COORDWR( IWR, natm, ratma )
c
        wantlbl = 'efinal'
        read(IDAT,8000,err=1300,end=1300)  label
        if( printon ) write(IWR,9000)  label
      else
        goto 1301
      endif
c
c Total energy of final end-product
      wantlbl = 'efinal'
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  eb_neb
        if( printon ) write(IWR,9030)  eb_neb
      else
        goto 1301
      endif
c
c End-product coordinates
      wantlbl = 'coordi'
      read(IDAT,8000,err=1300,end=1300)  label
      if( printon ) write(IWR,9000)  label
c
      call ATFMTIN( jatmfmt, label )
      if( jatmfmt.ge.0 )then
        call COORDRD( IDAT, jatmfmt, natm, ratmb )
        if( printon ) call COORDWR( IWR, natm, ratmb )
        imgnstat = 1
      else
        goto 1301
      endif
c
c Coordinates for *all* intermediate images
c
      do 200 img=1,nimg_neb
c
        wantlbl = 'image '
        read(IDAT,8000,err=990,end=990)  label
        if( printon ) write(IWR,9000)  label
        if( readlbl.ne.wantlbl ) goto 1301
c
c       Which intermediate image (do *not* assume images are in order)?
        read(IDAT,*)  image
        if( printon ) write(IWR,'(1x,i6)')  image
        if( image.lt.1 .or. image.gt.nimg_neb )then
c         Index of image is outside the number in the NEB
          call STOPXERR( 'NEB restart: # of image is out-of-bounds' )
        endif
c
c       Determine format type of atom input line ...
        wantlbl = 'coordi'
        read(IDAT,8000)  label
        if( printon ) write(IWR,9000)  label
        call ATFMTIN( jatmfmt, label )
        if( jatmfmt.lt.0 ) goto 1301
c
c        ... and extract coordinates from input:
        call COORDRD( IDAT, jatmfmt, natm, ratm_neb(1,1,image) )
        if( printon ) call COORDWR( IWR, natm, ratm_neb(1,1,image) )
        imgstat(image) = 1
c
  200 continue
c
  990 continue
c Check that picked up all the images:
      do  img=1,nimg_neb
        if( imgstat(img).ne.1 )
     $   call STOPXERR( 'NEB restart: bad image set' )
      enddo
c
      RETURN
c
c >>>> ERROR handling
c
 1300 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '>>>>> error/end reading neb restart'
 1301 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '>>>>> Want label >',wantlbl,'<,
     $                            found >',readlbl,'<'
 1302 continue
      call STOPXERR( 'NEB restart: format error in input file' )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBRSAV
c
c
      entry NEBRSAV( IDAT,
     $ neb_step, nimg_neb, ea_neb,eb_neb,
     $ natm, imgstat, ratma,ratmb, ratm_neb )
c---------------------------------------------------------------
c Purpose: write geometries to NEB configuration file
c
c Written: Peter Schultz, 16-January-2002, for v2.54
c
c Revision history:
c  16Jan02-PAS/2.54: adapted from Aidan Thompson original code
c---------------------------------------------------------------
c
c Header for NEB section
      write(IDAT,9000)  'neb data:'
c
c The sequence number for this step:
      write(IDAT,9000)  'neb step number'
      write(IDAT,'(2(1x,i6))')  neb_step
c
c Number of images and atoms
      write(IDAT,9000)  'images, and number of atoms'
      write(IDAT,'(2(1x,i6))') nimg_neb, natm
c
c Total energy of initial start-point
      write(IDAT,9000)  'einitial - energy for start-point'
      write(IDAT,9030)  ea_neb
c
c Start-point coordinates
      write(IDAT,9000)  'coordinates for start-point'
      call COORDWR( IDAT, natm, ratma )
c
c Total energy of final end-product
      write(IDAT,9000)  'efinal - energy for end-product'
      write(IDAT,9030)  eb_neb
c
c End-product coordinates
      write(IDAT,9000)  'coordinates for end-product'
      call COORDWR( IDAT, natm, ratmb )
c
c Coordinates for *all* intermediate images:
c
      do 1200 img=1,nimg_neb
        image = img
c
        write(IDAT,9000)  'image number in NEB'
        write(IDAT,'(1x,i6)')  image
c
        write(IDAT,9000)  'coordinates for this image'
        call COORDWR( IDAT, natm, ratm_neb(1,1,image) )
c
 1200 continue
c
c Mark the end of this file to this point
      write(IDAT,9000)  'end neb data'
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFGET
c
c
      subroutine NEBFGET( IWR, IDAT, nimg_neb,
     $ natm, imgstat, engyneb,
     $ fatmx, fatm_neb )
c---------------------------------------------------------------
c Purpose: retrieve NEB images energies and (raw) forces
c
c Written: Peter Schultz, 18-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c status of images:
      DIMENSION  imgstat(*)
c total energies of images:
      DIMENSION  engyneb(*)
c forces for everybody:
      DIMENSION  fatm_neb(3,natm,*)
c scratch declaration:
      DIMENSION  fatmx(3,natm)
c
c >>>> EXECUTABLE CODE:
c
c Count up the number of images we should have:
      nimg_f = 0
      do  img=1,nimg_neb
        if( imgstat(img).eq.2 ) nimg_f = nimg_f + 1
      enddo
c
c Get the images:
      do  img=1,nimg_f
        call NEBGETF( IWR,IDAT, image, eimage, natm, fatmx, ierror )
        if( ierror.ne.0 ) goto 1301
        if( imgstat(image).ne.2 ) goto 1302
c
c       Take in the energy and the forces:
        engyneb(image) = eimage
        do  iatm=1,natm
          fatm_neb(1,iatm,image) = fatmx(1,iatm)
          fatm_neb(2,iatm,image) = fatmx(2,iatm)
          fatm_neb(3,iatm,image) = fatmx(3,iatm)
        enddo
c
      enddo
c
      RETURN
c
 1301 call STOPXERR( 'NEBFGET: error reading force file' )
 1302 call STOPXERR( 'NEBFGET: image force from image without force' )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFPUT
c
c
      entry NEBFPUT( IDAT, nimg_neb,
     $ natm, img_stat, engyneb,
     $ fatm_neb )
c---------------------------------------------------------------
c Purpose: Write out all image forces.
c
c Written: Peter A. Schultz,  8-February-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  img=1,nimg_neb
        image = img
        eimage = engyneb(img)
        call NEBPUTF( IDAT, image,eimage, natm, fatm_neb(1,1,image) )
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFWIPE
c
c
      entry NEBFWIPE( IDAT )
c---------------------------------------------------------------
c Purpose: Reset NEB energy/force file to null
c
c Written: Peter A. Schultz,  8-February-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      REWIND( unit=IDAT )
      write(IDAT,'(a)') 'no image results yet for this NEB step'
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBFSTAT
c
c
      entry NEBFSTAT( IWR, IDAT, nimg_neb,  natm, imgstat,
     $ fatmx )
c---------------------------------------------------------------
c Purpose: find number of complete images (forces) in NEB restart
c
c Written: Peter Schultz, 18-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      REWIND( IDAT )
c
C     nimgf_neb = 0
      do  img=1,nimg_neb
        call NEBGETF( IWR,IDAT, image, eimage, natm, fatmx, ierror )
        if( ierror.ne.0 ) goto 900
c
C       nimgf_neb = nimgf_neb + 1
        imgstat(image) = 2
      enddo
c
  900 continue
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBGETF
c
c
      subroutine NEBGETF( IWR,IDAT, image, eimage, natm, fatm, ierror )
c---------------------------------------------------------------
c Purpose: read energy/force of image to NEB configuration file
c
c Written: Peter Schultz, 16-January-2002, for v2.54
c
c Revision history:
c  10Jul03-PAS/2.55a: fix an output format
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c  coordinates for initial, final, and intermediate NEB images
      DIMENSION  fatm(3,natm)
c
c local declarations:
      CHARACTER    label*72, readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
      LOGICAL      printon
      DATA         lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a72)
 9000 format(a)
 9011 format(1x,f20.8,a,i5)
 9030 format(1x,3f20.8)
c
      ierror = 0
      printon = .false.
      if( IWR.ge.0 .and. lstout.gt.0 ) printon = .true.
c
c Image number
      wantlbl = 'imagef'
      read(IDAT,8000,err=1310,end=1310)  label
      if( printon ) write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  image
        if( printon ) write(IWR,'(1x,i6)')  image
      else
        call STOPXERR( 'NEB restart: imagef' )
      endif
c
c Energy of image
      wantlbl = 'eimage'
      read(IDAT,8000)  label
      if( printon ) write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  eimage
        if( printon ) write(IWR,9030)  eimage
      else
        call STOPXERR( 'NEB restart: eimage' )
      endif
c
c Forces for all atoms in image
      wantlbl = 'forcei'
      read(IDAT,8000)  label
      if( printon ) write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
        jatmfmt = 0
        call COORDRD( IDAT, jatmfmt, natm, fatm )
        if( printon ) call COORDWR( IWR, natm, fatm )
      else
        call STOPXERR( 'NEB restart: forcei' )
      endif
c
      RETURN
c
c >>>> ERROR handling
c
 1310 ierror = 1
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBPUTF
c
c
      entry NEBPUTF( IDAT, image, eimage, natm, fatm )
c---------------------------------------------------------------
c Purpose: write energy/force of image to NEB configuration file
c
c Written: Peter Schultz, 16-January-2002, for v2.54
c
c Revision history:
c  10Jul03-PAS/2.55a: add grep-able string to output image energy
c---------------------------------------------------------------
c
c Image number
      write(IDAT,9000)  'imagef - image number for E/F output'
      write(IDAT,'(1x,i6)')  image
c
c Energy of image
      write(IDAT,9000)  'eimage - energy of image (Ry)'
      write(IDAT,9011) eimage,' Ry = CONVERGED ENERGY, NEB image#',image
c
c Forces for all atoms in image
      write(IDAT,9000)  'forcei - forces (Ry/bohr) on atoms'
      call COORDWR( IDAT, natm, fatm )
c
      RETURN
c
c    That's all Folks!
c
      END
