c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBDAT
c
c
      subroutine NEBDAT( IWR,IDAT,
     $ neb_step, nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
C     $ dt_neb,nsteps_neb,
     $ Lantikink_neb,Loptall_neb, 
     $ Lclimb_neb, Lnoninteract, Lrestart_neb,
     $ natm, imgstat, ratma,ratmb, ratm_neb )
c---------------------------------------------------------------
c Purpose: read definition of nudged elastic band
c
c Written: Aidan Thompson, Art Edwards, and Peter A. Schultz
c          December-2001
c
c Revision history:
c   9Apr15-PAS/2.66: lnwrite lines
c  11May05-RMV/2.59: add flags: climb-NEB, noninteracting images
c  20Jan03-PAS/2.54e: change keywords to set NEB flag options
c  10Jan02-PAS/2.54: adapted for lcao production code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL    Lantikink_neb, Loptall_neb
      LOGICAL    Lclimb_neb, Lnoninteract
      LOGICAL    Lrestart_neb
c  coordinates for initial, final, and intermediate NEB images
      DIMENSION  ratma(3,natm),ratmb(3,natm), ratm_neb(3,natm,*)
      DIMENSION  imgstat(*)
c local declarations:
      PARAMETER  ( nline = 128 )
      CHARACTER    label*(nline), readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a128)
 9000 format(1x,a)
 9010 format(1x,3f20.8)
c
c We use setdata input as initial start-point:
      img0stat = 1
c  ... and we do not yet have a final end-product:
      imgnstat = 0
c
c >>>> Required data for NEB run:
c
c Header for NEB section:
      wantlbl = 'neb da'
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
      if( readlbl.ne.wantlbl )then
        write(IWR,*) '>>>>> looking for neb data section'
        write(IWR,*) '>>>>> incorrect label'
        write(IWR,*) '>>>>> looking for: >neb data<'
        goto 1301
      endif
c
      wantlbl = 'images'
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
c The number of this NEB step
c (optional)
      wantlbl = 'neb st'
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  neb_step
        write(IWR,'(1x,i6)')  neb_step
        if( neb_step.lt.0 ) call STOPXERR( 'NEB restart: bad step' )
c
        wantlbl = 'images'
        read(IDAT,8000,err=1300,end=1300)  label
        call LNWRIT( IWR, label, nline )
      endif
c
c Number of images (and possibly atoms - ignored):
c (required)
      wantlbl = 'images'
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  nimg_neb
        write(IWR,'(1x,i6)')  nimg_neb
c
        if( nimg_neb.gt.nimg_max )then
          write(IWR,*) '>>>>> max number of neb images=',nimg_max
          call STOPXERR( 'nebdat - too many neb images' )
        endif
      else
        write(IWR,*) 'images - number of NEB images?'
        goto 1301
      endif
c
c
c Total energy of initial start-point
c (required)
      wantlbl = 'einiti'
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  ea_neb
        write(IWR,9010)  ea_neb
      else
        write(IWR,*) 'einitial - energy of initial image for NEB?'
        goto 1301
      endif
c
c Start-point coordinates:
c (optional - defaults to setdata input)
      wantlbl = 'coordi'
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
      call ATFMTIN( jatmfmt, label )
      if( jatmfmt.ge.0 )then
        call COORDRD( IDAT, jatmfmt, natm, ratma )
        call COORDWR( IWR ,          natm, ratma )
        img0stat = 1
c
        wantlbl = 'efinal'
        read(IDAT,8000,err=1300,end=1300)  label
        call LNWRIT( IWR, label, nline )
C     else
C       write(IWR,*) 'coordinates (x y z)?'
C       call STOPXERR( 'NEB restart: start-pt coords missing' )
      endif
c
c Total energy of end-product
c (required)
      wantlbl = 'efinal'
      if( readlbl.eq.wantlbl )then
        read(IDAT,*)  eb_neb
        write(IWR,9010)  eb_neb
      else
        write(IWR,*) 'efinal - energy of final image for NEB?'
        goto 1302
      endif
c
c End-product coordinates
c (required)
      wantlbl = 'coordi'
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
      call ATFMTIN( jatmfmt, label )
      if( jatmfmt.ge.0 )then
c       Get atomic coordinates for end-product:
        call COORDRD( IDAT, jatmfmt, natm, ratmb )
        call COORDWR( IWR ,          natm, ratmb )
        imgnstat = 1
      else
        write(IWR,*) 'atom, type, coordinates for endpoint of NEB?'
        goto 1302
      endif
c
c  ... and then parse *optional* inputs
c (required: an "end neb data" to signal end of neb data)
c
  100 continue
      wantlbl = 'end_ne'
      read(IDAT,8000,err=1301,end=1301)  label
      call LNWRIT( IWR, label, nline )
c
  101 continue
      if    ( readlbl.eq. 'end ne' .or.
     $        readlbl.eq. 'end_ne' )then
c       End of neb parameter input
        goto 990
c
C      elseif( readlbl.eq. 'timest' )then
Cc       Timestep size for nudged elastic band damped dynamics
C        read(IDAT,*)  dt_neb
C        write(IWR,'(1x,f15.8)')  dt_neb
Cc
C      elseif( readlbl.eq. 'nsteps' )then
Cc       Number of steps for nudged elastic band damped dynamics
C        read(IDAT,*)  nsteps_neb
C        write(IWR,'(1x,i8)')  nsteps_neb
Cc
      elseif( readlbl.eq. 'nebspr' )then
c       Spring constant for nudged elastic band
        read(IDAT,*)  spring_neb
        write(IWR,'(1x,f15.8)')  spring_neb
c
      elseif( readlbl.eq. 'antiki' )then
c       Flag to use anti-kink method (uphill-only tangent)
cpas:
        write(IWR,*) 'Antikink turned ON, to turn off use: no_antikink'
        Lantikink_neb = .true.
      elseif( readlbl.eq. 'no_ant' )then
c       Flag to turn off anti-kink method
        Lantikink_neb = .false.
c
cpas: note change of keyword structure for NEB
      elseif( readlbl.eq. 'optall' )then
c       Flag to drive minimal/maximal points downhill/uphill
        write(IWR,*) 'optall keyword changed, use: extreme|no_extreme'
        call STOPXERR( 'change optall => extreme|no_extreme' )
c
      elseif( readlbl.eq. 'extrem' )then
c       Push nearest images to minimum/maximum points along band
        Loptall_neb = .true.
      elseif( readlbl.eq. 'no_ext' )then
c       Do not push images to minimum/maximum points along band
        Loptall_neb = .false.
c
      elseif( readlbl.eq.'climb ' )then
c       Push high images to maxima (climbing image NEB)
        Lclimb_neb = .true.
c      elseif( readlbl.eq.'climbn' )then
c         Do not push high images to maxima
c         Lclimb_neb = .false.
c
      elseif( readlbl.eq.'nointer' )then
c       Non-interacting images
         Lnoninteract = .true.
c
      elseif( readlbl.eq. 'restar' )then
c       Use restrart of nudged elastic band
        Lrestart_neb = .true.
cpas:
        write(IWR,*) 'Restart turned ON, to turn off use: no_restart'
      elseif( readlbl.eq. 'no_res' )then
c       Prevent restart of nudged elastic band
        Lrestart_neb = .false.
c
      elseif( readlbl(1:5).eq. 'image' )then
c       Get geometry of one of the NEB images
        if( readlbl(6:6).eq. ' ' .or.
     $      readlbl(6:6).eq. 'I' )then
c         Coordinates for an intermediate image, which one?
          read(IDAT,*)  image
          write(IWR,'(1x,i4)')  image
          if( image.lt.1 .or. image.gt.nimg_neb )then
c           Number of image is outside the number in the NEB
            call STOPXERR( 'NEBDAT: image# is out-of-bounds' )
          endif
c
c       Following are non-standard, included only for completeness
        elseif( readlbl(6:6).eq. '0' )then
          image = 0
        elseif( readlbl(6:6).eq. 'N' )then
          image = -1
        else
          wantlbl = 'image '
          goto 1301
        endif
c
c       Find out format type of atom input line ...
        wantlbl = 'coordi'
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
        call ATFMTIN( jatmfmt, label )
        if( jatmfmt.lt.0 )then
          write(IWR,9000) 'coordinates ( x y z )?'
          goto 1301
        endif
c
c        ... and extract coordinates from input:
        if( image.gt.0 .and. image.le.nimg_neb )then
c         Intermediate image
          call COORDRD( IDAT, jatmfmt, natm, ratm_neb(1,1,image) )
          call COORDWR( IWR ,          natm, ratm_neb(1,1,image) )
          imgstat(image) = 1
        elseif( image.eq.0 )then
c         Initial-reactant positions
          call COORDRD( IDAT, jatmfmt, natm, ratma )
          call COORDWR( IWR ,          natm, ratma )
          img0stat = 1
        elseif( image.lt.0 )then
c         End-product positions
          call COORDRD( IDAT, jatmfmt, natm, ratmb )
          call COORDWR( IWR ,          natm, ratmb )
          imgnstat = 1
        endif
c
      else
        call STOPXERR( 'NEBDAT - unknown command' )
c
      endif
c
      goto 100
c
c >>>> Provide summary of NEB data:
c
  990 continue
      write(IWR,9991)  nimg_neb,  ea_neb,eb_neb,spring_neb,
C     $ dt_neb,nsteps_neb,
     $ Lantikink_neb,Loptall_neb,Lclimb_neb,Lrestart_neb
 9991 format(/1x,'NEB parameter summary:',
     $    ' Nimage, Estart, Eend, spring',
C     $    ', dt, nsteps',
     $ /   1x,i3,  3(1x,f15.8),
C     $     1x,f15.8,1x,i3,
     $ /   1x,'NEB options: antikink, optall, climb, restart=',
     $     4(1x,l8) / )
c
c Check for noninteracting image calculation
      if( Lnoninteract )then
         write(IWR,*)  '>>>> CAUTION: Noninteracting images selected.'
      endif
c Check status of endpoints:
      if( img0stat.eq.0 .or. imgnstat.eq.0 )
     $ call STOPXERR( 'NEB: either initial or final end-point missing' )
c
      RETURN
c
c >>>> ERROR handling
c
 1300 continue
      write(IWR,*)  'neb data: error/end reading neb data'
 1301 continue
      write(IWR,*) '>>>>> Want label >',wantlbl,'<, found >',readlbl,'<'
 1302 continue
      call STOPXERR( 'NEBDAT: format error in input file' )
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End
