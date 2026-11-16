      PROGRAM BONDS
c***********************************************************************
c Program: BONDS
c
c Purpose: bond distance analysis for general unit cell
c
c Written: Peter A. Schultz
c
c Revision history:
c  17Apr14-PAS/1.2 : vrs markings, harmonized with "cells"
c   3Dec04-PAS/1.12: fix bug when have 0-coordinated atom
c  26Jul04-PAS/1.11: free-format coord system input
c   8May03-PAS/B1.1: skip from title to dimension on tape read
c   7Nov02-PAS/2.54: general read
c  19Jul00-PAS: updated for new code.
c   3Jan00-PAS: general scaling factors and non-cubic unit cell
c***********************************************************************
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Maximum number of atoms in unit cell:
      PARAMETER  (maxatm=2001)
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'1.2  '/, lastup/'17Apr14-PAS'/
c Atomic coordinates:
      DIMENSION  ratm(3,maxatm), ratm0(3,maxatm)
c Primitive unit cell vectors:
      DIMENSION  rprim(3,4), orig(3)
c Logical
      LOGICAL    dolvatm
c Unified scaling factors:
      DIMENSION  rscale(3)
c
      CHARACTER  opt*1
      CHARACTER*80  geomfil,rngfil
c I/O unit numbers:
      DATA  IRD,IWR,IGIN,ILST / 5,6,7,8 /
c Conversion to Angstrom from bohr:
      DATA  Aperau  / 0.5291771d0 /
c Default unit cell vectors:
      DATA  rprim / 1.d0, 0.d0, 0.d0,
     $              0.d0, 1.d0, 0.d0,
     $              0.d0, 0.d0, 1.d0,
     $              1.d0, 1.d0, 1.d0 /
c Default origin
      DATA  orig / 0.d0,0.d0,0.d0 /
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Attribute the blame:
      write(IWR,9000)  version, lastup,   maxatm
 9000 format(1x,72('=')
     $     /15x,'Program BONDS: atomic bond distance program'
C     $     /15x,'Written by: Peter A. Schultz'
     $     /15x,'Version:  ',a,'     Last revised: ',a
     $      /1x,72('=')
     $     /15x,'Maximum number of atoms = ',i5/ )
c
c Set a bunch of unit cell defaults:
      isetuc = 0
      ndim = 3
      dolvatm = .false.
      do  i=1,3
        rscale(i) = one
      enddo
c
c  Try to get a tape file with cell info.
c  If not tape file given, go to options.
c  If get a tape file, skip options and go directly to coordinates
c
      iopt = 1
   90   continue
        write(IWR,*)  'Enter name of tape file containing cell info:'
        read(IRD,8000)  geomfil
        if( geomfil .eq. ' ' .or. geomfil .eq. 'opt' ) goto 100
        OPEN(IGIN,file=geomfil,form='formatted',status='old',err=100)
        write(IWR,*) 'tape file found'
        call GETCELL(IWR,IGIN,ierr,ndim,dolvatm,rscale,rprim)
        CLOSE(unit=IGIN)
        if( ierr.ne.0 )then
          write(IWR,*) '***** ERROR: reading cell info=',ierr
          goto 100
        else
          write(IWR,*) 'Cell info successfully read from tape file'
        endif
        isetuc = 1
c
      if( iopt.eq.1 )then
c       Skip options, go get file of coordinates
        opt = 'c'
        goto 101
      endif
c
c >>>>>> Select options:
c
      mdim = -1
  100 continue
      iopt = 0
      write(IWR,*)
      write(IWR,*)  'Enter option (t,c,h=help,g,n,s,u,z,v,r,x):'
      read(IRD,'(a1)')  opt
c
  101 continue
      if(     opt .eq. 't' )then
c       Go get cell info from tape file
        goto 90
      elseif( opt .eq. 'c' )then
c       Go get file of atomic coordinates
        goto 200
      elseif( opt .eq. 'g' )then
c       Go get geometry file with coordinates
        goto 300
      elseif( opt .eq. 'n' )then
        write(IWR,*) 'Enter ndim:'
        read(IRD,*)  ndim
        mdim = ndim
      elseif( opt .eq. 'l' )then
        dolvatm = .not. dolvatm
        write(IWR,*) 'lattice coordinates toggled to',dolvatm
      elseif( opt .eq. 's' )then
        write(IWR,*) 'Enter pscale:'
        read(IRD,*)  pscale
        do  id=1,3
          rscale(i) = pscale*rscale(i)
        enddo
        isetuc = 1
      elseif( opt .eq. 'u' )then
        if( mdim.lt.0 )then
          write(IWR,*) '**** NEED NDIM SET TO USE THIS'
          goto 100
        endif
        write(IWR,*) 'Enter uscale:'
        read(IRD,*)  uscale
        do  id=1,ndim
          rscale(i) = uscale*rscale(i)
        enddo
        isetuc = 1
      elseif( opt .eq. 'x' )then
        write(IWR,*) 'Enter xscale:'
        read(IRD,*)  xscale
        rscale(1) = xscale*rscale(1)
        isetuc = 1
      elseif( opt .eq. 'z' )then
        write(IWR,*) 'Enter yscale:'
        read(IRD,*)  yscale
        rscale(2) = xscale*rscale(2)
        isetuc = 1
      elseif( opt .eq. 'z' )then
        write(IWR,*) 'Enter zscale:'
        read(IRD,*)  zscale
        rscale(3) = xscale*rscale(3)
        isetuc = 1
      elseif( opt .eq. 'v' )then
        write(IWR,*) 'Enter unit cell primitive vectors:'
        read(IRD,*)  rprim
        isetuc = 1
      elseif( opt .eq. 'r' )then
        write(IWR,'(a,i1)') 'ndim=',ndim
        write(IWR,'(a,3f15.8)')  'rscale=',rscale
        write(IWR,'(a,3f15.8/(6x,3f15.8))') 'rprim=',rprim
      elseif( opt .eq. 'x' )then
        goto 999
      else
        if( opt .ne. 'h' .or. opt .ne. '?' )
     $   write(IWR,*) '***** ERROR: no such option: ',opt
        write(IWR,9001) 
 9001   format(1x,'Available options:',
     $        /3x,'t - get unit cell info from tape file',
     $        /3x,'c - input coordinate file',
     $        /3x,'g - read geometry data',
     $        /3x,'n - set ndim (default=3)',
     $        /3x,'s - set scale factor for problem',
     $        /3x,'u - set uscale factor',
     $        /3x,'x,y,z - set x,y,z scale factor',
     $        /3x,'v - set unit cell primitive vectors',
     $        /3x,'r - recap current cell parameters',
     $        /3x,'h - help: get this message',
     $        /1x,'Standard path:',
     $        /3x,'t (to enter a tape5 for cell info)',
     $        /3x,'c (to enter coordinate list)',
     $        /3x,'You will be prompted for file names'/)
      endif
      goto 100
c
c >>>> Get coordinate file name:
c
  200 continue
      write(IWR,*)
      write(IWR,*) 'Enter number of atoms:'
      read(IRD,*)  natom
      if( natom .gt.maxatm )then
        write(IWR,*) 'ERROR: maximum atoms, maxatm=',maxatm
        STOP
      endif
c
      write(IWR,*) 'Enter name of geometry file:'
      read(IRD,8000)  geomfil
c
c Read in coordinates of atomic system:
      OPEN(IGIN,file=geomfil,form='formatted',status='old',err=1302)
      call RDCOORD( IGIN, natom, ratm0 )
      CLOSE(unit=IGIN)
c
      goto 400
c
c >>>> Get geometry file name:
c
  300 continue
      write(IWR,*)
      write(IWR,*) 'Enter geometry file name:'
      read(IRD,8000)  geomfil
 8000 format(a)
      if( geomfil.eq.' ' ) geomfil = 'atoms.in'
c
      OPEN(IGIN,file=geomfil,form='formatted',status='old',err=1302)
c
c Read in coordinates of atoms in system:
      call GETGEOM( IWR,IGIN,ierr, maxatm,natom,gscale,rcube,ratm0 )
c
      CLOSE(unit=IGIN)
      if( ierr.ne.0 ) goto 1302
      write(IWR,*) 'Geometry successfully read'
c
c Check if need to set lattice stuff:
      if( isetuc.eq.0 )then
c       Default is a cubic cell:
        pscale = gscale
        rprim(1,1) = rcube
        rprim(2,2) = rcube
        rprim(3,3) = rcube
      endif
      write(IWR,*)
c
c >>>> Scale coordinates and lattice vectors:
c
  400 continue
c  Convert into angstrom:
      do  i=1,3
        rscale(i) = Aperau*rscale(i)
        rprim(i,1) = Aperau*rprim(i,1)
        rprim(i,2) = Aperau*rprim(i,2)
        rprim(i,3) = Aperau*rprim(i,3)
      enddo
c
      call RINTRNL( ndim,natom, dolvatm, rprim,orig,rscale,
     $ ratm0,ratm  )
c
  500 continue
      write(IWR,*)
      write(IWR,*) 'Enter bond range cutoff (in Angstrom):'
      read(IRD,*) bond
      if( bond .lt. 0.10 .or. bond .gt. 10.00 ) bond = 3.20
c
  600 continue
      write(IWR,*) 'Enter output file name:'
      read(IRD,8000)  rngfil
      if( rngfil.eq.' ' ) rngfil = 'bonds.out'
c
      OPEN(ILST,file=rngfil,form='formatted',status='new',err=600)
c
      call RNGALL( IWR,ILST, natom, ratm, ndim,rprim, bond )
c
  700 continue
      write(IWR,*) 'Enter atom pair (0 0 to exit):'
      read(IRD,*,err=990)  iatm,jatm
      if( iatm.le.natom .and. jatm.le.natom 
     $  .and. iatm.ge.1 .and. jatm.ge.1 )then
       call RNG2AT( rij, ndim,rprim, ratm(1,iatm),ratm(1,jatm) )
       write(IWR ,9100)  iatm,jatm,rij
       write(ILST,9100)  iatm,jatm,rij
 9100  format(1x,'Atoms (',i5,'-',i5,'), Rij=',f8.4)
       goto 700
      endif
c
  990 continue
      CLOSE(unit=ILST,status='keep')
      write(IWR,*)
  999 continue
      STOP 'all done'
c
 1302 continue
      write(IWR,*)  '***** ERROR: no such file'
      STOP 'tryagain'
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETCELL
c
c
      SUBROUTINE GETCELL(IWR,IDAT,ierr,ndim,dolvatm,rscale,rprim)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER (nline=36)
      CHARACTER*(nline)  label2
      CHARACTER  readlbl*6,wantlbl*6,  label*6
      EQUIVALENCE  (readlbl,label,label2)
      LOGICAL    dolvatm
      DIMENSION  rprim(3,4),rscale(3)
      DIMENSION  strain(3,3),tmp(3,3)
c
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>>>> Find starting point for cell information in input deck:
c
      write(IWR,*) '******* DIGEST INPUT FILE *******'
   10 continue
      ierr = 1
      read(IDAT,8000,end=1300,err=1300) readlbl
      if(     readlbl(1:5) .eq. 'title' )then
        goto 100
      elseif( readlbl.eq. 'dimens' .or. readlbl.eq.'lattic')then
        goto 200
      endif
      write(IWR,9000)  label
      goto 10
c
 8000 format(a6)
 8002 format(a36)
 9000 format(a6)
 9001 format(a)
 9020 format(1x,i5)
 9010 format(1x,3f16.8)
c
c >>>>>> Check for title lines in data deck:
c (optional)
  100 continue
C     read(IDAT,8000)  label
c
      if( readlbl(1:5) .eq. 'title' )then
        write(IWR,9000)  label
        read( readlbl(6:6), '(i1)' )  ntitl
        if( ntitl.eq.0 ) ntitl = 1
        write(IWR,*) 'Title line(s) for this run:', ntitl
        if( ntitl.gt.0 )then
          do  ititl=1,ntitl
            read(IDAT,8000)  label
            write(IWR,9000)  label
          enddo
        endif
        read(IDAT,8000)  label
      else
        BACKSPACE(unit=IDAT)
      endif
c
c >>>>>> PROBLEM DIMENSION (0=cluster thru 3=bulk):
c (required)
  200 continue
      wantlbl = 'dimens'
      write(IWR,9000)  label
      if( readlbl.eq.wantlbl .or. readlbl.eq.'lattic' )then
        read(IDAT,*   )  ndim
        write(IWR,9020)  ndim
        if( ndim.lt.0 .or. ndim.gt.3 ) STOP 'ndim err'
c
        read(IDAT,8000)  label
      else
        read(IDAT,8000)  label
c       Cycle through until we find a dimension statement:
        goto 200
C        write(IWR,*) 'dimension of problem (0=cluster thru 3=bulk)?'
C        ierr = 2
C        goto 1300
      endif
c
c >>>>>> COORDINATES: cartesian or lattice
c (optional)
c  Set default to be Cartesian, not lattice vector, coordinates:
      dolvatm = .false.
c
      if( readlbl .eq. 'coordi' )then
        write(IWR,9000)  label
c
        read(IDAT,8002)  label2
        write(IWR,9001)  label2
c
c       Parse out the type of atomic coordinates
        call STRPARS( label2, 1,nline, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'blank coordinate type' )
        if( n1 .eq. 0 ) STOP 'coordi'
        wantlbl = '      '
        nlbl = MIN( n2-n1+1, 6 )
        wantlbl(1:nlbl) = label2(n1:n2)
c
        if( wantlbl .eq. 'lattic' .or. wantlbl .eq. 'Lattic' .or.
     $      wantlbl .eq. 'LATTIC' )then
          dolvatm = .true.
c
          read(IDAT,8000)  label
c
          if( readlbl.eq.'lvscal' )then
            write(IWR,9000)  label
            read(IDAT,*   )  (rprim(i,4),i=1,3)
            write(IDAT,9010) (rprim(i,4),i=1,3)
            read(IDAT,8000)  label
          else
            rprim(1,4) = one
            rprim(2,4) = one
            rprim(3,4) = one
          endif
          do  i=1,3
            if( rprim(i,4).le.zero ) rprim(i,4) = one
          enddo
c
        elseif( wantlbl .eq. 'cartes' .or. wantlbl .eq. 'Cartes' .or.
     $          wantlbl .eq. 'CARTES' )then
          dolvatm = .false.
c
          read(IDAT,8000)  label
        else
          write(IWR,*)  '***** ERROR: must be lattice or cartesian'
          STOP 'coord'
        endif
c
      else
        write(IWR,9001)  'coordinates given in (lattice,cartesian):'
        write(IWR,9001)  'cartesian'
      endif
c
c >>>>>> SCALING FACTORS
c (optional)
c  set up defaults
      do  id=1,3
        rscale(id) = one
      enddo
      scalep = one
      scaleu = one
      scalex = one
      scaley = one
      scalez = one
c
c  scalep = Global problem scaling factor:
      wantlbl = 'scalep'
      if( readlbl.eq.wantlbl .or. readlbl.eq.'scale ' )then
        write(IWR,9000)  label
        read(IDAT,*   )  scalep
        write(IWR,9010)  scalep
        read(IDAT,8000)  label
      else
        write(IWR,9001) 'scalep: scale factor for problem'
        write(IWR,9010)  scalep
      endif
c
c  scaleu = Global scaling factor, in periodic directions:
c    NB: using "scaleu" (periodic directions) makes major assumption
c    about directions of unit cell vectors, an assumption that is
c    *not* checked in this code.  Use at own risk.
      wantlbl = 'scaleu'
      if( readlbl.eq.wantlbl .or. readlbl.eq.'uc sca' )then
        write(IWR,9000)  label
        read(IDAT,*   )  scaleu
        write(IWR,9010)  scaleu
        read(IDAT,8000)  label
C      else
C        write(IWR,9001) 'scaleu factor: scale periodic directions'
C        write(IWR,9010)  scaleu
      endif
c
c  scalex = Global scaling fector in x-direction:
      wantlbl = 'scalex'
      if( readlbl.eq.wantlbl .or. readlbl.eq.'xscale' )then
        write(IWR,9000)  label
        read(IDAT,*   )  scalex
        write(IWR,9010)  scalex
        read(IDAT,8000)  label
C      else
C        write(IWR,9001) 'scalex factor: scale x-direction'
C        write(IWR,9010)  scalex
      endif
c
c  scaley = Global scaling fector in y-direction:
      wantlbl = 'scaley'
      if( readlbl.eq.wantlbl .or. readlbl.eq.'yscale' )then
        write(IWR,9000)  label
        read(IDAT,*   )  scaley
        write(IWR,9010)  scaley
        read(IDAT,8000)  label
C      else
C        write(IWR,9001) 'scaley factor: scale y-direction'
C        write(IWR,9010)  scaley
      endif
c
c  scalez = Global scaling fector in z-direction:
      wantlbl = 'scalez'
      if( readlbl.eq.wantlbl .or. readlbl.eq.'zscale' )then
        write(IWR,9000)  label
        read(IDAT,*   )  scalez
        write(IWR,9010)  scalez
        read(IDAT,8000)  label
      else
        write(IWR,9001) 'scalez factor: scale z-direction'
        write(IWR,9010)  scalez
      endif
c
c Unify scaling factors:
c
      do  id=1,3
        rscale(id) = scalep*rscale(id)
      enddo
      if( ndim.lt.0 .or. ndim.gt.3 ) STOP 'ndim-sc'
      do  id=1,ndim
        rscale(id) = scaleu*rscale(id)
      enddo
      rscale(1) = scalex*rscale(1)
      rscale(2) = scaley*rscale(2)
      rscale(3) = scalez*rscale(3)
c
c >>>>>> STRAIN (test code)
c (optional)
      istrain = 0
      wantlbl = 'strain'
      if( readlbl.eq.wantlbl )then
        if( .not. dolvatm )then
          write(IWR,*) '***** ERROR: strain only in lattice coords'
          STOP 'str-nolv'
        endif
        istrain = 1
        write(IWR,9000)  label
        do  j=1,3
          read(IDAT,*   )  (strain(i,j),i=1,3)
          write(IWR,9010)  (strain(i,j),i=1,3)
        enddo
c
        read(IDAT,8000)  label
c
        strfac = one
        wantlbl = 'strfac'
        if( readlbl.eq.wantlbl )then
          write(IWR,9000)  label
          read(IDAT,*   )  strfac
          do  j=1,3
            do  i=1,3
              strain(i,j) = strfac*strain(i,j)
            enddo
          enddo
c
          read(IDAT,8000)  label
        else
          write(IWR,*   )  'strfac: scale factor for strain'
        endif
        write(IWR,9010)  strfac
c
      endif
c
c >>>>>> PRIMITIVE LATTICE VECTORS FOR UNIT CELL
c (required)
      wantlbl = 'primit'
      write(IWR,9000)  label
      if( readlbl.eq.wantlbl )then
        do  idim=1,3
          read(IDAT,*   )  (rprim(i,idim), i=1,3)
          write(IWR,9010)  (rprim(i,idim), i=1,3)
          mdim = idim
        enddo
c Scale primitive lattice vectors, periodic directions *ONLY*:
c This needs to be done here since k-vectors will depend on
c rprim() if input for k-vectors done in lattice coordinates.
        if( ndim.lt.0 .or. ndim.gt.3 ) STOP 'ndim-lv'
        do  j=1,ndim
          rprim(j,1) = rscale(j)*rprim(j,1)
          rprim(j,2) = rscale(j)*rprim(j,2)
          rprim(j,3) = rscale(j)*rprim(j,3)
        enddo
        if( istrain.eq.1 )then
c         Install strain into lattice vectors
          do  j=1,3
            tmp(i,j) = strain(1,j)*rprim(i,1)
     $               + strain(2,j)*rprim(i,2)
     $               + strain(3,j)*rprim(i,3)
          enddo
          do  j=1,3
            rprim(i,j) = rprim(i,j) + tmp(i,j)
          enddo
        endif
      else
        write(IWR,*) 'primitive vectors?'
        ierr = 6
        goto 1300
      endif
c
      ierr = 0
c
 1300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETGEOM
c
c
      SUBROUTINE GETGEOM(IWR,IGIN,ierr,maxatm,natom,gscale,rcube,ratm)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  ratm(3,*)
      CHARACTER  wantlbl*6,readlbl*6,line*72
      EQUIVALENCE  (readlbl,line)
c parameters for atC cells:
C      DATA  r64 / 14.3052324648 /,r216 / 10.874762 /
c parameters for Si cells:
      DATA  r64 / 20.40d0 /,r216 / 30.60d0 /
c
 8000 format(a)
      ierr = 0
      read(IGIN,8000)  line
      write(IWR,8000)  line
      if( readlbl .ne. 'iterat' ) goto 1310
      read(IGIN,*) iter
c
      iscale = 0
      scale = 1.d0
      read(IGIN,8000)  line
      write(IWR,8000)  line
      wantlbl = 'number'
      if( readlbl .ne. wantlbl )then
        wantlbl = 'center'
        if( readlbl .ne. wantlbl ) goto 1310
c        I have a QUEST format geometry
        iform = 1
        read(IGIN,8000)  line
        read(line,*,err=20)   natom,lenlbl,rcube
        write(*,*)  'a0 reset to ',rcube
        goto 30
   20   read(line,*)  natom,lenlbl
        if( natom.eq.64 )then
          rcube = r64
        elseif( natom.eq.216 )then
          rcube = r216
        else
          ierr = 3
          write(IWR,*) '***** ERROR: need a lattice parameter'
          RETURN
        endif
   30   continue
        write(IWR,*)  'quest: number of atoms=',natom
        read(IGIN,8000)  line
        wantlbl = 'scale'
        if( readlbl .eq. wantlbl )then
          iscale = 1
          read( line(6:72) ,*) scale
          write(IWR,*) 'scale set to ',scale
        else
          BACKSPACE(unit=IGIN)
          write(IWR,*) 'no scale factor read'
        endif
      else
c       FEIBELMAN format
        iform = 2
        read(IGIN,'(i2)') natom
        lenlbl = 8
        read(IGIN,8000)  line
        write(IWR,*) 'serial: number of atoms=',natom
      endif
c
c get cube dim and scale, in Angstrom
      rcube = rcube/scale
c
      if( natom.gt.maxatm )then
        write(IWR,*) '***** ERROR: too many atoms'
        write(IWR,*) '>>>>> have ',natom,', maxatm=',maxatm
        STOP 'max atom'
      endif
c
      write(*,'(a,2f12.8)') 'rcube,scale=',rcube,scale
      lrd = lenlbl + 1
      do  iatm=1,natom
        read(IGIN,8000)  line
        read(line(lrd:72),*)  (ratm(i,iatm),i=1,3)
      enddo
c
      write(IWR,*) 'Atomic coordinates read successfully'
c
      RETURN
c
 1310 write(IWR,*)  '***** ERROR: looking for label:',wantlbl
      write(IWR,*)  '                         found:',readlbl
      ierr = 1
      RETURN
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDCOORD
c
c
      SUBROUTINE RDCOORD( IDAT, natom, ratm )
c---------------------------------------------------------------
c Purpose: read in coordinates from an unknown format geoemtry file
c Written: Peter A. Schultz, 7-November-2002
c Revision history: none
c---------------------------------------------------------------
c
c The structure of the file is presumed as follows:
c (A) The atom lines might have one of two formats:
c     (1) Quest format: atom_name atom_type x_atom y_atom z_atom
c  or (2) Pure coordinates: x_atom y_atom z_atom
c (B) The geometry might have a heading line, presumed to be
c     either: "atom, type, position" or "coordinates".
c If one of these heading lines is present, it is used to determine
c determine the format of the atom lines.
c If neither of these is the heading, the code guesses that the last
c three fields are the coordinates, and attempts to extract the
c coordinates from there.
c 
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output coordinates
      DIMENSION  ratm(3,natom)
c Local declarations:
      CHARACTER  line*80, label*6
      EQUIVALENCE  (line,label)
      DATA  ncline / 80 /
c
c >>>> EXECUTABLE CODE:
c
   10 continue
c 
      read(IDAT,'(a)')  line
      if( label.eq.'coordi' )then
c       Coordinate format, just: x_atom y_atom y_atom
        iatmfmt = 0
        nfields = 3
        read(IDAT,'(a)') line
      elseif( label .eq. 'atom, ' )then
c       Quest format: atom_name atom_type x_atom y_atom y_atom
        iatmfmt = 2
        nfields = 5
        read(IDAT,'(a)') line
      else
c       This is assumed to be an atom line.
c       Count up total number of fields in first line
c       Use this for all lines, taking last three fields for coords
        nfields = 0
        ncl1 = 1
  100   continue
        call STRPARS( line, ncl1,ncline, ic1,icn )
        if( ic1.eq.0 ) goto 101
        nfields = nfields + 1
        ncl1 = icn + 1
        if( ncl1 .le. ncline ) goto 100
c
  101   continue
        if( nfields .lt. 3 )then
          write(*,'(a)') line
          write(*,'(a)') 'this line cannot be atom line, try next'
          goto 10
        endif
      endif
c
      if( nfields .lt. 3 )then
        write(*,'(a)') 'error reading atom line, not enough fields:'
        write(*,'(a)') line
        STOP 'atomline'
      endif
c
c Parse first atom line:
      nskip = nfields - 3
      iatm = 1
      call RDATMLN( line,ncline, nskip, ratm(1,iatm) )
      do  iatm=2,natom
        read(IDAT,'(a80)') line
        call RDATMLN( line,ncline, nskip, ratm(1,iatm) )
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDATMLN
c
c
      subroutine RDATMLN( line,ncline, nskip, ratm )
c---------------------------------------------------------------
c Purpose: get atom coordinates out of line
c Written: Peter A. Schultz, 7-November-2002
c Revision history: none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
      DIMENSION  ratm(3)
c
c >>>> EXECUTABLE CODE:
c
      icl1 = 1
      ncl1 = 1
      do  iskip=1,nskip+1
c       Find beginning of first field, after skipping nskip fields:
        if( ncl1 .gt. ncline ) STOP 'bad-atmA'
        call STRPARS( line, ncl1,ncline, ic1,icn )
        if( ic1 .eq. 0 ) STOP 'bad-atmB'
        ncl1 = icn + 1
        icl1 = ic1
      enddo
c
c Take the coordinates from the rest of line, free format:
      read( line(icl1:ncline), * )  ratm
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRPARS
c
c
      subroutine STRPARS( line, nl1,nl2, n1,n2 )
c---------------------------------------------------------------
c Purpose: parse bounds of first character string within line
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
c
c local declarations:
      CHARACTER*(1)  blank
      DATA           blank / ' ' /
c
c >>>> EXECUTABLE CODE:
c
      n1 = 0
      n2 = nl2
      do  il=nl1,nl2
        if( line(il:il) .ne. blank )then
c         Found the first non-blank character ...
          n1 = il
          n2 = n1
c          ... and find the last character in string:
          do  jl=n1+1,nl2
            if( line(jl:jl) .eq. blank ) goto 999
            n2 = jl
          enddo
          goto 999
        endif
      enddo
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RINTRNL
c
c
      subroutine RINTRNL( ndim,nvec, dolvatm, rprim,orig,rscale,
     $ rextr,rintr )
c---------------------------------------------------------------
c Purpose: take atom coordinates, shift, scale to/from user input
c
c Written: Peter A. Schultz, 24-September-1997, v2.21
c
c Revision history:
c   7Jul00-PAS/2.45: unify scaling functions; lattice coords
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Note: For lattice coordinates, I *enforce* the assumption,
c       *not* checked, that for reduced dimension systems,
c       non-periodic axes are in special directions.
c       E.g., slab normal in a 2D calculation is the z-axis
c Note: I might consider pulling ruc() up
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input declarations:
      LOGICAL    dolvatm
      DIMENSION  rprim(3,4),orig(3),rscale(3)
c input/output arrays:
      DIMENSION  rextr(3,nvec),rintr(3,nvec)
c local declarations:
      DIMENSION  ruc(3,3), vect(3)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( dolvatm )then
c
c       Use reduced dimension unit cell vectors:
        do  j=1,3
          do  i=1,3
            ruc(i,j) = zero
          enddo
          ruc(j,j) = one
        enddo
        do  j=1,ndim
          do  i=1,ndim
            ruc(i,j) = rprim(i,j) / rprim(j,4)
          enddo
        enddo
c       I choose to install scaling factor here.
        do  i=ndim+1,3
          ruc(i,i) = ruc(i,i) * rscale(i)
        enddo
c
c       Expand the vectors into coordinate space:
        do  ivec=1,nvec
          call LVXPND( ruc, rintr(1,ivec), rextr(1,ivec) )
        enddo
c
      else
c
        do  ivec=1,nvec
          rintr(1,ivec) = rextr(1,ivec)*rscale(1)
          rintr(2,ivec) = rextr(2,ivec)*rscale(2)
          rintr(3,ivec) = rextr(3,ivec)*rscale(3)
        enddo
c
      endif
c
      do  ivec=1,nvec
        rintr(1,ivec) = rintr(1,ivec) - orig(1)
        rintr(2,ivec) = rintr(2,ivec) - orig(2)
        rintr(3,ivec) = rintr(3,ivec) - orig(3)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVXPND
c
c
      subroutine LVXPND( a, r, v )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  a(3,3), r(3)
c output arrays:
      DIMENSION  v(3)
c
c >>>> EXECUTABLE CODE:
c
      r(1) = v(1)*a(1,1) + v(2)*a(1,2) + v(3)*a(1,3)
      r(2) = v(1)*a(2,1) + v(2)*a(2,2) + v(3)*a(2,3)
      r(3) = v(1)*a(3,1) + v(2)*a(3,2) + v(3)*a(3,3)
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RNGALL
c
c
      SUBROUTINE RNGALL( IWR,ILST, natom, ratm,ndim,rprim, bond )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (maxnn=15)
      DIMENSION  ratm(3,natom),rprim(3,*)
      DIMENSION  nfold(maxnn), jnn(maxnn),rnn(maxnn)
c
      n0fold = 0
      do  i=1,maxnn
        nfold(i) = 0
      enddo
      bndrms = 0.d0
      bnderr = 0.d0
c
      write(ILST,*)  ' Atom  #nn  [Inn  Rnn(Angstrom)]'
c
      do 200 iatm=1,natom
        nn = 0
        do 100 jatm=1,natom
          if( jatm.eq.iatm ) goto 100
c
          call RNG2AT( rng, ndim,rprim, ratm(1,iatm),ratm(1,jatm) )
          if( rng .le. bond )then
            nn = nn + 1
            if( nn.le.maxnn )then
c             Record this atom
              jnn(nn) = jatm
              rnn(nn) = rng
            else
c             Keep maxnn nearest atoms:
              imx = 1
              rmx = rnn(1)
              do  mm=2,maxnn
                if( rnn(mm).gt.rmx )then
                  imx = mm
                  rmx = rnn(mm)
                endif
              enddo
              do  mm=imx,maxnn-1
                jnn(mm) = jnn(mm+1)
                rnn(mm) = rnn(mm+1)
              enddo
              jnn(maxnn) = jatm
              rnn(maxnn) = rng
            endif
          endif
  100   continue
        nntot = nn
        nn = MIN( nntot, maxnn )
        bndord = 0.d0
C        do  i=1,nn
C          if( rnn(i) .lt. 1.75 )then
C            bndord = bndord + (1.75-rnn(i))/.20
C          endif
C        enddo
        if( nn.gt.0 )then
          nfold(nn) = nfold(nn) + 1
        else
          n0fold = n0fold + 1
        endif
        if( nn.le.5 )then
          write(ILST,9106)  iatm,nntot,(jnn(i),rnn(i),i=1,nn)
        else
          write(ILST,9107)  iatm,nntot,(jnn(i),rnn(i),i=1,nn)
        endif
        bnderr = bnderr + ABS(bndord-4.d0)
        bndrms = bndrms + (bndord-4.d0)**2
 9106   format(1x,i3,i5, 5(i6,f8.4) )
 9107   format(1x,i3,i5, 5(i6,f8.4) / (9x,5(i6,f8.4)) )
  200 continue
c
      if( n0fold.gt.0 )then
          write(IWR ,'(6(i3,"-fold atoms:",i4))') 0,n0fold
          write(ILST,'(6(i3,"-fold atoms:",i4))') 0,n0fold
      endif
      do  i=1,maxnn
        if( nfold(i).ne.0 )then
          write(IWR ,'(6(i3,"-fold atoms:",i4))') i,nfold(i)
          write(ILST,'(6(i3,"-fold atoms:",i4))') i,nfold(i)
        endif
      enddo
c
      bnderr = bnderr/DBLE(natom)
      bndrms = SQRT(bndrms/DBLE(natom))
      write(ILST,*) 'bond err, rms err=',bnderr,bndrms
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RNG2AT
c
c
      SUBROUTINE RNG2AT( rij, ndim,rprim, ri,rj )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  rprim(3,3)
      DIMENSION  ri(3),rj(3)
      DIMENSION  rt(3), ru(3)
      DATA  half / 0.5d0 /
c
      do  i=1,3
        rt(i) = ri(i) - rj(i)
      enddo
      rij = SQRT( rt(1)**2 + rt(2)**2 + rt(3)**2 )
c
c Set up which lattice vectors to search over:
      na = 0
      nb = 0
      nc = 0
      if( ndim.gt.0 )then
        na = 1
        if( ndim.gt.1 )then
          nb = 1
          if( ndim.eq.3 )then
            nc = 1
          endif
        endif
      endif
c
c Search for closest image:
      do  ic=-nc,nc
        do  ib=-nb,nb
          do  ia=-na,na
            do  i=1,3
              ru(i) = ia*rprim(i,1) + ib*rprim(i,2) + ic*rprim(i,3)
              ru(i) = ru(i) + rt(i)
            enddo
            rijx = SQRT( ru(1)**2 + ru(2)**2 + ru(3)**2 )
            if( rijx.lt.rij ) rij = rijx
          enddo
        enddo
      enddo
c
      RETURN
      END
