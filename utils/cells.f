      PROGRAM CELLS
c
c###############################################################
c
c Program: CELLS
c
c Written: Peter A. Schultz, January 1999, March 2014
c
c Purpose: generate NxMxL larger unit cells from basis cell
c          Secondary: convert from LAMMPS and VASP formats
c
c Revision history:
c  26Mar14-PAS/1.50: CELL, VASP poscar converter installed
c  20Nov12-PAS/1.43: LAMMPS-compat: turn off numeric order req.
c  21Jun06-PAS/1.4: set up for larger dimensions, naming opts.
c   3Aug04-PAS/1.3: bias distance sort, add greeting
c  11Jan04-PAS/1.2: offsets
c  19Nov03-PAS/1.1: general Quest geom reader
c  15Sep03-PAS/1.0a: fix bug in rprim expansion (thx A. Thompson)
c  18Nov02-PAS/1.0: add sorting options, input general atomline
c  19Jul00-PAS: add skewed cells, origins
c  21Apr00-PAS: add dimensionality
c###############################################################
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c natmd=maximum number of total atoms
      PARAMETER  ( natmd = 5000 )
c ntypd=maximum count of types of atoms
      PARAMETER  ( ntypd = 10 )
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'1.51 '/, lastup/'29Apr15-PAS'/
c
      DIMENSION  rprim(3,3), rprimb(3,3)
      DIMENSION  ratmuc(3,natmd),itypuc(natmd)
      DIMENSION  ratm(3,natmd),itypa(natmd)
      DIMENSION  rsqat(natmd), rsqbat(natmd), rbias(3)
      DIMENSION  r1(3), nlat(3), rlat(3,27)
      DIMENSION  rscale(3)
      DIMENSION  orig(3),offset(3),rvec(3),rvecb(3)
      DIMENSION  nlsc(3,3)
      CHARACTER  filnm*80, line*80, sorttyp*1
      CHARACTER  opt*1
      LOGICAL    incell, lstcell, lvatm
c
c natmnm=string length for atom/type names (from Quest)
      PARAMETER  ( natmnm = 12 )
      CHARACTER*(24)  typnm(ntypd), atmnm(natmd), atmnmt
      CHARACTER*(24)  atmnmuc(natmd)
      CHARACTER  attyp*2,atnum*4,atname*6
c This determines output atom name stylet
      DATA       iatnames / 2 /
c
      DATA  IRD,IWR, IGFL, IDT / 5,6, 7,8 /
      DATA  zero,one / 0.d0,1.d0 /
      DATA  rbias / 0.000035d0, 0.000027d0, 0.000021d0 /
      DATA  A2au / 0.5291771d0 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,9000)  version, lastup
 9000 format(
     $ 1x,'Program CELLS: generates Quest supercells'
     $/1x,'Version:  ',a,'     Last revised: ',a
     $ /)
c
c Set up some defaults ...
c
c Set coordinate type and system dimension to "undefined"
      icoordtyp = -1
      ndim = -1
c Initialize origin, offset to zero's, scaling factors to one
      ioffset = 0
      do  i=1,3
        orig(i) = zero
        offset(i) = zero
        rscale(i) = one
      enddo
c Set output of atom ranges to OFF
      lstcell = .false.
c Set lattice/cartesian coordinates to CARTESIAN
      lvatm = .false.
c
c How do we default order atoms? 1=block,2=interleave,3=R(atm) ...
      iatmorder = 3
c  ... and do we project them into the principal unit cell?
      incell = .false.
c
  100 continue
      write(IWR,*) 'Enter input type: (Geom/Cell/Vasp/Quest):'
      read(IRD,'(a1)')  opt
 8010 format(1x,'Configuration input options:',
     $      /3x,'C,c = CELL format file',
     $      /3x,'G,g = {G,g}eom file input',
     $      /3x,'V,v = {V,v}asp {P,p}oscar format file',
     $      /3x,'  ... and the following not ready yet:',
     $      /3x,'Q,q = {Q,q}uest {T,t}ape5 data file',
     $     //3x,'Enter (single character) option:')
c
      write(IWR,*)
      if( opt.eq. 'C' .or. opt.eq. 'c' )then
c       Read from "CELL" format
        ndim = 3
        write(IWR,*)  'Enter CELL file name'
        read(IRD,'(a)')  filnm
        OPEN(unit=IGFL,file=filnm,status='old',form='formatted',err=100)
        write(IWR,*) 'CELL file opened, processing ...'
        goto 200
c
      elseif( opt.eq. 'V' .or. opt.eq. 'v' .or.
     $        opt.eq. 'P' .or. opt.eq. 'p' )then
c       Read from VASP poscar file
        ndim = 3
        write(IWR,*)  'Enter VASP poscar file name:'
        read(IRD,'(a)')  filnm
        OPEN(unit=IGFL,file=filnm,status='old',form='formatted',err=100)
        write(IWR,*) 'VASP POSCAR file opened, processing ...'
        goto 300
c
      elseif( opt.eq. 'Q' .or. opt.eq. 'q' .or.
     $        opt.eq. 'T' .or. opt.eq. 't' )then
c       Read from Quest tape5 input file
        write(IWR,*)  'Enter Quest input file name:'
        read(IRD,'(a)')  filnm
        OPEN(unit=IGFL,file=filnm,status='old',form='formatted',err=100)
        write(IWR,*) 'Quest input file opened, processing ...'
        write(IWR,*) 'ERROR/DEV: not implemented yet'
        CLOSE(unit=IGFL)
        goto 100
C        goto 400
c
      elseif( opt.eq. 'G' .or. opt.eq. 'g' )then
c       Read atomic positions from Quest geometry/configuration file
        write(IWR,*)  'Get atomic data from geometry file ...'
        goto 500
c
      else
c       Unknown option, tell people what is available:
        write(IWR,8010)
        goto 100
      endif
c
c ***** CELL format input
c
  200 continue
      call RDCELL( ierr, IGFL, rprim, iatmfmt, natmnm,
     $ ntyp, typnm,   natmuc, atmnmuc, itypuc, ratmuc )
      CLOSE( unit=IGFL )
c
      if( ierr.ne.0 )then
        write(IWR,*) 'ERROR processing cell file, ierr=',ierr
        goto 100
      endif
      goto 1000
c
c ***** VASP poscar input
c
  300 continue
      call RDPOSCAR( ierr, IWR, IRD, IGFL,
     $ rprim, iatmfmt,natmnm, lvatm, rscale,
     $ ntyp,ntypd, typnm, natmuc,natmd, atmnmuc,
     $ itypuc, ratmuc )
      CLOSE( unit=IGFL )
c
      icoordtyp = 0
      if( lvatm ) icoordtyp = 3
      write(IWR,*) 'Atom output order set to BLOCK ...'
      iatmorder = 1
      if( ierr.ne.0 )then
        write(IWR,*) 'ERROR processing poscar file, ierr=',ierr
        goto 100
      endif
c
c Install Angstrom-to-bohr scaling factor into *scales*
      do  i=1,3
        rscale(i) = rscale(i) / A2au
      enddo
      goto 1000
c
c ***** Quest input file
c
  400 continue
C      call GETQIN(
C      - read until "number"
C      read(IGFL,*) natm
C      goto 600
      write(IWR,*) '***** QUEST TAPE5 NOT YET IMPLEMENTED *****'
      CLOSE(unit=IGFL)
      goto 100
cdev      goto 1000
c
c ***** interacive user input
c
  500 continue
      write(IWR,*) 'Enter dimensionality:'
      read(IRD,*)  ndim
c
      write(IWR,*) 'Enter lattice vectors of basis cell:'
      read(IRD,*)  rprim
c
      write(IWR,*) 'Enter number of atoms in basis cell:'
      read(IRD,*)  natmuc
      if( natmuc.eq.0 ) goto 999
c
c Ask user for geometry file name:
      IRDR = IRD
      write(IWR,*) 'Enter geometry file name (i=direct input):'
      read(IRD,'(a)') filnm
      OPEN(unit=IGFL,file=filnm,status='old',form='formatted',err=550)
      IRDR = IGFL
c
c     Determine format of file from lead line:
      read(IRDR,'(a80)') line
      write(IWR,'(a80)') line
      call ATFMTIN( iatmfmt, line )
      write(IWR,*) 'format style=',iatmfmt
c     Tell it to pick up atom types without checking if have it:
      if( iatmfmt.eq.2 ) iatmfmt = 3
      goto 600
c
  550 continue
      write(IWR,*) 'Enter {atom type x y z} for atoms in basis cell:'
      iatmfmt = 3
c
  600 continue
c Use Quest configuration reader to extract data:
c  turn off write:
      IWRX = -1
      call CONFIGRD( IRDR,IWRX, iatmfmt, natmnm, ntyp,typnm,
     $ natmuc, atmnmuc, itypuc, ratmuc )
      if( iatmfmt .eq. 3 ) iatmfmt = 2
C     do  iatm=1,natmuc
C       read(IRDR,*)  jatm,itypuc(iatm),(ratmuc(i,iatm),i=1,3)
C     enddo
c
      if( IRDR.eq.IGFL ) CLOSE(unit=IRDR)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Setup options
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1000 continue
c Regurgitate defaults
      if( iatmorder.eq.1 )then
        write(IWR,*) 'Default atom order is by BLOCK'
      elseif( iatmorder.eq.2 )then
        write(IWR,*) 'Default atom order is INTERLEAVE'
      elseif( iatmorder.eq.3 )then
        write(IWR,*) 'Default atom order is distance to origin'
      else
        iatmorder = 3
      endif
      if( incell )then
        write(IWR,*) 'Projection into primary cell is ON/'
      else
        write(IWR,*) 'Projection into primary cell is OFF/'
      endif
c
 1100 continue
      write(IWR,*)
      write(IWR,*) 'Enter option: (g-o,s-ort,o-ffset,h-elp,c-enter):'
      read(IRD,'(a1)')  opt
c
      if(     opt .eq. 'g' .or. opt .eq. 'G' )then
        goto 2000
c
      elseif( opt .eq. 'p' .or. opt .eq. 'P' )then
c       Toggle projection into cell
        incell = .not. incell
        if( incell )then
          write(IWR,*) 'Projection into cell is toggled ON/',incell
        else
          write(IWR,*) 'Projection into cell is toggled OFF/',incell
        endif
c
      elseif( opt .eq. 's' .or. opt .eq. 'S' )then
c       Sort atoms using what protocol?
        write(IWR,*) 'Atom order (1=block,2=interleave,3=R(atm):'
ccc        read(IRD,*) iopt
        read(IRD,'(a80)') line
        do  i=1,80
          if( line(i:i).ne.' ' )then
            sorttyp = line(i:i)
            goto 1121
          endif
        enddo
        sorttyp = ' '
 1121   continue
        if(     sorttyp .eq. 'b' .or. sorttyp .eq. '1' )then
          iopt =1
        elseif( sorttyp .eq. 'i' .or. sorttyp .eq. '2' )then
          iopt =2
        elseif( sorttyp .eq. 'r' .or. sorttyp .eq. '3' )then
          iopt =3
        else
          write(IWR,*) 'Sort type >',sorttyp,'< invalid '//
     $     'revert back to options'
          goto 1100
        endif
        if( iopt.gt.0 .and. iopt.le.3 ) iatmorder = iopt
c
      elseif( opt .eq. 'c' .or. opt .eq. 'C' )then
        write(IWR,*) 'Enter center for R-sort:'
        read(IRD,*)  orig
c
      elseif( opt .eq. 'o' .or. opt .eq. 'O' )then
        write(IWR,*) 'Enter position offset:'
        read(IRD,*)  offset
        ioffset = 1
c
      elseif( opt .eq. 'f' .or. opt .eq. 'F' )then
        write(IWR,'(a,3f16.8)') 'Current scaling factors=',rscale
 1140   continue
        write(IWR,*) 'Change or apply factors (c,a)?'
        read(IRD,'(a1)')  opt
        if( opt.eq. 'c' .or. opt.eq. 'C' )then
          write(IWR,*) 'Enter new scaling factors (3:x,y,&z):'
          read(IRD,*)  rscale
        elseif( opt.eq. 'a' .or. opt.eq. 'A' )then
          write(IWR,*) 'Apply scaling, lvatm,dim=',lvatm,ndim
          call SCALECELL( ndim, lvatm, rscale,
     $                    natmuc, rprim, ratmuc, orig,offset )
c         Applied, so remove from scaling:
          do  i=1,3
            rscale(i) = one
          enddo
        endif
c
      elseif( opt .eq. 'l' .or. opt .eq. 'L' )then
c       Modify lattice vector setting and dimenstionality
        write(IWR,*) 'Currently, lattice units =',lvatm, ', Dim=',ndim
 1150   continue
        if( ndim.lt.0 .or. ndim.gt.3 )then
          write(IWR,*) 'Enter current problem dimension (e.g., 2=slab):'
          read(IRD,*)  ndim
          goto 1150
        endif
c
 1155   continue
        write(IWR,*) 'Lattice units=', lvatm,
     $          ' --- Convert to ',.not.lvatm,' (y/n/reset)?'
        read(IRD, '(a1)' )  opt
      
        if( opt.eq. 'r' .or. opt.eq. 'R' )then
c         Relabel lattice unit convention
          write(IWR,*)
          write(IWR,*)'Reset (NOT convert) lattice unit setting(ynlcd)?'
          read(IRD,'(a1)')  opt
          if( opt.eq. 'y' .or. opt.eq. 'Y' )then
            lvatm = .not. lvatm
            write(IWR,*) 'Lattice units toggled to ', lvatm
          elseif( opt.eq. 'l' .or. opt.eq. 'L' )then
            lvatm = .true.
            write(IWR,*) 'Lattice units changed to LATTICE=', lvatm
          elseif( opt.eq. 'c' .or. opt.eq. 'C' )then
            lvatm = .false.
            write(IWR,*) 'Lattice units changed to CARTESIAN=', lvatm
          elseif( opt.eq. 'd' .or. opt.eq. 'D' )then
c           Reset the problem dimension
            ndim = -1
            goto 1150
          endif
        elseif( opt.eq. 'y' .or. opt.eq. 'Y' )then
          write(IWR,*) 'Converting coords to lvatm=', .not.lvatm
          if( lvatm )then
c           Expand into CARTESIAN
            do  iatm=1,natmuc
              call LVXPND( rprim, rvec, ratmuc(1,iatm) )
              do  i=1,3
                ratmuc(i,iatm) = rvec(i)
              enddo
            enddo
            lvatm = .false.
          else
c           Porject onto LATTICE
            do  iatm=1,natmuc
              call LVPROJ( rprim, ratmuc(1,iatm), rvec   )
              do  i=1,3
                ratmuc(i,iatm) = rvec(i)
              enddo
            enddo
            lvatm = .true.
          endif
        endif
        if( lvatm )then
          write(IWR,*) 'Coordinates: LATTICE  , dimension=',ndim
        else
          write(IWR,*) 'Coordinates: CARTESIAN, dimension=',ndim
        endif
c
      elseif( opt .eq. 'r' .or. opt .eq. 'R' )then
c       Toggle output of atom ranges ...
        lstcell = .not. lstcell
c
      elseif( opt .eq. 'n' .or. opt .eq. 'N' )then
c       Set how the names of atoms are done
        write(IWR,*) 'Enter naming option(0=none,1=nums,2=lbls):'
        read(IRD,*) iatnms
        if( iatnms.ge.0 .and. iatnms.le.2 ) iatnames = iatnms
c
      elseif( opt .eq. 'd' .or. opt .eq. 'D' )then
c       Display the cell vectors and some atom coords:
        write(IWR,*) 'Primitive lattice parameters:'
        do  idim=1,3
          write(IWR,'(3x,3f16.8)')  (rprim(i,idim),i=1,3)
        enddo
        write(IWR,*) 'Number of atoms=',natmuc
        write(IWR,*) 'Atom, type, position for lead atoms:'
        natmwr = MIN( natmuc, 5 )
        do  iatm=1,natmwr
          ityp = itypuc(iatm)
          write(IWR,9040) atmnmuc(iatm),typnm(ityp),
     $                    (ratmuc(i,iatm),i=1,3)
        enddo
 9040   format( 1x, a, 3x, a, 3x, 3f16.8 )
c
      else
        write(IWR,9050)
 9050   format(/1x,'Available options:'
     $   /5x,'h(elp)   = this message'
     $   /5x,'o(ffset) = offset coordinate origin'
     $   /5x,'c(enter) = set alternate center for R-sort'
     $   /5x,'l(attice)= lattice.cartesian units, dimension'
     $   /5x,'p(roj)   = project atom coords into primary cell?'
     $   /5x,'f(actors)= scaling *F*actors'
     $   /5x,'n(aming) = atom naming style'
     $   /5x,'d(isplay)= display latvecs and lead coords'
     $   /5x,'s(ort)   = specify sort of atom list -'
     $   /5x,'         (1) 1-N,1a-Na,...  [blocked by uc]'
     $   /5x,'         (2) 1,1a,...,2,2a,... [interleaved]'
     $   /5x,'         (3) R1<R2<R3... [distance from origin]'
     $   /)
      endif
      goto 1100
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Completed setup, now go on to generating output
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 2000 continue
c  First, install all scaling functions.
c
      call SCALECELL( ndim, lvatm, rscale,
     $                natmuc, rprim, ratmuc, orig,offset )
        natmwr = MIN( natmuc, 5 )
        do  iatm=1,natmwr
          ityp = itypuc(iatm)
          write(IWR,9040) atmnmuc(iatm),typnm(ityp),
     $                    (ratmuc(i,iatm),i=1,3)
        enddo
c     Not coming back here, so skip reset rscale=1
c
 2010 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Select expanded cell to generate from primitive cell
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write(IWR,*)
      write(IWR,*) 'Enter nla,nlb,nlc (dims of big cell):'
      read(IRD,*)  nla,nlb,nlc
c
c Set non-periodic directions to single repeat
      if( ndim.lt.3 ) nlc = 1
      if( ndim.lt.2 ) nlb = 1
      if( ndim.lt.1 ) STOP 'ndim=0'
      write(IWR,'(a,3i4)') 'nla,nlb,nlc=',nla,nlb,nlc
      if( nla .le. 0 )then
        write(IWR,*) '>>>> nla.le.0 - all done'
        goto 999
      endif
c
c
      nabig = nla*nlb*nlc*natmuc
      if( nabig.gt.natmd )then
        write(IWR,*) '***** ERROR: generates too many atoms'
        write(IWR,*) '>>>>> # atoms, max=',nabig,natmd
        goto 2010
      endif
c
      if( ilatx.eq.1 )then
        write(IWR,*) 'Enter supercell vectors, lv units:'
        write(IWR,*) '***** BE AFRAID, BE VERY AFRAID *****'
        read(IRD,*) nlsc
        write(IWR,*) '>>>> DISABLED - diagonal expansion ONLY'
      else
        do  j=1,3
          do  i=1,3
            nlsc(i,j) = 0
          enddo
        enddo
        nlsc(1,1) = nla
        nlsc(2,2) = nlb
        nlsc(3,3) = nlc
      endif
c
 2100 continue
      write(IWR,*) 'Enter output file name (or /quit/):'
      read(IRD,'(a)')  filnm
c
      if( filnm(1:4).eq.'quit' ) goto 999
c
c Perform expansion based on input lattice expansion matrix
c  (error in original code discovered by A. Thompson
C      do  idim=1,3
C        do  j=1,3
C          rprimb(j,idim) = DBLE( nlsc(1,idim) )*rprim(j,1)
C     $                   + DBLE( nlsc(2,idim) )*rprim(j,2)
C     $                   + DBLE( nlsc(3,idim) )*rprim(j,3)
C        enddo
C      enddo
c The atom-coord code currently cannot use general form, so ...
c I will do a strictly diagonal expansion A'=lA,B'=mB,C'=nC
      do  idim=1,3
        do  j=1,3
          rprimb(j,idim) = DBLE( nlsc(idim,idim) ) * rprim(j,idim)
        enddo
      enddo
c
      nlat(1) = nla
      nlat(2) = nlb
      nlat(3) = nlc
c
c Copy primary cell over directly:
c
c If interleaved order specified, change stride:
      jdel = 1
      if( iatmorder.eq.2 ) jdel = nla*nlb*nlc
c
      jatm = 1
      do  iatm=1,natmuc
        atmnm(jatm) = atmnmuc(iatm)
        itypa(jatm) = itypuc(iatm)
        do  i=1,3
          ratm(i,jatm) = ratmuc(i,iatm)
        enddo
        jatm = jatm + jdel
      enddo
c
      natm = natmuc
      ncopys = 1
      jcell = 1
      do 2300 idim=1,3
        jatm = natm
        r1(1) = zero
        r1(2) = zero
        r1(3) = zero
        ncopy = ncopys
        do 2250 ilat=2,nlat(idim)
          r1(1) = r1(1) + rprim(1,idim)
          r1(2) = r1(2) + rprim(2,idim)
          r1(3) = r1(3) + rprim(3,idim)
c
          iatm = 1
          jatm = natm + 1
          do 2210 icell=1,ncopy
            jcell = jcell + 1
            if( iatmorder.eq.2 )then
c             Interleaved, start at copy number
              iatm = icell
              jatm = jcell
            endif
            do 2200 iatmuc=1,natmuc
              atmnm(jatm) = atmnm(iatm)
              itypa(jatm) = itypa(iatm)
              ratm(1,jatm) = ratm(1,iatm) + r1(1)
              ratm(2,jatm) = ratm(2,iatm) + r1(2)
              ratm(3,jatm) = ratm(3,iatm) + r1(3)
              jatm = jatm + jdel
              iatm = iatm + jdel
 2200       continue
            natm = natm + natmuc
            ncopys = ncopys + 1
 2210     continue 
c
 2250   continue
c
 2300 continue
c
c Have raw cell, now install desired offset
c
      if( ioffset.ne.0 )then
        write(IWR,*) 'Install offset=',offset
        do  iatm=1,natm
          ratm(1,iatm) = ratm(1,iatm) - offset(1)
          ratm(2,iatm) = ratm(2,iatm) - offset(2)
          ratm(3,iatm) = ratm(3,iatm) - offset(3)
        enddo
      endif
c
      if( iatmorder.eq.3 )then
c       We will be sorting atoms by distance from origin
c        so set up local lattice (search) vectors
c
c       Get nearby lattice vectors ...
c        ... first, restrict dimensionality
        jlc = 1
        jlb = 1
        jla = 1
        if( ndim.lt.3 )then
          jlc = 0
          if( ndim.lt.2 ) jlb = 0
        endif
c        ... and then create lattice vectors:
        ils = 1
        do 2303 ilc=-jlc,jlc
          do 2302 ilb=-jlb,jlb
            do 2301 ila=-jla,jla
c             Put origin at beginning of list
              if( ila.eq.0 .and. ilb.eq.0 .and. ilc.eq.0 )then
                il = 1
              else
                ils = ils + 1
                il = ils
              endif
              do i=1,3
                rlat(i,il) = ila*rprimb(i,1) + ilb*rprimb(i,2)
     $                     + ilc*rprimb(i,3)
              enddo
 2301       continue
 2302     continue
 2303   continue
      endif
c
c  Sort by distance from origin:
c
      do 2400 iatm=1,natm
c
        if( incell .or. iatmorder.eq.3 )then
c
c         Project atom position onto cell vectors:
C         call PROJLV( ratm(1,iatm), rprimb(1,1),rprimb(1,2),rprimb(1,3), r1 )
          call LVPROJ( rprimb, ratm(1,iatm), r1 )
          edge = 0.0001d0
c
c         Put atom position inside primary cell:
          do 2310 idim=1,ndim
            vchk = r1(idim) + edge
            ilat = ABS( vchk )
            if( r1(idim)+edge .lt. zero ) ilat = -ilat - 1
c
            if( ilat.ne.0 )then
              clat = DBLE( ilat )
              do  i=1,3
                ratm(i,iatm) = ratm(i,iatm) - clat*rprimb(i,idim)
              enddo
            endif
 2310     continue
        endif
c
        if( iatmorder.eq.3 )then
c         Sort atoms in big cell by distance to origin
c
c         Find lattice image of atom iatm closest to origin:
c          (biased by small shift to get canonical order)
          do  i=1,3
            rvec(i) = ratm(i,iatm) - orig(i)
            rvecb(i) = ratm(i,iatm) - orig(i) - rbias(i)
          enddo
          rsqmn = zero
          rsqbmn = zero
          do  id=1,ndim
            rsqmn = rsqmn + rvec(id)**2
            rsqbmn = rsqbmn + rvecb(id)**2
          enddo
          ilmn = 1
          do  il=2,ils
            rsq = zero
            rsqb = zero
            do  id=1,ndim
              rsq = rsq + ( rvec(id) - rlat(id,il) )**2
              rsqb = rsqb + ( rvecb(id) - rlat(id,il) )**2
            enddo
            if( rsqb.lt.rsqbmn )then
              rsqmn = rsq
              rsqbmn = rsqb
              ilmn = il
            endif
          enddo
c
c         Shift atom coordinates to image closest to origin:
          ityp = itypa(iatm)
          do  i=1,3
            r1(i) = ratm(i,iatm) - rlat(i,ilmn)
          enddo
c
c         Insert atom iatm into list ordered by (biased) range:
          katm = iatm
          atmnmt = atmnm(iatm)
          do 2350 jatm=iatm,2,-1
            if( rsqbmn.lt.rsqbat(jatm-1) )then
              atmnm(jatm) = atmnm(jatm-1)
              itypa(jatm) = itypa(jatm-1)
              ratm(1,jatm) = ratm(1,jatm-1)
              ratm(2,jatm) = ratm(2,jatm-1)
              ratm(3,jatm) = ratm(3,jatm-1)
              rsqat(jatm) = rsqat(jatm-1)
              rsqbat(jatm) = rsqbat(jatm-1)
              katm = jatm - 1
            else
              goto 2351
            endif
 2350     continue
 2351     continue
          atmnm(katm) = atmnmt
          itypa(katm) = ityp
          ratm(1,katm) = r1(1)
          ratm(2,katm) = r1(2)
          ratm(3,katm) = r1(3)
          rsqat(katm) = rsqmn
          rsqbat(katm) = rsqbmn
c
c         End sorting this atom by distance to origin
        endif
c
c       end loop iatm=1,natm
 2400 continue
c
      if( iatmorder .eq. 3 .and. lstcell )then
c       Atoms sorted by distance to origin, tell user
        do  iatm=1,natm
          rsqat(iatm) = SQRT( rsqat(iatm) )
        enddo
        write(IWR,917)
        write(IWR,*) '>>>> R for range-sorted atoms:'
  915   format(a,5(i6,f8.3))
  916   format(5x,5(i7,f8.3))
        ncol = 5
        nline = ( natm + ncol - 1 ) / ncol
        nrow = nline
        do  iline=1,nline
          iatm1 = iline
C         write(IWR,915) 'range:', ( iatm,rsqat(iatm),
          if( iline.le.nrow .or. iline.eq.nline )then
            if( iline.eq.nline .and.
     $          nline.gt.(nrow+1) ) write(IWR,*) ' ...'
            write(IWR,916)  ( iatm,rsqat(iatm),
     $        iatm=iatm1,natm,nline )
          endif
        enddo
        write(IWR,917)
  917   format(1x, 72('-') )
      endif
c
  910 format(3f12.8)
  920 format(i5)
  930 format(2x,i5,4x,i2,4x,4f12.8)
c
c Install scaling:
C      do  j=1,ndim
C        do  i=1,ndim
C          rprimb(i,j) = rprimb(i,j) / rscale(i)
C        enddo
C      enddo
C      do  iatm=1,natm
C        do  i=1,3
C          ratm(i,iatm) = ratm(i,iatm) / rscale(i)
C        enddo
C      enddo
c
      OPEN(unit=IDT,file=filnm,status='unknown',form='formatted')
      write(IDT,'(a,3i3)') 'primitive lattice vectors',nlat
      write(IDT,910)  rprimb
      write(IDT,'(a)') 'number of atoms in unit cell'
      write(IDT,920)  natm
      write(IDT,'(a)') 'atom, type, position vector'
C     do iatm=1,natm
C       write(IDT,930) iatm,itypa(iatm),(ratm(i,iatm),i=1,3)
C     enddo
c
c renumber/label the atoms:
C        write(IWR,*) 'Labeling the atoms ...'
        attyp = 'AT'
        iatmnum = 0
        do  iatm=1,natm
          iatmnum = iatmnum + 1
          atmnm(iatm) = ' '
          if( iatnames.eq.1 )then
c           Number the atoms and nothing else
            write(atmnm(iatm),'(i6)') iatm
          elseif( iatnames.eq.2 )then
c           Use a label with a number 2(str)+4(num)=name
c           Note: atnum must be equal to chars needed for maxnum
            nstrnum = 4
            maxnum = natmd
            call STRNUM( nstrnum,atnum,iatmnum,maxnum )
            atname = attyp//atnum
            write(atmnm(iatm),'(a6)') atname
c         else
c           The default is to leave the name blank
          endif
        enddo
      write(IWR,*) 'Number of atoms = ',iatmnum
c
C      write(IWR,*) 'Outputing the atomic configuration ...'
      call CONFIGWR( IDT,      iatmfmt, natmnm, ntyp,typnm,
     $ natm  , atmnm  , itypa , ratm   )
      CLOSE(unit=IDT)
c
      goto 2010
c
  999 continue
      STOP 'bye guy'
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDCELL
c
c
      subroutine RDCELL( ierr, IGFL, rprim, iatmfmt, natmnm,
     $ ntyp, typnm,   natm  , atmnm  , itypa , ratm   )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/output declarations:
      INTEGER        IGFL, iatmfmt, natmnm, ntyp,natm
      INTEGER        itypa(*)
      CHARACTER*(*)  typnm(*)
      CHARACTER*(*)  atmnm(*)
      
      DIMENSION      rprim(3,3), ratm(3,*)
c Local declarations
      CHARACTER*(6)  label
c
c >>>> EXECUTABLE CODE
c
      ierr = 1
      write(*,*) 'RDCELL:, IGFL=',IGFL
c
 8010 format(a6)
      read(IGFL,8010)  label
      if( label .ne. 'primit' ) goto 990
      do  lv=1,3
        read(IGFL,*)  ( rprim(i,lv), i=1,3 )
      enddo
c
      read(IGFL,8010)  label
      if( label .ne. 'number' ) goto 990
      read(IGFL,*)  natm
      write(*,*) 'Number of atoms: ',natm
c
c
      read(IGFL,8010)  label
      if( label .ne. 'atom, ' ) goto 990
c     Tell it to pick up atom types without checking if have it:
      iatmfmt = 3
c
c     And do not output anything as it goes, either ...
      IWRX = -1
      call CONFIGRD( IGFL,IWRX, iatmfmt, natmnm, ntyp,typnm,
     $ natm  , atmnm  , itypa , ratm   )
      if( iatmfmt .eq. 3 ) iatmfmt = 2
c
      ierr = 0
  990 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDPOSCAR
c
c
      subroutine RDPOSCAR( ierr, IWR, IRD, IGFL,
     $ rprim, iatmfmt, natmnm, lvatm, rscalev,
     $ ntyp,ntypd, typnm, natm  ,natmd, atmnm  ,
     $ itypa , ratm   )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/output declarations:
      INTEGER        IGFL, iatmfmt, natmnm, ntyp,natm
      INTEGER        itypa(*)
      CHARACTER*(*)  typnm(*)
      CHARACTER*(*)  atmnm(*)
      LOGICAL        lvatm
      DIMENSION      rscalev(3)
      DIMENSION      rprim(3,3), ratm(3,*)
c
c Local declarations
      PARAMETER      ( nline=120 )
      INTEGER        nline
      CHARACTER*(nline)  rdline
      CHARACTER          label*(6), opt*(1)
      LOGICAL            printon
c
c >>>> EXECUTABLE CODE
c
      ierr = 0
      printon = .false.
      if( IWR.gt.0 ) printon = .true.
c
 8006 format(a6)
 8010 format(a)
c
 9021 format(a)
 9022 format(1x,f12.8,10x,a )
 9023 format(1x,3f12.8,10x,a )
c
c >>>> Read POSCAR line 1: title line
      read(IGFL,8010)  rdline
      call LNENDS( rdline, nline, nl1, lenline )
      if( printon ) write(IWR,9021)  rdline(1:lenline)
c
c >>>> Read POSCAR line 2: scale factor
      read(IGFL,*)  scalefac
      if( printon )then
        write(IWR,9022)  scalefac, ' = scale factor'
      endif
      do  i=1,3
        rscalev(i) = scalefac
      enddo
c
c >>>> Read POSCAR line 3: read supercell vectors
      do  ivec=1,3
        read(IGFL,*)  ( rprim(i,ivec), i=1,3 )
      enddo
      if( printon )then
        write(IWR,9023)  ( rprim(1,i),i=1,3 ), ' = A'
        write(IWR,9023)  ( rprim(2,i),i=1,3 ), ' = B'
        write(IWR,9023)  ( rprim(3,i),i=1,3 ), ' = C'
      endif
c
c >>>> (VASP5+) Read POSCAR line 4(a): read atom type names
c >>>> Read POSCAR line 4(b): read counts of each atom type
      read(IGFL,8010)  rdline
      call LNENDS( rdline, nline, nl1, lenline )
      if( printon )  write(IWR,9021)  rdline(1:lenline)
c     First, prompt use for the number of types, so we know what to look for:
  410 continue
      write(IWR,*) 'Enter number of atom types in this poscar file:'
      read(IRD,*)  ntyp
      if( ntyp.le.0 .or. ntyp.gt.ntypd )then
        write(*,*) 'ERROR: ntyp out-of-bounds, 1:ntypd=', ntypd
        ierr = 42
        goto 410
      endif
c     Preload (default) atom type names:
      label(1:6) = 'TY00  '
      do  ityp=1,ntyp
c       We have already check that ityp>0
        if( ityp.lt.10 )then
          write(label(4:4),'(i1)') ityp
        elseif( ityp.lt.100 )then
          write(label(3:4),'(i2)') ityp
        endif
        typnm(ityp) = label
      enddo
c
C     ... and then parse out the type names
      read( rdline, *, err=425 )  natmt
      goto 430
c
c     Failed integer read, hence, we read in the type *names*
  425 continue
      il1 = 1
      do  ityp=1,ntyp
c       parse of start/end of next string (type name) in line
        call STRPARS( rdline,il1,lenline, n1str,ndstr )
        lenstr = ndstr - n1str + 1
        if( lenstr .gt. natmnm ) lenstr = natmnm
        typnm(ityp) = rdline(n1str:ndstr)
        il1 = ndstr + 1
      enddo
c
      read(IGFL,8010)  rdline
      call LNENDS( rdline, nline, nl1, lenline )
      if( printon )  write(IWR,9021)  rdline(1:lenline)
c
  430 continue
      il1 = 1
      natm = 0
      do 440  ityp=1,ntyp
c       Parse start/end of next string (type counts) in line
c       parse of start/end of next string (type name) in line
        call STRPARS( rdline,il1,lenline, n1str,ndstr )
        lenstr = ndstr - n1str + 1
        il1 = ndstr + 1
        ierr = 43
C        write(*,*) 'reading atom count from string',
C     $    '(', n1str,':',ndstr, ')',
C     $    ': >',rdline(n1str:ndstr),'<'
        read(rdline(n1str:ndstr),*,err=1399)  natmt
        if( natmt .le. 0 )then
c         Houston, we have a problem
          write(*,*)  '# of atom types is out-of-bounds:',natmt
          ierr = 44
          goto 1399
        endif
        natm1 = natm + 1
        natm  = natm + natmt
        if( natm .le. natmd )then
          label = 'AT0000'
          do  iatm=natm1,natm
            itypa(iatm) = ityp
c           Preload a (default) atom label
            if(     iatm.lt.10 )then
              write(label(6:6),'(i1)')  iatm
            elseif( iatm.lt.100 )then
              write(label(5:6),'(i2)')  iatm
            elseif( iatm.lt.1000 )then
              write(label(4:6),'(i3)')  iatm
            elseif( iatm.lt.10000 )then
              write(label(3:6),'(i4)')  iatm
            endif
            atmnm(iatm) = label
          enddo
        endif
c
  440 continue
      write(*,*)  'Total number of atoms to be processed=',natm
      if( natm.gt.natmd )then
        write(*,*) 'ERROR: too many atoms, max natmd=',natmd
        ierr = 49
        goto 1399
      endif
c
c >>>> (optional) Read POSCAR line 5: "Selective dynamics"?
      read(IGFL,8010)  rdline
      call LNENDS( rdline, nline, nl1, lenline )
      if( printon ) write(IWR,9021)  rdline(1:lenline)
      opt = rdline(1:1)
      if( opt .eq. 'S' .or. rdline .eq. 's' )then
c       We do not need to use this info, but need to read next input line
        read(IGFL,8010)  rdline
        call LNENDS( rdline, nline, nl1, lenline )
        if( printon ) write(IWR,9021)  rdline(1:lenline)
        opt = rdline(1:1)
      endif
c
c >>>> Read POSCAR line 6: read kind of input (direct, or cartesian)
      if( opt .eq. 'D' .or. opt .eq. 'd' )then
c       Direct (or "lattice" ) coordinats
        lvatm = .true.
      elseif( opt .eq. 'C' .or. opt .eq. 'c' )then
c       Cartesian coordinates
        lvatm = .false.
      else
c       We have a format error:
        write(*,*)  'Error trying to read coord-type line'
        ierr = 60
        goto 1399
      endif
c
c     And do not output anything as it goes, either ...
      iatmfmt = 0
      IWRX = -1
      call CONFIGRD( IGFL,IWRX, iatmfmt, natmnm, ntyp,typnm,
     $ natm  , atmnm  , itypa , ratm   )
c
c We have preloaded type names and atom names, use labeled atom format
      iatmfmt = 2
c
  990 continue
      ierr = 0
      goto 999
c
 1399 continue
      write(*,*)  'ERROR: processing poscar file, err=',ierr
      if( ierr.eq.0 ) ierr = -1
c
c    That's all, Folks!
c
  999 continue
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LNENDS
c
c
      subroutine LNENDS( line, nline, n1line, ndline )
c
      IMPLICIT NONE
c
      INTEGER    nline, n1line, ndline
      CHARACTER*(*)  line
c Local declarations:
      INTEGER    ic
c
c Find last non-blank character in line:
      ndline = 0
      do 10 ic=nline,1,-1
        if( line(ic:ic) .ne. ' ' )then
          ndline = ic
          goto 11
        endif
   10 continue
c     Line is blank
      n1line = 0
      goto 999
   11 continue
c
c Find first non-blank character in line:
      n1line = 1
      do 20 ic=1,ndline
        if( line(ic:ic) .ne. ' ' )then
          n1line = ic
          goto 21
        endif
   20 continue
   21 continue
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SCALECELL
c
c
      subroutine SCALECELL( ndim, lvatm, rscale,
     $                natmuc, rprim, ratmuc, orig,offset )
c---------------------------------------------------------------
c Purpose: project vector r() onto linearly independent vectors
c          a(), b(), c(), putting result (Cramer's rule) in v().
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      LOGICAL    lvatm
      DIMENSION  rscale(3)
c Input/output arrays:
      DIMENSION  rprim(3,3), ratmuc(3,*), orig(3), offset(3)
c
c >>>> EXECUTABLE CODE
c
cxxx:
      write(*,*) 'SCALECELL: ndim,lvatm,rscale=',
     $ ndim,lvatm,rscale
      if( ndim.le.0 .or. ndim.gt.3 ) RETURN
c
      do  j=1,3
        do  i=1,3
          rprim(i,j) = rscale(i)*rprim(i,j)
        enddo
      enddo
c
c     Install scaling into atomic positions:
c
c Only project along cartesian dimensions
      mdim = 1
      if( lvatm ) mdim = ndim + 1
      if( mdim.lt.1 .or. mdim.gt.3 ) RETURN
c
cxxx:start
      write(*,*) 'SCALECELL: mdim,natmuc=',ndim,mdim,natmuc
  913 format(1x,a,i5,3f16.6)
cxxx:end
        do  iatm=1,natmuc
cxxx:
C      if( iatm.le.5 ) write(*,913)'SC-b4',iatm,(ratmuc(i,iatm),i=1,3)
          do  i=mdim,3
            ratmuc(i,iatm) = rscale(i)*ratmuc(i,iatm)
          enddo
cxxx:
C      if( iatm.le.5 ) write(*,913)'SC-aft  ',iatm,(ratmuc(i,iatm),i=1,3)
        enddo
        do  i=mdim,3
          orig(i) = rscale(i)*orig(i)
          offset(i) = rscale(i)*offset(i)
        enddo
c
c    That's all FOlks
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVPROJ
c
c
      subroutine LVPROJ( a, r, v )
c---------------------------------------------------------------
c Purpose: project vector r() onto linearly independent vectors
c          a(), b(), c(), putting result (Cramer's rule) in v().
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  a(3,3), r(3)
c Output arrays:
      DIMENSION  v(3)
c
c >>>> EXECUTABLE CODE:
c
      denom =  a(1,1)*( a(2,2)*a(3,3) - a(3,2)*a(2,3) ) +
     $         a(1,2)*( a(2,3)*a(3,1) - a(3,3)*a(2,1) ) +
     $         a(1,3)*( a(2,1)*a(3,2) - a(3,1)*a(2,2) )
c
      v(1) = ( r(1)  *( a(2,2)*a(3,3) - a(3,2)*a(2,3) ) +
     $         a(1,2)*( a(2,3)*r(3)   - a(3,3)*r(2)   ) +
     $         a(1,3)*( r(2)  *a(3,2) - r(3)  *a(2,2) ) ) / denom
      v(2) = ( a(1,1)*( r(2)  *a(3,3) - r(3)  *a(2,3) ) +
     $         r(1)  *( a(2,3)*a(3,1) - a(3,3)*a(2,1) ) +
     $         a(1,3)*( a(2,1)*r(3)   - a(3,1)*r(2)   ) ) / denom
      v(3) = ( a(1,1)*( a(2,2)*r(3)   - a(3,2)*r(2)   ) +
     $         a(1,2)*( r(2)  *a(3,1) - r(3)  *a(2,1) ) +
     $         r(1)  *( a(2,1)*a(3,2) - a(3,1)*a(2,2) ) ) / denom
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVXPND
c
c
      subroutine LVXPND( a, r, v )
c---------------------------------------------------------------
c Purpose: expand vector v() into Cartesian r() using basis a()
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  a(3,3), v(3)
c Output arrays:
      DIMENSION  r(3)
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
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATFMT
c
c
      subroutine ATFMTIN( iatmfmt, line )
c---------------------------------------------------------------
c Purpose: connect type of atom-line format to atom-line label
c
c Written: Peter A. Schultz, 11-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  "coordinates"
c         = 1  "atom, type, positions"
c
c Warning: routine assumes (does not check) that input "line"
c is at least six characters long, and that the output "line"
c has enough space to load the identifying label.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
      CHARACTER*(*)  readlbl*6
c
c >>>> EXECUTABLE CODE:
c
      iatmfmt = -1
c
      readlbl = line(1:6)
      if( readlbl.eq. 'coordi' )then
c       This will be position-only input (x,y,z)
        iatmfmt = 0
      elseif( readlbl.eq. 'atom, ' )then
c       This denotes atom,type,position (at,ty,x,y,z) input
        iatmfmt = 2
      endif
      RETURN
c
      entry ATFMTOUT( iatmfmt, line )
c
      if( iatmfmt.eq.0 )then
        line = 'coordinates (x y z) for atoms'
      elseif( iatmfmt.eq.1 .or. iatmfmt.eq.2 )then
        line = 'atom, type, position (at,ty,x,y,z) for atoms'
      else
C        call STOPXERR( 'ATFMTOUT: unknown atom output format' )
        write(*,*)     'ATFMTOUT: unknown atom output format'
        STOP 'atom-fmt'
      endif
      RETURN
c
c >>>> That's all, Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGRD
c
c
      subroutine CONFIGRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  10Jan02-PAS/2.54: bare coord read enabled
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  just x,y,z
c         = 1  atom, type(as number), coords
c         = 2  atom name, type(as name), coords
c         = 3  types previously unknown
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
      DIMENSION  itypa(*)
      DIMENSION  ratm(3,natm)
c
c >>>> EXECUTABLE CODE:
c
      do  iatm=1,natm
        call ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGWR
c
c
      entry CONFIGWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        call ATOMWR( IDT, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDRD
c
c
      subroutine COORDRD( IDT, iatmfmt, natm,ratm )
c---------------------------------------------------------------
c Purpose: extract atomic coordinates from geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c atomic coordinates:
      DIMENSION  ratm(3,natm)
c Local
      DIMENSION  ratmi(3)
c
c >>>> EXECUTABLE CODE:
c
      jatmfmt = iatmfmt
      do  iatm=1,natm
        call RATMRD( IDT, jatmfmt, ratmi )
        ratm(1,iatm) = ratmi(1)
        ratm(2,iatm) = ratmi(2)
        ratm(3,iatm) = ratmi(3)
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDWR
c
c
      entry COORDWR( IDT,  natm,ratm )
c---------------------------------------------------------------
c Purpose: dump raw atomic coordinates to geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        ratmi(1) = ratm(1,iatm)
        ratmi(2) = ratm(2,iatm)
        ratmi(3) = ratm(3,iatm)
        call RATMWR( IDT, ratmi )
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMRD
c
c
      subroutine ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Oh, the ineffable joys of character manipulation in fortran ...
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  itypa(*)
      DIMENSION  ratm(3)
c
c local declarations:
      PARAMETER  ( nline=80 )
      CHARACTER  line*(nline), strtmp*(nline)
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = ( IWR .ge. 0 )
c Global format for output of atomic coordinates:
 9100 format(1x,3f16.10)
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
        if( printon ) write(IWR,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 .or.
     $        iatmfmt .eq. 3 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        if( printon ) write(IWR,'(a)' )   line
        nl1 = 1
        nl2 = nline
c
c       Parse the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        if( n1 .eq. 0 ) STOP 'atm-name'
        nmlen = n2 - n1 + 1
C       if( nmlen .gt. natmnm ) call STOPXERR( 'atom name too long' )
        if( nmlen .gt. natmnm ) STOP 'longatom'
        strtmp = ' '
        strtmp(1:nmlen) = line(n1:n2)
        atmnm = strtmp(1:natmnm)
        nl1 = n2 + 1
c
c       First, check if atoms input as number sequence.
c       If they are, check that atoms are indeed in sequence.
        read( atmnm, * , err=120 ) jatm
C        if( jatm .ne. iatm ) call STOPXERR( 'atom out of sequence' )
CCCC        if( jatm .ne. iatm )  STOP 'atom-seq'
c       Shift atom name to right aligned 6-char deep
        atmnm = ' '
        write( atmnm, '(i6)' ) iatm
c
c       Parse the atom type:
c
  120   continue
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        if( n1 .eq. 0 )  STOP 'no type'
        nmlen = n2 - n1 + 1
C        if( nmlen .gt. natmnm ) call STOPXERR( 'type name too long' )
        if( nmlen .gt. natmnm ) STOP 'longtype'
        nl1 = n2 + 1
c
c       And now try and get the atom type:
c
        if( iatmfmt .eq. 2 .or. iatmfmt .eq. 3 )then
          itypa(iatm) = 0
          do  ityp=1,ntyp
            strtmp = typnm(ityp)
            call STRPARS( strtmp, 1,nline, nt1,nt2 )
            ntlen = nt2 - nt1 + 1
            if( ntlen .eq. nmlen .and.
     $          line(n1:n2) .eq. strtmp(nt1:nt2) )then
              itypa(iatm) = ityp
              goto 200
            endif
          enddo
        endif
        if( iatmfmt .eq. 3 )then
c         Declare a new type, and name it
c         **** WARNING **** ntyp dim not checked to be adequate
          ntyp = ntyp + 1
          write(*,9125) iatm,ntyp,line(n1:n2)
 9125     format(1x,'>>>> Atom#',i4,' - new type#',i2,' named >',a,'<')
          typnm(ntyp) = line(n1:n2)
          itypa(iatm) = ntyp
          goto 200
        endif
c       Did not find type in list of type names, perhaps as type #?
        read( line(n1:n2), * , err=130 ) ityp
C        if( ityp.lt.1 .or. ityp.gt.ntyp ) call STOPXERR( 'bad type' )
        if( ityp.lt.1 .or. ityp.gt.ntyp ) STOP 'bad type'
        itypa(iatm) = ityp
c
c       Atom type are input as numbers, not names:
        iatmfmt = 1
c
        goto 200
c
C  130   call STOPXERR( 'atom type unknown' )
  130   continue
        write(*,'(a,a)') 'type input=',line(n1:n2)
        write(*,*) 'number of types:',ntyp
        do  ityp=1,ntyp
           write(*,'(a,i2,1x,a,a)') 'type#',ityp,'=',line(n1:n2)
        enddo
        STOP 'unk-type'
c
  200   continue
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (input) format' )
        STOP 'geomform'
      endif
c
      RETURN
Cc
Cc >>>> ERROR HANDLING
Cc
C1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
C     RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMWR
c
c
       entry ATOMWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Write out raw coordinates without atomic id's
        write(IDT,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c       Build line with the atomic information:
        line = ' '
        nl1 = 2
        nl2 = nline
c
c       Parse and incorporate the atom name:
        call STRPARS( atmnm, 1,natmnm, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'nameless atom in atom output' )
        if( n1 .eq. 0 ) STOP 'no_atmnm'
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = atmnm(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Parse and incorporate the atom type name:
        ityp = itypa(iatm)
        strtmp = ' '
        if( iatmfmt .eq. 1 )then
c         Output type as type #, not type name
          n2 = 6
          write( strtmp(1:6), '(i6)' ) ityp
        else
c         Output type as type name
          call STRPARS( typnm(ityp), 1,natmnm, n1,n2 )
C          if( n1 .eq. 0 ) call STOPXERR( 'nameless type in output' )
          if( n1 .eq. 0 ) STOP 'no_typnm'
          strtmp(1:natmnm) = typnm(ityp)
        endif
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = strtmp(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Incorporate the atom coordinates:
        nl2 = nl1 + 1 + 3*16
C        if( nl2 .gt. nline ) call STOPXERR( 'atom type/coord too long' )
        if( nl2 .gt. nline ) STOP 'longatll'
        write( line(nl1:nl2), 9100 )  ratm
c
        write(IDT,'(a)')  line(1:nl2)
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (output) format' )
        STOP 'geom-fmt'
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMRD
c
c
      entry RATMRD( IDT, iatmfmt, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        nl1 = 1
        nl2 = nline
c
c       Skip the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        if( n1 .eq. 0 ) STOP 'no-atmnm'
        nl1 = n2 + 1
c
c       Skip the atom type (second string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
C        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        if( n1 .eq. 0 ) STOP 'no-typin'
        nl1 = n2 + 1
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
C        call STOPXERR( 'unknown geometry (input) format' )
        STOP 'geom-fmt'
c
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMWR
c
c
      entry RATMWR( IDT, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      write(IDT,9100)  ratm
c
      RETURN
c
c >>>> ERROR HANDLING
c
C 1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
 1301 STOP 'err-ratm'
      RETURN
c
c    That's all Folks!
c
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
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRNUM
c
c
      subroutine STRNUM( nstr, str, num, maxnum )
c---------------------------------------------------------------
c Purpose: put integer into string of appropriate length
c
c Written: Peter A. Schultz, 11-July-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Integer "num" to be written must be greater than or equal to zero.
c Integer "num" is written to character string of minimum length
c     needed to write largest possible input integer "maxnum"
c     with leading zeroes as needed.
c If error (num<0.or.maxnum<num), return null, i.e. 0-len, string
c If num=maxnum=0, return null string
c
c E.g. STRNUM( nstr, str, 13, 1000 ) returns str='0013' (nstr=4)
c      STRNUM( nstr, str,  0,   13 ) returns str='00'   (nstr=2)
c      STRNUM( nstr, str,  0,    0 ) returns nstr=0 (zero-len str)
c      STRNUM( nstr, str,113,  113 ) returns str='113'  (nstr=3)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output: character string with number
      CHARACTER*(*)  str
c
c Local declarations
      CHARACTER*(1)  chzero
      DATA           chzero / '0' /
c
c >>>> EXECUTABLE CODE:
c
      nstr = 0
c
c Check that number is within its allowed bounds:
      if( num.lt.0 .or. num.gt.maxnum .or. maxnum.eq.0 ) RETURN
c
c Determine length of full string:
      istr = 0
      mx = 1
   10 if( mx.gt.maxnum ) goto 19
      istr = istr + 1
      str(istr:istr) = chzero
      mx = mx*10
      goto 10
   19 nstr = istr
c
c Write in characters for this number:
      n0 = num
      istr = nstr
   20 nx = n0 / 10
      id = n0 - 10*nx
      write( str(istr:istr), '(i1)' ) id
      n0 = nx
      istr = istr - 1
      if( n0.gt.0 ) goto 20
c
      RETURN
c
c    That's all Folks!
c
      END
