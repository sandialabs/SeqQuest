c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SETDATA
c
c
      subroutine SETDATA( IDAT,IWR, do_gga,do_spin,
     $ c6vdw,r0vdw, c6ijvdw,r0ijvdw,
     $ do_field, efield, dielec,
     $ dokdefct, nedefct,
     $ ndim, n1r,n2r,n3r, nearopt,near2c,
     $ iatmfmt, natmnm,typnm,atmnm,
     $ natm,natmd,ntyp,ntypd,nshld,nald, noad,
     $ itypa, numshl,lshel,nala,ala,cala,  occsh, norba, znuc,
     $ atmass,atengy,
     $ dolvatm, ratm, rprim,dorig,rscale,
     $ nk,nkd, veckg,wtk,wtkscal,nku,iktype,dokgrid,
     $ nsymd,nsymi,ndimsy, rsym,isymtyp,syminfo,
     $ spinpol, elchrg,rchrg,ion_opt,
     $ any1ctr, lmx1ctr,almnv1c,lmxnlp1,almnnl,
     $ nrd,nrad,nrps,radmsh,radwt,vpsrad,
     $ anycore, nrcor,corden,cordrv,
     $ anyxc1c,nfityp,alxcmin,nalfxc,nrxcfit,nalfxcd,alfxc,alsq,ipvt,
     $ wkrd )
c---------------------------------------------------------------
c Purpose: get all data necessary for setup phase of lcao code
c
c Written: Peter A. Schultz, 15-January-2001,  for SeqQ 2.46
c
c Revision history:
c   8Apr15-PAS/2.66: character limited line-writes; lcao.geom_in!!!
c   6Apr13-PAS/2.65: add vdW options
c  12May11-PAS/2.62: add flags for slab double-layer countercharge
c  23Sep08-PAS/2.62: cell volume optimization
c   5Dec07-PAS/2.60c: to_lattice->rsym,dorig uninit fixed
c  30Jun06-PAS/2.60: expand format for data echo
c  14Feb05-PAS/2.59: defect k-sampling
c  20Oct03-PAS/2.58: cosmetic change; 19Feb04: clean-up hex vectors
c  17Aug03-JTK+PAS/2.57: add external electric field; dielectric
c   8Aug03-PAS/2.56: add special offset hexagonal grid;
c                    k-vectors renamed and put in lattice units
c   5Aug03-PAS/2.55a: fix k-vector/lattice-units input
c   5Jan02-PAS/2.54: call to "configwr" altered
c  24May02-PAS/2.53e: hack: cartesian<->lattice atom conversion
c  16May02-PAS/2.53: cutoff input moved to "rundata"
c   4Apr02-PAS/2.52b: more string handling modifications
c   7Mar02-PAS/2.52: N/U occij cleaned out/isymi init
c  26Nov01-PAS/2.51: configuration files; file/unit mgr
c   3Oct01-PAS/2.50: input atom masses and reference energies
c  24Jul01-PAS/2.49: spin-polarized dft/input polarization/wtk
c  21JUn01-PAS/2.48: replace STOPs
c   5May01-PAS/2.47: fix error handling for file open, clean echo
c  15Jan01-PAS/2.46: adapted from old "indata" routine
c---------------------------------------------------------------
c
c Notes:
c  This routine reads the data which defines the problem - dimension,
c  lattice parameters, atoms and positions, BZ sampling, etc.  The
c  input is driven by six-letter "labels".  This aids several things.
c  Checks that input file is correctly assembled, allows for the input
c  of optional data, and allows the input file to be self-documenting
c  for the user.  Most of the data is read in free-format, and then
c  (almost) perfectly echoed (formatted) as it comes in to aid in
c  debugging.  For the most part, the input of data is in a very
c  strict order (the input of run phase, geometry phase, and neb input
c  in "rundata", "geomdat", and "nebdat"), with optional data.
c  There is required input (latvecs, atoms and positions, etc.), and
c  optional input.  Some of the optional input involves required
c  data but for which hopefully sensible defaults are in place, and
c  cutoffs), and others which fuel capabilities (e.g. BZ samples,
c  and others which simply make input easier (scaling).
c  Some of the key optional inputs are echoed even if not in input
c  deck, simply to illustrate their use to a user.  The trickiest
c  part of the input is keeping track of where the last read was
c  done so that you know whether to read a new label, or check one
c  you have already read previously.  This should get cleaned up.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c flags to evaluate gga gradients, and spin-polarized:
      LOGICAL    do_gga, do_spin
      LOGICAL    do_field
c flag to select k-sampling for defect states
      LOGICAL    dokdefct
c flags for auto-kgrid, lattice coordinates, clean up vectors
      LOGICAL    dokgrid, dolvatm, do_clean
cvdw vdw-defining information
      LOGICAL    DO_VDW
      DIMENSION  c6vdw(*),r0vdw(*)
      DIMENSION  c6ijvdw(*),r0ijvdw(*)
c electric field
      DIMENSION  efield(3)
c number of defect states of up and of dn spin
      DIMENSION  nedefct(2)
c scaling parameters, primitive lattice vectors, origin offset
      DIMENSION  rscale(3),rprim(3,4),dorig(3)
c atom coordinates, charge/wigner-seitz coordinates, center of symmetry
      DIMENSION  ratm(3,natmd), rchrg(3,*), rsym(3)
c
c symmetry definitions
      DIMENSION  isymtyp(*),syminfo(6,*)
c Bloch vectors, wts, and defining params for auto-kgrid
      DIMENSION  veckg(3,nkd),wtk(nkd),nku(3)
c
c basis and potential information for atoms:
      CHARACTER*(*)  typnm(ntypd), atmnm(natmd)
      DIMENSION  itypa(natmd), numshl(ntypd),lshel(nshld,ntypd),
     $ nala(nshld,ntypd),ala(nald,nshld,ntypd),cala(nald,nshld,ntypd)
      DIMENSION  norba(ntypd)
      DIMENSION  znuc(ntypd),occsh(nshld,ntypd)
      DIMENSION  atmass(ntypd),atengy(ntypd)
      DIMENSION  lmxnlp1(ntypd),almnnl(ntypd)
      DIMENSION  nrad(ntypd),nrps(ntypd),
     $ radmsh(nrd,ntypd),radwt(nrd,ntypd),vpsrad(nrd,4,ntypd)
      LOGICAL    anycore
      DIMENSION  nrcor(ntypd), corden(nrd,ntypd),cordrv(nrd,ntypd)
c dynamic 1-center definitions(non-standard):
      LOGICAL    any1ctr
      DIMENSION  lmx1ctr(ntypd),almnv1c(ntypd)
c xc-fit variables (nonstandard/non-functional developmental code):
      LOGICAL    anyxc1c
      DIMENSION  nfityp(ntypd),alxcmin(ntypd)
      DIMENSION  nalfxc(ntypd),nrxcfit(ntypd),alfxc(nalfxcd,ntypd),
     $ alsq(nalfxcd,nalfxcd,3,ntypd),ipvt(nalfxcd,3,ntypd)
c
c Scratch space:
      DIMENSION  wkrd(*)
c
c Local declarations:
      PARAMETER  ( nline = 128 )
      CHARACTER*(nline)  label
      CHARACTER*128 atmfil,filenm
      CHARACTER*6  wantlbl,readlbl
      CHARACTER*1  ktype
      EQUIVALENCE  (label,readlbl),(label,atmfil)
      DIMENSION    vtmp(4)
      DIMENSION    gvec(3,3), strain(3,3),tmp(3,3)
c
      DATA  toev / 13.605805d0 /, toang / 0.5291771d0 /
      DATA  zero,half,one,two,three / 0.d0,0.5d0,1.d0,2.d0,3.d0 /
      DATA  vsmall / 1.d-9 /
      DATA  degpi / 180.d0 /
c
c >>>> EXECUTABLE CODE:
c
      do_clean = .true.
      rtpid2 = rtpi/two
c
      write(IWR,'(a)') 'setup data:'
c
c >>>>>> TITLE LINES/NOTES
c (optional)
      read(IDAT,8000)  label
c
      if( readlbl(1:5) .eq. 'title' )then
        call LNWRIT( IWR, label, nline )
        read( readlbl(6:6), '(i1)' )  ntitl
        if( ntitl.eq.0 ) ntitl = 1
        write(IWR,*) 'Title line(s) for this run:'
        if( ntitl.gt.0 )then
          do  ititl=1,ntitl
            read(IDAT,8000)  label
            call LNWRIT( IWR, label, nline )
          enddo
        endif
        read(IDAT,8000)  label
      endif
c
      if( readlbl(1:5) .eq. 'notes' )then
        call LNWRIT( IWR, label, nline )
        ntitl = 20
        do 100 ititl=1,ntitl
          read(IDAT,8000)  label
          call LNWRIT( IWR, label, nline )
          if( readlbl .eq. 'end_no' ) goto 101
  100   continue
        write(IWR,*) 'end_notes not found after line limit=',ntitl
        write(IWR,9013) 'endnotes - need end tag for notes: end_notes'
        goto 1310
  101   continue
        read(IDAT,8000)  label
      endif
c
c >>>>>> DFT FCNAL TYPE (capz=lda,pw91=gga,pbe,blyp)
c (optional, default=lda)
      wantlbl = 'functi'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
c
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
c
        call STRPARS( label, 1,nline, n1,n2 )
        if( n1 .eq. 0 )then
           write(IWR,9013) 'nofcnal - functional label is blank'
           goto 1310
        endif
        wantlbl = '      '
        nlbl = MIN( n2-n1+1, 6 )
        wantlbl(1:nlbl) = label(n1:n2)
c
        read(IDAT,8000)  label
      else
        write(IWR,*) 'functional type, default=CAPZ'
        wantlbl = 'LDA   '
        write(IWR,*)  '  ',wantlbl
      endif
c
c Set functional type:
      call DFTSET( wantlbl )
      call DFTGETGGA( do_gga )
      call DFTGETSPIN( do_spin )
c
c >>>>>> SPIN POLARIZATION (net spin, strictly non-negative)
c (required, iff we have spin polarized calculation)
      spinpol = zero
      if( do_spin )then
        wantlbl = 'spin p'
        if( readlbl.eq.wantlbl )then
          call LNWRIT( IWR, label, nline )
          read(IDAT,*   )  spinpol
          write(IWR,9010)  spinpol
          if( spinpol .lt. zero )then
            write(IWR,9013) 'negspin - spinpol must be non-negative'
            goto 1310
           endif
c
          read(IDAT,8000)  label
        else
          write(IWR,*) 'spin polarization (strictly non-negative)'
          write(IWR,*   )  ' X.XXXXXX - free format'
          write(IWR,9013) 'nospin - provide spin polarization'
          goto 1310
        endif
      endif
c
c >>>>>> VDW TYPE (default=OFF, ulg=ulg12,dftd2=D2)
c (optional, default=lda)
      wantlbl = 'vdw_po'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
c
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
c
        call STRPARS( label, 1,nline, n1,n2 )
        if( n1 .eq. 0 )then
           write(IWR,9013) 'no vdw - vdw label input is blank'
           goto 1310
        endif
        wantlbl = '      '
        nlbl = MIN( n2-n1+1, 6 )
        wantlbl(1:nlbl) = label(n1:n2)
c
        read(IDAT,8000)  label
      else
        write(IWR,*) 'vdw_potential type, default=OFF '
        wantlbl = 'OFF   '
        write(IWR,*)  '  ',wantlbl
      endif
c     Set vdW type:
      call VDWSET( wantlbl )
c
c >>>>>> DEFECT STATES FOR GAMM_POINT SAMPLING
c (optional) AND experimental - method deemed inappropriate
      wantlbl = 'defect'
      if( readlbl.eq.wantlbl )then
        dokdefct = .true.
        call LNWRIT( IWR, label, nline )
c       Get number of states, of each spin, to treat at gamma:
        if( do_spin)then
          read(IDAT,*)     nedefct
          write(IWR,9020)  nedefct
        else
          read(IDAT,*)     nedefct(1)
          write(IWR,9020)  nedefct(1)
        endif
c
        read(IDAT,8000)  label
      else
        dokdefct = .false.
        nedefct(1) = 0
        nedefct(2) = 0
      endif
c
c >>>>>> EXTERNAL APPLIED ELECTRIC FIELD
c (optional)
      wantlbl = 'efield'
      if( readlbl.eq.wantlbl )then
c       Get external field, in Ry/au
        call LNWRIT( IWR, label, nline )
        read(IDAT,*   )  (efield(id),id=1,3)
        toRy = ( 1.0d+06 * toang) / (toev * 1.0d+08)
        do  id=1,3
          vtmp(id) = efield(id) / toRy
        enddo
        write(IWR,'(1x,3f10.6,a,3f10.3,a)')
     $   ( efield(id), id=1,3 ), ' Ry/au =' ,
     $   ( vtmp(id)  , id=1,3 ), ' MV/cm'
        do_field = .true.
c
        read(IDAT,8000)  label
      else
        call MKZERO( 3, efield )
        do_field = .false.
      endif
c
c >>>>>> BULK DIELECTRIC CONSTANT (for defect bulk polarization)
c (optional)
      wantlbl = 'dielec'
      if( readlbl.eq.wantlbl )then
c       Get external dielectric constant
        call LNWRIT( IWR, label, nline )
        read(IDAT,*   )  dielec
        write(IWR,9010)  dielec
c
        read(IDAT,8000)  label
      else
c       Default to turn off: set to negative value:
        dielec = -1.d0
      endif
c
c >>>>>> PROBLEM DIMENSION (0=cluster thru 3=bulk):
c (required)
      wantlbl = 'dimens'
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl )then
c       Retired "lattice dimension" as alternative keyword - 21Jun11
        read(IDAT,*   )  ndim
        write(IWR,9020)  ndim
        if( ndim.lt.0 .or. ndim.gt.3 )then
          write(IWR,9013) 'badndim - ndim must be 0,1,2,or 3'
          goto 1310
        endif
c
        read(IDAT,8000)  label
      else
        write(IWR,*) 'dimension of problem (0=cluster thru 3=bulk)?'
        goto 1300
      endif
c
c >>>>>> COORDINATES: cartesian or lattice
c (optional)
c  Set default to be Cartesian, not lattice vector, coordinates:
      dolvatm = .false.
c
      if( readlbl .eq. 'coordi' )then
        call LNWRIT( IWR, label, nline )
c
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
c
c       Parse out the type of atomic coordinates
        call STRPARS( label, 1,nline, n1,n2 )
        if( n1 .eq. 0 )then
          write(IWR,9013) 'coordtyp - choose units: LATTICE/CARTESIAN'
          goto 1310
        endif
        wantlbl = '      '
        nlbl = MIN( n2-n1+1, 6 )
        wantlbl(1:nlbl) = label(n1:n2)
c
        if( wantlbl .eq. 'lattic' .or. wantlbl .eq. 'Lattic' .or.
     $      wantlbl .eq. 'LATTIC' )then
          dolvatm = .true.
c
          read(IDAT,8000)  label
c
          if( readlbl.eq.'lvscal' )then
            call LNWRIT( IWR, label, nline )
            read(IDAT,*   )  (rprim(i,4),i=1,3)
            write(IWR,9010)  (rprim(i,4),i=1,3)
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
c         At some point, we might have a Z-matrix option
          write(IWR,9013) 'badcoord - must be LATTICE or CARTESIAN'
          goto 1310
        endif
c
      else
        dolvatm = .false.
        write(IWR,9001)  'coordinates given in (lattice,cartesian):'
        write(IWR,9001)  ' cartesian'
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
        call LNWRIT( IWR, label, nline )
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
        call LNWRIT( IWR, label, nline )
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
        call LNWRIT( IWR, label, nline )
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
        call LNWRIT( IWR, label, nline )
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
        call LNWRIT( IWR, label, nline )
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
      do  id=1,ndim
        rscale(id) = scaleu*rscale(id)
      enddo
      rscale(1) = scalex*rscale(1)
      rscale(2) = scaley*rscale(2)
      rscale(3) = scalez*rscale(3)
c
c >>>>>> STRAIN
c (optional)
      istrain = 0
      wantlbl = 'strain'
      if( readlbl.eq.wantlbl )then
        if( .not. dolvatm )then
          write(IWR,9013) 'str-nolv/strain only with lattice coords'
          goto 1310
        endif
        istrain = 1
        call LNWRIT( IWR, label, nline )
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
          call LNWRIT( IWR, label, nline )
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
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl )then
        do  idim=1,3
          read(IDAT,*   )  (rprim(i,idim), i=1,3)
C          write(IWR,9010)  (rprim(i,idim), i=1,3)
          mdim = idim
        enddo
        if( do_clean ) call LVCLEAN( rprim )
        do  idim=1,3
          write(IWR,9010)  (rprim(i,idim), i=1,3)
        enddo
c Scale primitive lattice vectors, periodic directions *ONLY*:
c This needs to be done here since k-vectors will depend on
c rprim() if input for k-vectors done in lattice coordinates.
        do  j=1,ndim
          rprim(j,1) = rscale(j)*rprim(j,1)
          rprim(j,2) = rscale(j)*rprim(j,2)
          rprim(j,3) = rscale(j)*rprim(j,3)
        enddo
        if( istrain.eq.1 )then
c         Install strain into lattice vectors
          do  j=1,3
            do  i=1,3
              tmp(i,j) = strain(1,i)*rprim(1,j)
     $                 + strain(2,i)*rprim(2,j)
     $                 + strain(3,i)*rprim(3,j)
            enddo
          enddo
          do  j=1,3
            rprim(1,j) = rprim(1,j) + tmp(1,j)
            rprim(2,j) = rprim(2,j) + tmp(2,j)
            rprim(3,j) = rprim(3,j) + tmp(3,j)
          enddo
c         Leave record of modified rprim:
          write(IWR,*) 'Strained primitive lattice vectors:'
          do  idim=1,3
            write(IWR,9010)  (rprim(i,idim),i=1,3)
            tmp(1,idim) = rprim(1,idim)
            tmp(2,idim) = rprim(2,idim)
            tmp(3,idim) = rprim(3,idim)
          enddo
          do  j=1,ndim
            tmp(j,1) = tmp(j,1) / rscale(j)
            tmp(j,2) = tmp(j,2) / rscale(j)
            tmp(j,3) = tmp(j,3) / rscale(j)
          enddo
          write(IWR,*) 'Descaled primitive lattice vectors:'
          do  idim=1,3
            write(IWR,9010)  (tmp(i,idim),i=1,3)
          enddo
        endif
      else
        write(IWR,*) 'primitive vectors?'
        goto 1300
      endif
c
c     Leave a record of the cell volume as given:
      call CROSS( rprim(1,2),rprim(1,3), vtmp )
      volum = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $           + rprim(3,1)*vtmp(3) )
      write(IWR,'(1x,a,f16.6)') '> Input Cell volume{bohr^3)=',volum
c
c >>>>>> GRID DIMENSIONS ALONG UNIT CELL SIDES
c (required)
      wantlbl = 'grid d'
      read(IDAT,8000) label
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl .or. readlbl.eq.'points' )then
        read(IDAT,*   )  n1r, n2r, n3r
        write(IWR,9120)  n1r, n2r, n3r
      else
        write(IWR,*) 'grid dimensions?'
        goto 1300
      endif
c
c >>>>>> NEARBY FUNCTION TYPE
c (optional) - expert only
c  Options: nearopt=(0=none; 1=near pt; 2=weighted sampling):
c      <0 => 2ctr and 1ctr;  >0 => 1ctr corrections only
c      Default = 2/weighted sampling, 1ctr ONLY (no 2ctr)
c
      nearopt = 2
      near2c = 1
      wantlbl = 'nearby'
      read(IDAT,8000)  label
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IDAT,*   )  nropt
        if( IABS(nropt).lt.3 ) nearopt = nropt
c
        read(IDAT,8000)  label
      else
        write(IWR,9001)  'nearby option: 0=0,1=pt,2=wt'
      endif
      write(IWR,9020)  nearopt
      if( nearopt.lt.0 )then
        near2c = 2
        nearopt = -nearopt
      elseif( nearopt.eq.0 )then
        near2c = 0
      endif
c
c >>>>>> READ NUMBER OF ATOM TYPES
c (required)
      wantlbl = 'atom t'
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl )then
        read(IDAT,*   )  ntyp
        write(IWR,9020)  ntyp
        if( ntyp.gt.ntypd )then
          write(IWR,*) '>>>>> atom types=',ntyp,', but max ntypd=',ntypd
          write(IWR,9013) 'typ#2big/ max number of atom types exceeded'
          goto 1310
        endif
      else
        write(IWR,*) 'atom types (number of)?'
        goto 1300
      endif
c
c  Initialize number of types for which xc-fit to be done:
      itypxc = 0
c
      do 2000 ityp=1,ntyp
c
c >>>>>> ATOM TYPE NUMBER/LABEL
c
c (one or the other is required)
        wantlbl = 'atom f'
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
        if( readlbl.eq.wantlbl )then
          read(IDAT,8000)     atmfil
          call LNWRIT( IWR, atmfil, nline )
c         Parse an atom name along with the atom file name:
          call ATOMNM( atmfil,nline, filenm, typnm(ityp),natmnm )
C          print*, atmfil
C          print*, filenm
C          print*, typnm(ityp)
C          print*, natmnm
          write(IWR,'(a,a)') '>>>>> atom type name=',typnm(ityp)
c         Try open as fully qualified path (e.g. in remote library):
          call FLOPEN( IATMFL, filenm, 'OLD', 'FOR', iflerr )
c         If we failed, try to open in working directory:
          if( iflerr.ne.0 ) call FLOPENF( IATMFL, filenm, 'OLD', 'FOR' )
          IPOT = IATMFL
          IBAS = IATMFL
        else
          typnm(ityp) = ' '
          IPOT = IDAT
          IBAS = IDAT
          wantlbl = 'type n'
          goto 210
        endif
c
        wantlbl = 'type n'
        read(IPOT,8000)  label
        call LNWRIT( IWR, label, nline )
  210   continue
        if( readlbl.ne.wantlbl ) goto 1300
        read(IPOT,8040)  jtyp, label(1:24)
        write(IWR,9040)  ityp, label(1:24)
        if( jtyp.ne.ityp .and. IPOT.eq.IDAT )then
          write(IWR,9013) 'typ-ordr/ atom type out of sequence'
          goto 1310
        endif
c
        call ATMDAT1( IWR, IPOT,IBAS,
     $   ityp, ntyp,nshld,nald, noad,
     $   numshl,lshel,nala,ala,cala, occsh,  norba, znuc,
     $   atmass,atengy,
     $   lmx1ctr,almnv1c,lmxnlp1,almnnl,
     $   nrd,nrad,nrps,radmsh,radwt,vpsrad, nrcor,corden,cordrv,
     $   itypxc,nfityp,alxcmin,nalfxc,nrxcfit,nalfxcd,alfxc,alsq,ipvt,
     $   wkrd )
c
        if( IBAS.ne.IDAT )then
          call FLCLOSE( IATMFL )
          IPOT = IDAT
          IBAS = IDAT
        endif
c
 2000 continue
c
      read(IDAT,8000)  label
c
c >>>>>> VDW pair potential corrections
c (optional)
      if( DO_VDW() )then
c       Attempt to load vdw parameter from internal atom type data
        call LNWRIT( IWR, label, nline )
        write(IWR,'(a)') '# DEV/SET/VDW: load internal Z-based c_vdw'
        call VDWNULL( ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
        IWRF = -1
        call VDWLOADZ( IWRF, ntyp, znuc,typnm,natmnm, c6vdw,r0vdw )
        do  ityp=1,ntyp
          write(IWR,'(a,i5,5x,a,5x,2f16.5)')
     $     '# vdwloadz: ityp,c6(Ry),r0(Bohr)=',
     $      ityp,typnm(ityp),c6vdw(ityp),r0vdw(ityp)
        enddo
c
        wantlbl = 'vdw_da'
        if( readlbl.eq.wantlbl )then
c         Obtain vdW specifications from user input
          call LNWRIT( IWR, label, nline )
          IDATV = IDAT
          call VDWDATA( IWR ,IDATV, ntyp, typnm,natmnm, c6vdw,r0vdw,
     $                  c6ijvdw,r0ijvdw )
c
          read(IDAT,8000)  label
        else
c         Use internal vdW specifications according to vdW-ff selected
          call VDWCHECK( IWR, ntyp, typnm, c6vdw,r0vdw, ierr )
          if( ierr .ne. 0 )then
            write(IWR,*) 'Total atom types lacking vdW parameters=',ierr
            call STOPXERR( 'Some atoms missing vdW parameters' )
          endif
          call VDWCROSS( ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
        endif
      endif
c
c >>>>>> ATOM MASSES
c (optional)
      wantlbl = 'masses'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        do  ityp=1,ntyp
          read(IDAT,*)  atmassi
          if( atmassi .gt. zero ) atmass(ityp) = atmassi
          write(IWR,9011)  atmass(ityp)
        enddo
        read(IDAT,8000)  label
      endif
c
c >>>>>> ATOM REFERENCE ENERGIES
c (optional)
      wantlbl = 'energi'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        do  ityp=1,ntyp
          read(IDAT,*)  atengyi
          if( atengyi .le. zero ) atengy(ityp) = atengyi
          write(IWR,9011)  atengy(ityp)
        enddo
        read(IDAT,8000)  label
      endif
c
c >>>>>> READ CHARGE STATE INFORMATION
c (optional)
cpas: ion_opt default needs to be thought out more
      ion_opt = 0
      if( ndim.lt.3 ) ion_opt = 2
      wantlbl = 'ionopt'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
c
        read(IDAT,*   )  iopt
        if( IABS( iopt ).le.3 )then
          ion_opt = iopt
        else
          write(IWR,*) 'Charge option is not valid: ',iopt
          write(IWR,9013) 'ionopterr: invalid ion option'
          goto1310
        endif
        write(IWR,'(i3,a)') ion_opt,
     $   ' = ion/dipole option: 1=flat,2=gaus,3=qslab,<0=no dipole'
c
        read(IDAT,8000) label
      else
        write(IWR,9001)  'ionopt (0=none,2=gaussian lmcc,3=dble layer):'
        write(IWR,9120)  ion_opt
      endif
c
      wantlbl = 'charge'
      elchrg = zero
      call MKZERO( 3, rchrg )
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
c
        read(IDAT,*   )  elchrg
        write(IWR,9010)  elchrg
        if( ndim.gt.0 .and. elchrg.ne.zero )then
          if( IABS( ion_opt ).eq.1 )then
c           If charged *periodic* system and jellium,
c           must use uniform jellium potl, not spherical flat potential
            ion_opt = 0
          elseif( IABS( ion_opt ) .gt. 1 )then
            read(IDAT,8000)  label
            call LNWRIT( IWR, label, nline )
            wantlbl = 'locati'
            if( readlbl.ne.wantlbl )then
              write(IWR,*) '>>>>> periodic countercharge needs location'
              goto 1300
            endif
c           *** NB ***: rchrg is used as Wigner-Seitz origin
            if( IABS( ion_opt ) .eq. 3 )then
              nchrg = 2
            else
              nchrg = 1
            endif
            read(IDAT,*   )  ((rchrg(i,ichrg),i=1,3), ichrg=1,nchrg)
            write(IWR,9010)  ((rchrg(i,ichrg),i=1,3), ichrg=1,nchrg)
          endif
        endif
        read(IDAT,8000)  label
      else
        write(IWR,9001) 'charge state for system'
        write(IWR,9010)  elchrg
      endif
c
c >>>>>> NUMBER OF ATOMS
c (required)
      wantlbl = 'number'
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl )then
        read(IDAT,*   )  natm
        write(IWR,9120)  natm
        if( natm.gt.natmd )then
          write(IWR,*) '>>>>> # of atoms=',natm,', but max natmd=',natmd
          write(IWR,9013) 'natm2big/max number of atoms exceeded'
          goto 1310
        endif
      else
        write(IWR,*) 'number of centers/atoms in system?'
        goto 1300
      endif
c
      IDATG = IDAT
  300 continue
      read(IDATG,8000)  label
c
c >>>>>> ATOM TYPES AND POSITIONS
c (required)
      wantlbl = 'atom, '
      call LNWRIT( IWR, label, nline )
      if( readlbl.eq.wantlbl )then
        iatmfmt = 2
      elseif( IDAT.eq.IDATG .and. readlbl .eq. 'geomfi' )then
c       (NB: needed to ensure that geom file does not point to a file)
        call FLOPEN( IDATG, 'lcao.geom_in', 'OLD', 'FOR', ierror )
        if( ierror .ne. 0 )
     $   call STOPXERR( 'geometry input (lcao.geom_in) not found' )
c       We found an input geometry file, go get the data
        goto 300
c      elseif( readlbl.eq.'poscar' )then
cc       Read a VASP poacar
c        call FLOPEN( IDATG, 'lcao.poscar', 'OLD', 'FOR', ierror )
c        if( ierror .ne. 0 )
c     $   call STOPXERR( 'error opening lcao.poscar file' )
c        iatmfmt = (poscar format)
c        
      else
        wantlbl = 'atom t'
        if( readlbl.ne.wantlbl )then
          write(IWR,*) 'atom, types, and positions?'
          goto 1300
        endif
      endif
      call CONFIGRD( IDATG,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c
      read(IDAT,8000)  label
c
c >>>>>> TO_LATTICE/TO_CARTESIAN - convert to lattice/cartesian units
c (optional) ... and dangerous; expert only
      wantlbl = 'to_lat'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        if( ndim.ne.3 .or. .not. dolvatm )then
          write(IWR,9013) 'to_lattice - only for bulk + lattice units'
          goto 1310
        endif
c       Install the scaling factors:
        do  iatm=1,natm
          ratm(1,iatm) = ratm(1,iatm)*rscale(1)
          ratm(2,iatm) = ratm(2,iatm)*rscale(2)
          ratm(3,iatm) = ratm(3,iatm)*rscale(3)
        enddo
c       Eliminate scale factors (this is really weird and dangerous):
        rscale(1) = one
        rscale(2) = one
        rscale(3) = one
c       Project coordinates onto primitive lattice vectors:
        do  iatm=1,natm
          call LVPROJ( rprim, ratm(1,iatm), vtmp )
          ratm(1,iatm) = vtmp(1)
          ratm(2,iatm) = vtmp(2)
          ratm(3,iatm) = vtmp(3)
        enddo
c
        read(IDAT,8000)  label
c       Skip origins, since do not do correct conversion with this hack
c       but need to initialize origins to something!!!!! 4dec07-pas
        call MKZERO( 3, dorig )
        call MKZERO( 3, rsym  )
        goto 400
      endif
c
      wantlbl = 'to_car'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        if( ndim.ne.3 .or. dolvatm )then
          write(IWR,9013) 'to_cartesian - only for bulk + NOT lattice'
          goto 1310
        endif
c
c       Expand coordinates from primitive lattice vector units:
        do  iatm=1,natm
          call LVXPND( rprim, vtmp, ratm(1,iatm) )
          ratm(1,iatm) = vtmp(1)
          ratm(2,iatm) = vtmp(2)
          ratm(3,iatm) = vtmp(3)
        enddo
c
c       Extract the scaling factors:
        do  iatm=1,natm
          ratm(1,iatm) = ratm(1,iatm)/rscale(1)
          ratm(2,iatm) = ratm(2,iatm)/rscale(2)
          ratm(3,iatm) = ratm(3,iatm)/rscale(3)
        enddo
c
        read(IDAT,8000)  label
c       Skip origins, since do not do correct conversion with this hack
c       but need to initialize origins to something!!!!! 4dec07-pas
        call MKZERO( 3, dorig )
        call MKZERO( 3, rsym  )
        goto 400
      endif
c
c     call CONFIGWR( IWR ,     iatmfmt, natmnm, ntyp,typnm,
c    $ natm, atmnm, itypa, ratm )
c
c >>>>>> ORIGIN OFFSET OF SUPERCELL
c (optional)
      wantlbl = 'origin'
      if( readlbl.eq.wantlbl )then
c       Input origin (corner) of unit cell
        call LNWRIT( IWR, label, nline )
        read(IDAT,*   )  dorig
        write(IWR,9010)  dorig
        read(IDAT,8000)  label
      else
        write(IWR,9001) 'origin offset for unit cell'
c       Shift to -1/2*rprim in non-periodic dimensions
        call MKZERO( 3, dorig )
        write(IWR,9010)  dorig
      endif
c
c >>>>>> WIGNER-SEITZ ORIGIN
c (optional)
      if( elchrg.eq.zero )then
c
        wantlbl = 'wigner'
        if( readlbl.eq.wantlbl )then
c         Input Wigner-Seitz origin (for moments calculation)
          call LNWRIT( IWR, label, nline )
          read(IDAT,*   )  (rchrg(i,1),i=1,3)
          write(IWR,9010)  (rchrg(i,1),i=1,3)
c
          read(IDAT,8000)  label
        endif
      endif
c
c >>>>>> CENTER OF SYMMETRY
c (optional)
cpas: default rsym needs to be rethought.
      wantlbl = 'symmet'
      if( readlbl.eq.wantlbl )then
c       Input center of symmetry coordinate:
        call LNWRIT( IWR, label, nline )
        read(IDAT,*   )  rsym
        write(IWR,9010)  rsym
        read(IDAT,8000)  label
      else
        write(IWR,9001) 'symmetry center'
        rsym(1) = dorig(1)
        rsym(2) = dorig(2)
        rsym(3) = dorig(3)
        write(IWR,9010)  rsym
      endif
c
  400 continue
c
c >>>>>> READ BLOCH VECTOR INFORMATION, IF PERIODIC
c (optional)
      dokgrid = .false.
      if( ndim.eq.0 )then
        nk = 1
        wtk(1) = one
        wtkscal = one
        veckg(1,1) = zero
        veckg(2,1) = zero
        veckg(3,1) = zero
      else
        call LNWRIT( IWR, label, nline )
c
        nk = 0
        if( readlbl(1:5).eq.'kgrid' )then
          ktype = readlbl(6:6)
          read(IDAT,*   )  nku(1),nku(2),nku(3)
          write(IWR,9120)  nku(1),nku(2),nku(3)
          if( ktype.eq.' ' .or. ktype.eq.'0' )then
            iktype = 0
          elseif( ktype.eq.'1' )then
            iktype = 1
          elseif( ktype.eq.'2' .or. ktype.eq.'h' )then
            iktype = 2
          elseif( ktype.eq.'3' )then
            iktype = 3
          endif
          dokgrid = .true.
          goto 431
        elseif( readlbl(1:4).eq.'khex' )then
          ktype = readlbl(5:5)
          read(IDAT,*   )  nku(1),nku(2),nku(3)
          write(IWR,9120)  nku(1),nku(2),nku(3)
c         First, trivial check if this makes sense.
          alen = SQRT( rprim(1,1)**2 + rprim(2,1)**2 + rprim(3,1)**2 )
          blen = SQRT( rprim(1,2)**2 + rprim(2,2)**2 + rprim(3,2)**2 )
          adotb = rprim(1,1)*rprim(1,2) + rprim(2,1)*rprim(2,2)
     $          + rprim(3,1)*rprim(3,2)
          diflen = ABS( DBLE( nku(1) ) - DBLE( nku(2) )*blen/alen )
          difdot = ABS( one - two*ABS( adotb ) / (alen*blen) )
          if( diflen  .gt. (1.d-6) .or. difdot .gt. (1.d-6) )then
c           This is not hexagonal.
            write(IWR,*) 'Cell vectors A and B are not hexagonal'
            write(IWR,*) 'tests (should be<1d-6):',diflen,difdot
            write(IWR,*) 'Cannot use a hexagonal k-grid'
            write(IWR,9013) 'khex - not hex cell vectors'
            goto 1310
          endif
c         Right now, provide for single simple 1x1 hex option:
          iktype = 20
          dokgrid = .true.
          goto 431
        endif
c
c >>>>>> READ BLOCH VECTOR (K-POINT) SAMPLING POINTS (weight, vectors)
c
        wantlbl = 'n bloc'
        if( readlbl.ne.wantlbl ) goto 1300
        read(IDAT,*   )  nk
        write(IWR,9120)  nk
        if( nk.gt.nkd )then
          write(IWR,*) '>>>>> have ',nk,' k-points, but max nkd =',nkd
          write(IWR,9013) 'nk2big - max number of k-points exceeded'
          goto 1310
        endif
c
        wantlbl = 'scalin'
        read(IDAT,8000)  label
        kscale = 0
        wtkscal = one
        scalekx = one
        scaleky = one
        scalekz = one
        if( readlbl.eq.wantlbl )then
c         Read scaling factors for k-pt weights and vectors
          call LNWRIT( IWR, label, nline )
          read(IDAT,*   )  wtkscal,scalekx,scaleky,scalekz
          write(IWR,9030)  wtkscal,scalekx,scaleky,scalekz
          if( wtkscal.le.zero ) wtkscal = one
          if( scalekx.le.zero ) scalekx = one
          if( scaleky.le.zero ) scaleky = one
          if( scalekz.le.zero ) scalekz = one
          read(IDAT,8000)  label
          kscale = 1
        endif
c
        lvkvec = 0
        wantlbl = 'lattic'
        if( readlbl.eq.wantlbl )then
          lvkvec = 1
          call LNWRIT( IWR, label, nline )
          read(IDAT,8000)  label
c
        elseif( readlbl .eq. 'cartes' )then
          lvkvec = 0
          call LNWRIT( IWR, label, nline )
          read(IDAT,8000)  label
        endif
c
c       Generate the BZ primitive cell in gvec() (approp. to ndim):
        call KCELL( ndim, rprim, gvec )
c
        wantlbl = 'bloch '
        call LNWRIT( IWR, label, nline )
        if( readlbl.ne.wantlbl ) goto 1300
        wtksum = zero
        do  k=1,nk
          read(IDAT,*   )  kk, wtk(k), (veckg(j,k), j=1,3)
          write(IWR,9080)  kk, wtk(k), (veckg(j,k), j=1,3)
          if( kk.ne.k )then
            write(IWR,9013) 'kvecordr/k-vectors out of order'
            goto 1310
          endif
          wtk(k) = wtk(k) / wtkscal
          wtksum = wtksum + wtk(k)
c
c         Install local scaling factors into input k-vectors
          veckg(1,k) = scalekx*veckg(1,k)
          veckg(2,k) = scaleky*veckg(2,k)
          veckg(3,k) = scalekz*veckg(3,k)
          if( lvkvec .eq. 0 )then
c           Install global scaling into Cartesian input:
            vtmp(1) = veckg(1,k) / rscale(1)
            vtmp(2) = veckg(2,k) / rscale(2)
            vtmp(3) = veckg(3,k) / rscale(3)
c           Convert k-vectors from Cartesian to lattice units:
            call LVPROJ( gvec, vtmp, veckg(1,k) )
C          else
Cc           Got vectors in lattice, expand into Cartesian vecs:
C            call LVXPND( gvec, vtmp, veckg(1,k) )
C            veckg(1,k) = vtmp(1)
C            veckg(2,k) = vtmp(2)
C            veckg(3,k) = vtmp(3)
          endif
        enddo
c       At this point, all k-points should be in lattice units.
c
c       Check that k-point sampling weights add up to one.
        if( ABS( wtksum-one ) .gt. vsmall )then
          write(IWR,*) '>>>>> sum of kvec weights=',wtksum
          write(IWR,9013) 'kvec-wts/sum of wts must be one'
          goto 1310
        endif
c
  431   continue
c
        read(IDAT,8000)  label
c
      endif
c
c >>>>>> SYMMETRY INFORMATION
c (optional)
      ndimsy = 0
      nsymi = 0
c
      wantlbl = 'symops'
      call LNWRIT( IWR, label, nline )
c
      if( readlbl.eq. wantlbl )then
c       3D-SYMMETRY DEFINITIONS
        ndimsy = 3
        read(IDAT,*   )  nsymi
        write(IWR,9020)  nsymi
        if( nsymi.gt.(nsymd-1) )then
          write(IWR,*) '>>>>> have ',nsymi,' but max nsymd=',nsymd-1
          write(IWR,9013) 'sym#3def/max num sym definitions exceeded'
          goto 1310
        endif
c
        wantlbl = 'defini'
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
        if( readlbl.ne.wantlbl ) goto 1300
c
c Format for symmetry definitions:
c
c      isymtyp  syminfo(1-3)  syminfo(4-6)
c
c  isymtyp = order of rotation(+inversion if less than 0)
c     Valid inputs: -6,-4,-3,-2,-1,1,2,3,4,6
c     NB: 1=reflection ... i.e., equivalent to "-2"= 2-fold rotation+inv
c  syminfo(1-3) = axis of rotation (Cartesian)
c  syminfo(4-6) = translation vector (lattice units)
c
        do  isymi=1,nsymi
          read(IDAT,*   )  isymtyp(isymi),(syminfo(i,isymi),i=1,6)
          if( do_clean ) call SYMCLEAN( syminfo(1,isymi) )
          write(IWR,9071)  isymtyp(isymi),(syminfo(i,isymi),i=1,6)
          istyp = isymtyp(isymi)
          jstyp = IABS(istyp)
          if( jstyp.gt.6 .or. jstyp.eq.5 .or. jstyp.eq.0 )then
            write(IWR,9013) '3dsymtyp/bad 3D symmetry type'
            goto 1310
          endif
          if( istyp .ne. -1 )then
c           Normalize axis (if not inversion which needs no axis):
            snorm = syminfo(1,isymi)**2 + syminfo(2,isymi)**2
     $            + syminfo(3,isymi)**2
            if( snorm.le.zero )then
              write(IWR,9013) 'symmetry axis has zero length'
              goto 1310
            endif
            snorm = one / SQRT( snorm )
            syminfo(1,isymi) = snorm*syminfo(1,isymi)
            syminfo(2,isymi) = snorm*syminfo(2,isymi)
            syminfo(3,isymi) = snorm*syminfo(3,isymi)
          endif
        enddo
c
      elseif( readlbl.eq. 'symmet' )then
c       2D-SYMMETRY DEFINITIONS
        ndimsy = 2
c
c Need to step carefully here as the 2D symmetry input is a casualty
c of the reworking of the symmetry in v2.30. I have to guarantee the
c new convention is solid, which means killing the old input, but must
c leave a clean diagnostic of change.  Over time this can be dropped,
c especially if total input is completely reengineered.
c
        read(IDAT,8020)  nrot, nref
        write(IWR,9020)  nrot, nref
c
        nsymi = 0
        if( nrot.gt.0 ) nsymi = nsymi + 1
        if( nref.gt.0 ) nsymi = nsymi + nref
c
        if( nsymi.lt.1 .or. nsymi.gt.(nsymd-1) )then
          write(IWR,*) '>>>>> have ',nsymi,' but min,max =',1,nsymd-1
          write(IWR,9013) 'sym#2def/max num sym definitions exceeded'
          goto 1310
        endif
c
        isymi = 0
c
        if( nrot.gt.0 )then
c         ROTATIONS
c         Check that order of rotation is reasonable
          if( nrot.eq.1 .or. nrot.eq.5 .or. nrot.gt.6 )then
            write(IWR,9013) 'bad2Drot/only 2,3,4,6-fold rotations'
            goto 1310
          endif
c         Load into symmetry convention format:
          isymi = isymi + 1
          isymtyp(isymi) = nrot
c          ... z-axis is symmetry axis, zero translation
          call MKZERO( 6, syminfo(1,isymi) )
          syminfo(3,isymi) = one
        endif
c
        if( nref.gt.0 )then
c         REFLECTIONS
c         Read angles of reflection planes with z-x plane:
          read(IDAT,*   ,err=1302)  (wkrd(iref), iref=1,nref)
          write(IWR,9030)           (wkrd(iref), iref=1,nref)
c         Check that angles are reasonable:
          do  iref=1,nref
            ang = wkrd(iref)
            if( ang.lt.zero ) ang = ang + 360.
            if( ang.lt.zero .or. ang.gt.180. )then
              write(IWR,9013) 'bad-refl - allowed angle range 0 to 180'
              goto 1310
            endif
c           Load into symmetry convention format ...
            isymi = isymi + 1
c            ... as a two-fold rotation plus inversion ...
            isymtyp(isymi) = -2
c            ... with a rotation axis normal to the reflection plane:
            call MKZERO( 6, syminfo(1,isymi) )
            ang = ang*pi/degpi
            syminfo(1,isymi) = - SIN( ang )
            syminfo(2,isymi) =   COS( ang )
          enddo
        endif
c
      else
c
        goto 470
c
      endif
c
      read(IDAT,8000)  label
c
  470 continue
c
c >>>>>> END OF SETUP SECTION OF INPUT
c
      wantlbl = 'end se'
      if( readlbl.ne.wantlbl ) goto 1300
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Do some processing for non-standard usages
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      anycore = .false.
      any1ctr = .false.
      anyxc1c = .false.
      do  iatm=1,natm
        ityp = itypa(iatm)
        if( nrcor(ityp).ne.0 ) anycore = .true.
        if( lmx1ctr(ityp).ge.0 ) any1ctr = .true.
        if( nfityp(ityp).gt.0 ) anyxc1c = .true.
      enddo
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      RETURN
c
c >>>>>> ERROR HANDLING:
c
 9013 format(1x,'ERROR: ',a)
 1300 continue
      write(IWR,*) '>>>>> Want label >',wantlbl,'<, found >',readlbl,'<'
      call STOPXERR( 'inputlbl:  format error in input file' )
c
 1301 continue
      write(IWR,*) '>>>>> atom file not found: ',atmfil
      call STOPXERR( 'atom file not found' )
c
 1302 continue
      write(IWR,*) '>>>>> Note: 2D symmetry input changed as of v2.30'
      write(IWR,*) 'New format (keyword line + one or two data lines:'
      write(IWR,*) 'symmetry [keyword]'
      write(IWR,*) 'nrot nref [order of rotation,# of reflections]'
      write(IWR,*) '(ang(iref),iref=1,nref) [reflection plane angles]'
      call STOPXERR( '2d-syms:  2D symmetry input failure' )
c
 1310 continue
      call STOPXERR( 'inputval: input error in input file' )
c
c  Data input formats
c
 8000 format(a128)
 8010 format(3f12.8)
 8020 format(4(i2,1x))
 8040 format(i2,a)
 8080 format(i2,4x,4f12.8)
c
c  Data output formats
c
 9000 format(1x,a)
 9001 format(1x,a)
 9010 format(1x,3f12.8)
 9011 format(1x,3f16.8)
 9020 format(1x,4(i2,1x))
 9030 format(1x,1p4d16.8)
 9040 format(1x,i2,a)
 9050 format(1x,2(i4,4x),3f15.8)
 9071 format(1x,i2,1x,6f12.8)
 9080 format(1x,i3,4x,f12.8,4x,3f12.8)
 9120 format(1x,3i4)
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LVCLEAN
c
c
      subroutine LVCLEAN( rprim )
c----------------------------------------------------------------
c Purpose: numerical clean up lattice vectors, esp hex.
c Written: Peter A. Schultz, 23-February-2004, for v2.58
c Revision history:
c  none
c----------------------------------------------------------------
c
c  I will clean up only strictly canonical vectors, i.e.,
c  where the first vector A is X,0,0, and the other vectors B,C
c  are rotated 60 or 120 degrees from it.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (zero=0.d0, one=1.d0, two=2.d0, three=3.d0 )
c
c Input/output primitive vectors:
      DIMENSION  rprim(3,3)
c Local declarations:
      DIMENSION  abc(3,3)
c
c >>>> EXECUTABLE CODE:
c
      rt3 = SQRT( three )
      rt3d2 = rt3 / two
c
      do  idim=1,3
        do  i=1,3
          abc(i,idim) = rprim(i,idim)
        enddo
      enddo
c
c     Check that A = (X,0,Z)
      if( ABS( abc(1,1) ) .lt. 1.d-1 .or.
     $    ABS( abc(2,1) ) .gt. 1.d-6 ) goto 999
c
      do 200 idim=2,3
        alen = SQRT( abc(1,1)**2 + abc(2,1)**2 )
c
c       Check that |B| = |A|  (in x-y plane)
        blen = SQRT( abc(1,idim)**2 + abc(2,idim)**2 )
        btest = blen / alen
        if( ABS( btest - one ) .gt. 1.d-5 ) goto 200
c
c       Check that A.B = 1/2  (i.e., w/above, implies have hex)
        adotb = ABS( abc(1,1)*abc(1,idim) + abc(2,1)*abc(2,idim) )
        btest = adotb * two / (alen*blen)
        if( ABS( btest - one ) .gt. 1.d-6 ) goto 200
c
c       Got it, but confirm that B2 = A1*rt3/2 anyway
        b2 = abc(2,idim)
        btest =  ABS( b2 ) / ( rt3d2 * alen )
        if(  ABS( btest - one ) .gt. 3.d-5 ) goto 200
c
c       We have identified a hex vector, put in the rt3/2:
        if( abc(2,idim) .gt. 0 )then
          abc(2,idim) =  rt3d2*alen
        else
          abc(2,idim) = -rt3d2*alen
        endif
  200 continue
c
      do  idim=1,3
        do  i=1,3
          rprim(i,idim) = abc(i,idim)
        enddo
      enddo
c
  999 continue
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMCLEAN
c
c
      subroutine SYMCLEAN( svec )
c----------------------------------------------------------------
c Purpose: clean up numerics of vectors with rt3 in them
c Written: Peter A. Schultz, 23-February-2004, for v2.58
c Revision history:
c  none
c----------------------------------------------------------------
c
c  I will clean up only strictly canonical vectors, i.e.,
c  where the first vector A is X,0,0, and the other vectors B,C
c  are rotated 60 or 120 degrees from it.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (zero=0.d0, one=1.d0, two=2.d0, three=3.d0 )
c
c Input/output symmetry vector:
      DIMENSION  svec(6)
c Local declarations:
      DIMENSION  smagic(4)
c
c >>>> EXECUTABLE CODE:
c
      rt3 = SQRT( three )
      rt3d2 = rt3 / two
      smagic(1) = rt3
      smagic(2) = rt3d2
      smagic(3) = one / rt3
      smagic(4) = one / rt3d2
c
c     Clean up axis vectors (rt3-related):
c
      slen = svec(1)**2 + svec(2)**2 + svec(3)**2
      if( slen .lt. 1.d-2 ) goto 500
      slen = SQRT( slen )
c
      do 100 i=1,3
        symi = ABS( svec(i) )
        do 50 im=1,4
          stest = symi / smagic(im)
          if( ABS( stest - one ) .lt. 3.d-5 )then
c           Found a magic number!
            if( svec(i) .gt. zero )then
              svec(i) =  smagic(im)
            else
              svec(i) = -smagic(im)
            endif
            goto 100
          endif
   50   continue
  100 continue
c
c     Clean up offset vector (thirds)
c
  500 continue
      do  i=4,6
        symi = svec(i)
        if( ABS( symi ) .gt. 0.32d0 .and.
     $      ABS( symi ) .lt. 0.34d0 )then
c         Set to a clean 1/3:
          if( symi .gt. zero )then
            svec(i) =  one/three
          else
            svec(i) = -one/three
          endif
        elseif( ABS( symi ) .gt. 0.65d0 .and.
     $          ABS( symi ) .lt. 0.68d0 )then
c         Set to a clean 2/3:
          if( symi .gt. zero )then
            svec(i) =  two/three
          else
            svec(i) = -two/three
          endif
        endif
      enddo
c
  999 continue
      RETURN
c
c    That's all Folks!
c
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End
