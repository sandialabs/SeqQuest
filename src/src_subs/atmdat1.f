c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATMDAT1
c
c
      subroutine ATMDAT1( IWR, IPOT,IBAS,
     $ ityp, ntyp,nshld,nald, noad,
     $ numshl,lshel,nala,ala,cala, occsh,  norba, znuc,
     $ atmass,atengy,
     $ lmx1ctr,almnv1c,lmxnlp1,almnnl,
     $ nrd,nrad,nrps,radmsh,radwt,vpsrad, nrcor,corden,cordrv,
     $ itypxc,nfityp,alxcmin,nalfxc,nrxcfit,nalfxcd,alfxc,alsq,ipvt,
     $ wkrd )
c---------------------------------------------------------------
c Purpose: read atom file (potential and basis)
c
c Written: Peter A. Schultz,  9-May-2000, for v2.45
c          Extracted from "indata"
c
c Revision history:
c   9Apr15-PAS/2.66: lnwrit lines
c  21Jun07-PAS/2.60: removed dsifa call to clean dependencies
c  25Jun02-PAS/2.53e: FLOPEN bugfix on dormant path
c   5Mar02-PAS/2.51d: fix bug with null rfitmax/nrxc
c   6Dec01-PAS/2.51: file mgmt
c   2Oct01-PAS/2.50: input atom masses and reference energies
c  21Jun01-PAS/2.48: replace STOPs
c   5May01-PAS/2.47: clean echo; purge old format
c  28Dec00-PAS/2.46: gga code; add notes
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c basis information ...
c  ... number of shells on atom type, L(shell)
      DIMENSION  numshl(ntyp),lshel(nshld,ntyp),
c  ... gaussians/shell, gaussian exponents, contraction coefficients
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
c  ... shell occupancies for reference atom
      DIMENSION  occsh(nshld,ntyp)
c  ... total number of basis on a given atom type
      DIMENSION  norba(ntyp)
c potential information ...
c  ... effective nuclear charge
      DIMENSION  znuc(ntyp)
c  ... atomic masses (for md), and atomic reference energies
      DIMENSION  atmass(ntyp), atengy(ntyp)
c  ... Lmax+1 of pseudopotential, effective range of pot as gaussian exp.
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
c  ... number of points in radial mesh, number used in pseudopotential
      DIMENSION  nrad(ntyp),nrps(ntyp)
c  ... radial mesh, mesh weights, and pseudopotential on mesh (by L)
      DIMENSION  radmsh(nrd,ntyp),radwt(nrd,ntyp),vpsrad(nrd,4,ntyp)
c  ... number of radial pts. for core density, core density, and deriv.
      DIMENSION  nrcor(ntyp), corden(nrd,ntyp),cordrv(nrd,ntyp)
c dynamic 1-center params (non-standard) ...
c  ... Lmax for dynamic fast local density, effective gaussian range
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
c xc-fit variables (NOT production code - development only)
      DIMENSION  nfityp(ntyp),alxcmin(ntyp)
      DIMENSION  nalfxc(ntyp),nrxcfit(ntyp),alfxc(nalfxcd,ntyp),
     $ alsq(nalfxcd,nalfxcd,3,ntyp),ipvt(nalfxcd,3,ntyp)
c scratch array:
      DIMENSION  wkrd(noad)
c
c local declarations:
      PARAMETER  ( nline = 128 )
      CHARACTER    label*(nline)
      CHARACTER*6  wantlbl,readlbl
      EQUIVALENCE  (label,readlbl)
      CHARACTER*1  ndist
c
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Set up some defaults ...
c  ... no pseudopotential
      lmxp1 = -1
      almnps = one
c  ... no core density
      nrcor(ityp) = 0
c  ... no dynamic 1-center terms
      lmx1c = -1
      almn1c = 9.d9
c  ... no dynamic 1-center xc terms
      nalfxi = 0
      almnxc = 9.d9
      nfityp(ityp) = 0
c
C      wantlbl = 'type n'
C      read(IPOT,8000)  label
C      call LNWRIT( IWR, label, nline )
C  100 continue
C      if( readlbl.ne.wantlbl ) goto 1300
C      read(IPOT,8040)  jtyp, (label(il), il=1,4)
C      write(IWR,9040)  ityp, (label(il), il=1,4)
C      if( jtyp.ne.ityp .and. IPOT.eq.IDAT )then
C        write(IWR,*) '***** ERROR: atom type out of sequence',jtyp
C        call STOPXERR( 'typ ordr/atom type out of sequence' )
C      endif
c
c >>>>>> START READING FILE
c
      read(IPOT,8000)  label
c
c >>>>>> NOTES
c (optional)
      if( readlbl(1:5).eq.'notes' )then
        call LNWRIT( IWR, label, nline )
        read( readlbl(6:6), '(i1)' )  ntitl
        if( ntitl.eq.0 ) ntitl = 1
        write(IWR,*) 'notes for this atom file:'
        if( ntitl.gt.0 )then
          do  ititl=1,ntitl
            read(IPOT,8000)  label
            call LNWRIT( IWR, label, nline )
          enddo
        endif
c
        read(IPOT,8000)  label
c
      endif
c
c >>>>>> ATOMIC MASS
c (optional)
      wantlbl = 'mass  '
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,8030)  atmass(ityp)
        write(IWR,9030)  atmass(ityp)
c
        read(IPOT,8000)  label
      else
        atmass(ityp) = one
      endif
c
c >>>>>> ATOMIC REFERENCE ENERGY
c (optional)
      wantlbl = 'energy'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,8030)  atengy(ityp)
        write(IWR,9030)  atengy(ityp)
c
        read(IPOT,8000)  label
      else
        atengy(ityp) = zero
      endif
c
c >>>>>> EFFECTIVE NUCLEAR CHARGE
c
      wantlbl = 'effect'
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl ) goto 1300
      read(IPOT,8030)  znuc(ityp)
      write(IWR,9030)  znuc(ityp)
c
c If no nuclear charge/potential, skip directly to basis:
      if( znuc(ityp).eq.zero ) goto 1000
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c             have a nuclear charge, hence read potential
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      read(IPOT,8000)  label
c
c >>>>>> PSEUDOPOTENTIALS: Lmax, and effective gaussian range
c (optional)
      wantlbl = 'pseudo'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,*   )  lmxp1, almnps
        write(IWR,9025)  lmxp1, almnps
c
        read(IPOT,8000)  label
      endif
c
c >>>>>> DYNAMIC 1-CENTER PARAMS: Lmax, and maximum gaussian
c (optional)
      wantlbl = '1-cent'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,*   )  lmx1c, almn1c
        write(IWR,9025)  lmx1c, almn1c
c
        read(IPOT,8000)  label
      else
c       Must have 1-ctr stuff ON to use xc-fit stuff
        goto 200
      endif
c
c >>>>>> XC-FIT PARAMS (development code only - non-functional):
c (optional)
      wantlbl = 'xc-fit'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,8130)  nalfxi, ndist, rfitmax, almnxc
        write(IWR,9130)  nalfxi, ndist, rfitmax, almnxc
        if( nalfxi.gt.nalfxcd )then
C         write(IWR,*) '***** NOTE: setting nalfxc to max=',nalfxcd
C         write(IWR,*) '>>>>> too many fit alphas:',nalfxi,'/',nalfxcd
C         write(IWR,*) '>>>>> this option currently inactive in code'
          nalfxi = nalfxcd
        endif
c
        read(IPOT,8000)  label
      endif
c
  200 continue
c
c >>>>>> FUNCTIONAL TYPE USED TO GENERATE POTENTIAL
c (optional for now)
      wantlbl = 'functi'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
        read(IPOT,8000)  readlbl
        write(IWR,*   )  readlbl, ' = dft fcnal used to generate pot'
c
        read(IPOT,8000)  label
      endif
c
c >>>>>> NUMBER OF RADIAL POINTS FOR CENTER/NON-LOCAL
c
      wantlbl = 'radial'
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl ) goto 1300
      read(IPOT,*   )  nrl, nrnl
      write(IWR,9120)  nrl, nrnl
      if( 0.gt.nrl .or. nrnl.gt.nrl .or. nrl.gt.nrd )then
         write(IWR,*)
     $       '>>>>> Need nrnl=',nrnl,'.LE.nrl=',nrl,'.LE.nrd=',nrd
         call STOPXERR( 'nrlatom /Too many radial points for atom' )
      endif
      nrad(ityp) = nrl
      nrps(ityp) = nrnl
c
c >>>>>> RADIAL MESH POINT POSITIONS
c
      wantlbl = 'mesh p'
      read(IPOT,8000)  label
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl ) goto 1300
      read(IPOT,8065)  (radmsh(nr,ityp), nr=1,nrl)
      write(IWR,9065)  (radmsh(nr,ityp), nr=1,nrl)
c
c >>>>>> RADIAL MESH POINT WEIGHTS
c
      wantlbl = 'weight'
      read(IPOT,8000)  label
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl .and. readlbl.ne.'radwts' ) goto 1300
      read(IPOT,8065)  (radwt(nr,ityp), nr=1,nrl)
      write(IWR,9065)  (radwt(nr,ityp), nr=1,nrl)
      if( readlbl .eq. 'radwts' )then
        do  nr=1,nrl
          radwt(nr,ityp) = radwt(nr,ityp)*radmsh(nr,ityp)**2
        enddo
      endif
c
      locp = lmxp1 + 1
      if( locp.gt.0 )then
c
c >>>>>> READ NON-LOCAL POTENTIALS:
c
        do 230 lp=1,locp
c
          wantlbl = 'non-lo'
          read(IPOT,8000)  label
          call LNWRIT( IWR, label, nline )
          if( readlbl.ne.wantlbl ) goto 1300
          read(IPOT,8070)  lpot, (vpsrad(nr,lp,ityp), nr=1,MIN(6,nrl))
          if( nrl.gt.6 ) read(IPOT,8065)  (vpsrad(nr,lp,ityp), nr=7,nrl)
          write(IWR,9070)  lpot, (vpsrad(nr,lp,ityp), nr=1,nrl)
          if( lpot.ne.lp-1 ) call STOPXERR( 'lpotordr' )
  230   continue
c
c Highest non-local potential is local potential; process here.
c
        write(IWR,9240) ityp, lmxp1
 9240   format(/1x,'For atom type', i3, ', l=', i2,
     $   ' potential used as local nuclear potential' )
c
        if( locp.gt.1 )then
c         Extract local potential from non-local potentials:
          do 240 lp=1,locp-1
            do  nr=1,nrl
              vpsrad(nr,lp,ityp) = vpsrad(nr,lp,ityp)
     $                           - vpsrad(nr,locp,ityp)
            enddo
  240     continue
        endif
c
c  Remove weighting so can add locol pot directly to "vatrad()" in vlocpot
        rlogfac = LOG( radmsh(nrl,ityp) / radmsh(nrl-1,ityp) )
        do  nr=1,nrl
          vpsrad(nr,locp,ityp) = vpsrad(nr,locp,ityp) /
     $                        ( radmsh(nr,ityp)*rlogfac )
        enddo
      endif
c
c >>>>>> CORE DENSITY
c (optional)
      wantlbl = 'partia'
      read(IPOT,8000,err=1000,end=1000)  label
      if( readlbl.ne.wantlbl ) goto 1010
      call LNWRIT( IWR, label, nline )
      read(IPOT,8070)  lp, (corden(nr,ityp), nr=1,MIN(6,nrl))
      if( nrl.gt.6 ) read(IPOT,8065)  (corden(nr,ityp), nr=7,nrl)
      write(IWR,9070)  lp, (corden(nr,ityp), nr=1,nrl)
      if( lp.ne.-2 ) call STOPXERR( 'lpotord2' )
c
      if( nrl.gt.nrnl )then
        do  nr=nrnl+1,nrl
          corden(nr,ityp) = zero
          cordrv(nr,ityp) = zero
        enddo
      endif
      do 270 nr=nrnl,1,-1
        nrsc = nr
        if( corden(nr,ityp) .ne. zero ) goto 271
  270 continue
  271 nrcor(ityp) = nrsc
c
c >>>>>> CORE DENSITY DERIVATIVE
c (optional, and not used)
      wantlbl = 'deriva'
      read(IPOT,8000)  label
      if( readlbl.ne.wantlbl ) goto 1010
      call LNWRIT( IWR, label, nline )
      read(IPOT,8070)  lp, (cordrv(nr,ityp), nr=1,MIN(6,nrl))
      if( nrl.gt.6 ) read(IPOT,8065)  (cordrv(nr,ityp), nr=7,nrl)
      write(IWR,9070)  lp, (cordrv(nr,ityp), nr=1,nrl)
      if( lp.ne.-3 ) call STOPXERR( 'lpotord3' )
c
 1000 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                        get basis data from file
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c >>>>>> READ BASIS FOR THIS ATOM, FIRST GET NUMBER OF SHELLS
c
      read(IBAS,8000)  label
c
 1010 wantlbl = 'number'
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl ) goto 1300
      read(IBAS,8020)  nshl
      write(IWR,9020)  nshl
      if( nshl.gt.nshld )then
        write(IWR,*) '>>>>> input=',nshl,',   max=',nshld
        call STOPXERR( 'nshl2big/ too many shells for atom basis' )
      endif
      numshl(ityp) = nshl
c
      norbat = 0
      do 1070 ishl=1,nshl
c
c >>>>>> ANGULAR MOMENTUM, NUMBER OF GAUSSIANS FOR THIS SHELL
c
        wantlbl = 'angula'
        read(IBAS,8000)  label
        call LNWRIT( IWR, label, nline )
        if( readlbl.ne.wantlbl ) goto 1300
        read(IBAS,8020)  lshl, nal
        write(IWR,9020)  lshl, nal
        if( lshl.gt.2 .or. nal.gt.nald )then
          write(IWR,*) '>>>>> l,ng=',lshl,nal,',   max=',2,nald
          call STOPXERR( 'nala2big/atom basis out of dimensioning' )
        endif
        lshel(ishl,ityp) = lshl
        nala(ishl,ityp) = nal
        norbat = norbat + 2*lshl + 1
c
c >>>>>> ALPHAS OF GAUSSIANS FOR THIS SHELL
c
        wantlbl = 'alphas'
        read(IBAS,8000)  label
        call LNWRIT( IWR, label, nline )
        if( readlbl.ne.wantlbl ) goto 1300
        read(IBAS,8030)  (ala(ial,ishl,ityp), ial=1,nal)
        write(IWR,9030)  (ala(ial,ishl,ityp), ial=1,nal)
c
c >>>>>> COEFICIENTS FOR CONTRACTED GAUSSIAN SHELLS
c
        wantlbl = 'wave f'
        read(IBAS,8000)  label
        call LNWRIT( IWR, label, nline )
        if( readlbl.ne.wantlbl ) goto 1300
        read(IBAS,8030)  (cala(ial,ishl,ityp), ial=1,nal)
        write(IWR,9030)  (cala(ial,ishl,ityp), ial=1,nal)
c
c         Close loop over shells:
c
 1070 continue
c
      if( norbat.gt.noad )then
        write(IWR,*) '>>>>> typ=',ityp,', bf/atm=',norbat,',noad=',noad
        call STOPXERR( 'noad2big/too many orbitals per atom (noad)' )
      endif
      norba(ityp) = norbat
c
c >>>>>> OCCUPATION NUMBERS FOR REFERENCE ATOM OF THIS TYPE
c
      read(IBAS,8000)  label
      call LNWRIT( IWR, label, nline )
c
c EITHER: old by-orbital occupancies
      wantlbl = 'occupa'
      if( readlbl.eq.wantlbl )then
c       Read in occupancies by orbital
        read(IBAS,8065)  (wkrd(iorb), iorb=1,norbat)
        write(IWR,9065)  (wkrd(iorb), iorb=1,norbat)
c       Convert to shell occupancies:
        iorbat0 = 0
        do  ishl=1,nshl
          norbsh = 2*lshel(ishl,ityp) + 1
          occshl = zero
          do  iorbat=iorbat0+1,iorbat0+norbsh
            occshl = occshl + wkrd(iorbat)
          enddo
          occsh(ishl,ityp) = occshl
          iorbat0 = iorbat0 + norbsh
        enddo
        goto 1100
      endif
c
c OR: new by-shell occupancies
      wantlbl = 'shell '
      if( readlbl.ne.wantlbl ) goto 1300
c     Read in occupancies by shell
      read(IBAS,8065)  (occsh(ishl,ityp), ishl=1,nshl)
      write(IWR,9065)  (occsh(ishl,ityp), ishl=1,nshl)
c
c >>>>>> END OF ATOM DATA:
c
      read(IBAS,8000,end=1100,err=1100)  label
      call LNWRIT( IWR, label, nline )
c
      wantlbl = 'end at'
      if( readlbl.ne.wantlbl ) goto 1300
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                End input of atom data files
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1100 continue
c
      lmxnlp1(ityp) = lmxp1
      almnnl(ityp) = almnps
      if( almnnl(ityp).eq.zero )  almnnl(ityp) = one
c
      lmx1ctr(ityp) = lmx1c
      almnv1c(ityp) = almn1c
      if( almnv1c(ityp).eq.zero ) almnv1c(ityp) = 9.d9
c
      nalfxc(ityp) = nalfxi
      alxcmin(ityp) = almnxc
      if( alxcmin(ityp).eq.zero ) alxcmin(ityp) = 9.d9
c
      almina = ala(1,1,ityp)
      almaxa = ala(1,1,ityp)
      do  ishl=1,nshl
        lshl = lshel(ishl,ityp)
        nal = nala(ishl,ityp)
        call GNORM( lshl, nal, ala(1,ishl,ityp),cala(1,ishl,ityp) )
        do  ial=1,nal
          almina = MIN( almina, ala(ial,ishl,ityp) )
          almaxa = MAX( almaxa, ala(ial,ishl,ityp) )
        enddo
      enddo
      alamin = almina
      alamax = almaxa
c
c Determine whether this atom type will have a 1-center potential
      if( znuc(ityp).eq.zero .or. alamax.le.almnv1c(ityp) )
     $ lmx1ctr(ityp) = -1
c
cxcf:###############################################################
cxcf:start
cxcf: development code for testing schemes for fitting scf xc-pots
cxcf: rest of this code, until "cxcf:end", can be ignored.
cxcf:
        if( lmx1ctr(ityp).ge.0 .and. nalfxc(ityp).gt.0
     $   .and. alamax.gt.alxcmin(ityp) )then
          itypxc = itypxc + 1
          nfityp(ityp) = itypxc
          do 1210 nr=1,nrl
            nrxc = nr
            if( radmsh(nr,ityp).gt.rfitmax ) goto 1211
 1210     continue
 1211     nrxcfit(ityp) = nrxc
        else
          nrxcfit(ityp) = 0
        endif
c
c Determine alphas for xc-fit if this atom gives a 1-ctr potential
C        if( lmx1ctr(ityp).ge.0 .and. alamax.gt.alxcmin(ityp) )then
        if( nfityp(ityp).gt.0 )then
c         nalfxc(ityp) alphas geom. from alamin/3 to 2*alamax/3
          if( ndist.eq.'g' )then
            alfxc(1,ityp) = alamin/3.
            if( nalfxi.gt.1 )
     $       alrat = (2.*alamax/alamin)**(1./( DBLE(nalfxi)-1.))
            do  nal=2,nalfxi
              alfxc(nal,ityp) = alfxc(nal-1,ityp)*alrat
            enddo
CPAS: TEMP TO GET ALPHAS I WANT IN
          ELSEIF( NDIST.EQ.'r' )THEN
            CALL FLOPENF( IDTXC, 'ALF.XC' ,'OLD', 'FOR' )
            READ(IDTXC,8272)  NALFXI
            NALFXC(ITYP) = NALFXI
 8272       FORMAT(I2)
            READ(IDTXC,*) (ALFXC(IALF,ITYP),IALF=1,NALFXI)
            CALL FLCLOSE( IDTXC )
          else
c           space fit alphas arith. between alamin/3 and 2*alamax/3
            alfxc(1,ityp) = alamin/3.
            if( nalfxi.gt.1 )
     $       aldif = (2.*alamax-alamin)/(3.*( DBLE(nalfxi)-1.))
            do  nal=2,nalfxi
              alfxc(nal,ityp) = alfxc(nal-1,ityp) + aldif
            enddo
          endif
c         Also compute matrix for least squares fit
          do 1276 ial=1,nalfxi
            do 1275 jal=1,nalfxi
              alsq(ial,jal,1,ityp) = zero
              alsq(ial,jal,2,ityp) = zero
              alsq(ial,jal,3,ityp) = zero
              als = alfxc(ial,ityp) + alfxc(jal,ityp)
              do 1274 nr=1,nrxcfit(ityp)
                r = radmsh(nr,ityp)
                r2 = r*r
                r4 = r2*r2
                expfac = radwt(nr,ityp)*EXP( -als*r2 )
                alsq(ial,jal,1,ityp) = alsq(ial,jal,1,ityp) + expfac
                alsq(ial,jal,2,ityp) = alsq(ial,jal,2,ityp) + expfac*r2
                alsq(ial,jal,3,ityp) = alsq(ial,jal,3,ityp) + expfac*r4
 1274         continue
 1275       continue
 1276     continue
clib: disable calls to linpack/dsifa -
          call STOPXERR( 'DISABLED - removed linpack/dsifa calls' )
clib          call DSIFA(alsq(1,1,1,ityp),nalfxcd,nalfxi,ipvt(1,1,ityp),if)
          if( if.ne.0 )then
            write(IWR,*)'xcfit:',nalfxi,(alfxc(ial,ityp),ial=1,nalfxi)
            call STOPXERR( 'vfiterr1/dsifa-1 error' )
          endif
clib          call DSIFA(alsq(1,1,2,ityp),nalfxcd,nalfxi,ipvt(1,2,ityp),if)
          if( if.ne.0 ) call STOPXERR( 'vfiterr2/dsifa-2 error' )
clib          call DSIFA(alsq(1,1,3,ityp),nalfxcd,nalfxi,ipvt(1,3,ityp),if)
          if( if.ne.0 ) call STOPXERR( 'vfiterr3/dsifa-3 error' )
        endif
cxcf:
cxcf:end
cxcf:###############################################################
c
      RETURN
c
c >>>>>> ERROR HANDLING:
c
 1300 continue
      write(IWR,*) '>>>>> Want label >',wantlbl,'<, found >',readlbl,'<'
      call STOPXERR( 'atmdat1- format error in atom file' )
c
c  Data input formats
c
 8000 format(a)
 8020 format(4(i2,1x))
C9025 format(i2,1x,2f12.8)
 8030 format(4d16.8)
 8040 format(i2,4a6)
 8065 format(3x,6f12.8)
 8070 format(i2,1x,6f12.8)
C8120 format(i2,1x,i2,i3)
 8130 format(i2,a1,3x,2f12.8)
c
c  Data output formats
c
 9000 format(1x,a)
 9020 format(1x,4(i2,1x))
 9025 format(1x,i2,1x,2f12.8)
 9030 format(1x,1p4d16.8)
 9040 format(1x,i2,4a6)
 9065 format(4x,6f12.8)
 9070 format(1x,i2,1x,6f12.8 / (4x,6f12.8) )
 9120 format(1x,i3,1x,i3,1x,i3)
 9130 format(1x,i2,a1,3x,2f12.8)
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
