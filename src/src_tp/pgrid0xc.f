c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PGRID0XC
c
c
      subroutine PGRID0XC( IWR, do_gga, xcnrat, valatm,
     $ xcalcut,xcfac,xcrhocut,
     $ rhotot,cortot,exc0grd,exc0gsl,exc0pt,exc0wt,exccgrd, rhoatom,
     $ natm,ntyp,nshld,nald, nlat  , nptr,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,  nrd,nrcor,radmsh,corden, y2, u,  ratm,rlat,
     $ hh, weight, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ mxwkr, rv, wkr,
     $ rho0,vxc0,exc0,rhoc, drho )
c---------------------------------------------------------------
c Purpose: compute, on grid, combined isolated atom densities
c      and cutoff (short-range) xc-potentials as 0th reference
c      for iteration phase.  Also compute sum of isolated atom
c      xc-energies, integrated on grid.
c
c Written: Peter A. Schultz, 16-June-1997, for v2.20 (boxes)
c
c Revision history:
c  20Jun07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  17Jul06-PAS/2.59a: bug:"nproc"->"nprocl"
c  15Oct04-PAS/2.59-tp: task-parallelism
c  15Oct04-PAS/2.59: split grid-compute and summary
c   8Nov04-PAS/2.58a: cosmetics to match eslofrc
c  11Jun03-PAS/2.55: change DFTGGA arguments
c   8Apr02-PAS/2.53: gga stress
c   7Mar02-PAS/2.52: explicit constants to d0
c  26Jul01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs
c  24Mar01-PAS/2.46a: (gga+cntrch) bugfix
c  15Nov00-PAS/2.46: GGA; nearby split off to "potnr0"
c   4Jan00-PAS/2.42: cosmetic clean-up
c  31Aug99-DBR/2.38: bug fix for cross-shell occupancy cases
c  17Dec98-PAS/2.30: clean force/stress arrays; move force terms
c                    out of setup into new force rtn "eslofrc"
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  16Sep98-PAS/2.27: pseudo-ionic (for xc only) atoms
c  14Aug98-PAS/2.24: cutoffs
c  15Jul98-PAS/2.23: split basis
c  14May98-PAS/2.22: diagnostics/output altered
c  16Jun97-PAS/2.20: rewritten for boxes
c   3Jun97-PAS/2.17: naming conventions, new grid coord rtn,
c                    new exc correct, redo output yet again
c   1May96-PAS/2.15: fix bug in phantom force l=0
c  15Dec95-PAS/2.14: shell-by-shell grid orb ranges
c  28Aug95-PAS/2.13: partial core correction installed
c   7Nov94-PAS/2.11: allow fully non-orthogonal grid vectors
c  15Jun94-PAS/2.10: new cutoff function, install atom Exc
c  31Mar94-PAS: compute ghost xc grid force on isolated atoms
c  vers 2.00: complete rewrite, rename.
c   -installed sterner lattice cutoffs, nearby weighted sample
c---------------------------------------------------------------
c
c Notes:
c  The use of cutoffs in this routine should be compatible
c  with that used in calculating the grid wave functions.
c
c  *** NB **** "grid0xc" is master routine that the routines
c  "potnr0" (to compute nearby potentials for reference atoms),
c  and "eslofrc" must exactly copy in part.  Changes to this
c  this routine should only be done in tandem with corresponding
c  changes in these other routines as well.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
c Input arrays:
      LOGICAL    do_gga
      DIMENSION  xcnrat(4,natm)
      DIMENSION  itypa(natm)
      DIMENSION  numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $  nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  znuc(ntyp)
      DIMENSION  occsh(nshld,ntyp)
      DIMENSION  nrcor(ntyp),radmsh(nrd,*),corden(nrd,*)
      DIMENSION  ratm(3,natm),rlat(3,*), hh(3,3)
      DIMENSION  ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
c Output arrays:
      DIMENSION  rhoatom(natm,*)
      DIMENSION  rho0(*),vxc0(*),exc0(*),rhoc(*), drho(nptr,*)
c Scratch arrays:
      PARAMETER  ( nvalatm = 11 )
      DIMENSION  valatm(nvalatm,natm)
      DIMENSION  y2(*),u(*)
      DIMENSION  rv(3,*), wkr(*)
c
c Local declarations:
      DIMENSION  ratom(3), r0box(3),boxctr(3)
      LOGICAL  poorho, cntrch, anycntr
      DATA  gdrho / 0.001d0 /, badrho / 0.05d0 /, cutocc / 0.9d-8 /
      DATA  alcutfac / 1.33333333d0 /
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  three / 3.d0 /
c flag for additional output 0-none/1-some/2-more/3-full/
      DATA  lstout  / 0 /
c
c Statement function:
      LOGICAL    EVEN
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      fourpi = pi + pi + pi + pi
      argmx = half*cutgaus
c
c Acquire information about parallelism:
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( node0id )
      call MPCOMM( icomm )
c     nprocl: number of processors in this local image
      nprocl = nprocs
c  iprocl0: master/root node of set of processors in this image
      iprocl0 = node0id
      iprocl = iproc - iprocl0
c
c  Initialize sums of atom densities and xc potentials
c
      call MKZERO( 2*natm, rhoatom )
      call MKZERO( nvalatm*natm, valatm )
      rhotot = zero
      cortot = zero
      exc0grd = zero
      exc0gsl = zero
      exccgrd = zero
      exc0pt = zero
      exc0wt = zero
c Clear out derivatives in the event of a GGA calculation:
      if( do_gga )  call MKZERO( 9*nptr, drho )
c
      ipoor = 0
      ibad = 0
      if( lstout.gt.0 )  write(IWR,9010)
 9010 format(/1x,'GRID0XC: analysis on grid of isolated atoms')
c
      nlatmx = 1
c
c >>>>>>>>>> Cycle over atoms:
c
      anycntr = .false.
c  Initialize bitmap index for shell/box occupancy:
      irec00 = 0
      inbox00 = 0
      iatshl00 = 0
c
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
c
        nshli = nocc(ityp)
c
c       Initialize shell/box bitmap index for this atom ...
        irec0 = irec00
        inbox0 = inbox00
        iatshl0 = iatshl00
c        ... and prepare index for the *next* atom:
        do  ishl=1,numshl(ityp)
          inbox00 = 2*inbox00
          if( inrec*(iatshl00/inrec).eq.iatshl00 )then
            irec00 = irec00 + 1
            inbox00 = 1
          endif
          iatshl00 = iatshl00 + 1
        enddo
c
        noccsh = 0
        occatxc = zero
        occatm = zero
        almina = 9.d9
        do 38 ishl=1,nshli
          lshl = locc(ishl,ityp)
          occshl = occsh(ishl,ityp)
          if( occshl.eq.zero ) goto 38
c
          occatm = occatm + occshl
          if( occshl.gt.zero ) occatxc = occatxc + occshl
          noccsh = noccsh + 1
c
          almins = alsh(1,ishl,ityp)
          do  ial=2,nalsh(ishl,ityp)
            ali = alsh(ial,ishl,ityp)
            if( ali.lt.almins ) almins = ali
          enddo
          if( almins.lt.almina ) almina = almins
   38   continue
c
        if( noccsh.eq.0 ) goto 1000
c
        rhoatm = zero
        rhoatxc = zero
        coratm = zero
        vxcatm = zero
        vxcfst = zero
        vxcasl = zero
        excatm = zero
        excfst = zero
        excasl = zero
        exccor = zero
        corfst = zero
c
c       Bring out xc slow/fast function parameters
c
        alminxc = two*almina
        if( xcalcut.lt.zero ) alminxc = -xcalcut*almina
        alcutxc = alcutfac*alminxc
        if( xcalcut.gt.zero ) alcutxc = xcalcut
c
        nrsc = nrcor(ityp) + 1
        if( nrsc.ne.1 )then
c         Interpolate core density
          ddy0 = 1.d31
          ddy1 = 1.d31
          call SPLINE(nrsc,radmsh(1,ityp),corden(1,ityp),ddy0,ddy1,y2,u)
        endif
c
        rngatm = SQRT( cutgrid / (almina+almina) )
c
c       Cycle over boxes:
c
c       Set up task-parallel decomposition over boxes.  Rather than
c       doing MOD( ibox, nprocs ), will do block-wise distribution
c       to simplify gather of big grid vectors.
        nboxpr = nbox / nprocl
        nprocbig = nbox - nprocl*nboxpr
        if( iprocl .lt. nprocbig ) nboxpr = nboxpr + 1
        ibox1 = 1 + iprocl*nboxpr
        if( iprocl .ge. nprocbig ) ibox1 = ibox1 + nprocbig
        iboxn = ibox1 + nboxpr - 1
        ir1proc = 1
        nrproc = 0
c
        irpot = 0
ctp        do 300 ibox=1,nbox
        do 300 ibox=ibox1,iboxn
c
c         Get basic box info
          call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $     ira0,irb0,irc0, idra,idrb,idrc )
          nrbox = idra*idrb*idrc
c
          iv0 = irpot
          irpot = irpot + nrbox
c
          nrproc = nrproc + nrbox
c
          nobox = noboxs(ibox)
          if( nobox.le.0 ) goto 300
c
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          call MKRBOX( hh,r0box, idra,idrb,idrc, rv )
c
c         Set up for scratch space below
          ir0 = 0
          irsq  = ir0
          irngr = ir0 +  1*nrbox
          iatrh = ir0 +  2*nrbox
          icnrh = ir0 +  3*nrbox
          icore = ir0 +  4*nrbox
          idcor = ir0 +  5*nrbox
          ia0e  = icore
          inext = idcor + nrbox
          ia1e  = idcor
          ia2e  = ir0 +  6*nrbox
          idrho = ir0 +  7*nrbox
          iddrh = ir0 +  8*nrbox
          idcrho = ir0 +  9*nrbox
          iddcrh = ir0 + 10*nrbox
          inext = iddcrh + nrbox
c         Check memory sufficiency:
          if( inext.ge.mxwkr ) call STOPXERR( 'xc0-wkrbx' )
c
c         Lattice loop:
c
          do 900 ilat=1,nlat
c
            rngbox = zero
            do  jr=1,3
              ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
              rngbox = rngbox + (ratom(jr)-boxctr(jr))**2
            enddo
            rngbox = SQRT( rngbox ) - boxrad
            if( rngbox.gt.rngatm ) goto 900
c
            rsqbox = rngbox**2
            if( rngbox.lt.zero ) rsqbox = zero
c
            if( ilat.gt.nlatmx ) nlatmx = ilat
c
c           Generate r**2 on mesh ...
            do  nr=1,nrbox
              wkr(irsq+nr) = ( rv(1,nr) - ratom(1) )**2
     $                     + ( rv(2,nr) - ratom(2) )**2
     $                     + ( rv(3,nr) - ratom(3) )**2
              wkr(iatrh+nr) = zero
              wkr(icnrh+nr) = zero
            enddo
c
            do  nr=1,nrbox
              wkr(irngr+nr) = SQRT( wkr(irsq+nr) )
            enddo
c
            if( do_gga )then
c             Start with having no pseudo-ionic countercharges
              anycntr = .false.
c             Clear atom (xc) and global (es) reference density derivatives
              do  nr=1,nrbox
                wkr(idrho+nr) = zero
                wkr(iddrh+nr) = zero
                wkr(idcrho+nr) = zero
                wkr(iddcrh+nr) = zero
              enddo
            endif
c
c           Shell loop:
c
            irec = irec0
            inbox = inbox0
            iatshl = iatshl0
            do 800 ishl=1,nshli
              lshl = locc(ishl,ityp)
              twol = DBLE( 2*lshl )
c
              occshl = occsh(ishl,ityp)
              if( occshl.eq.zero ) goto 800
c
              cntrch = .false.
              if( occshl.lt.zero ) cntrch = .true.
              if( cntrch ) anycntr = .true.
c
c             Determine if shell occupies box.
c             NB: Box occupancy was set using *basis* shells, not the
c             occupied shells, and must convert to shell occupancy.
c
c             Determine min alpha for this *occupancy* shell
              almins = alsh(1,ishl,ityp)
              do  ial=2,nalsh(ishl,ityp)
                ali = alsh(ial,ishl,ityp)
                if( ali.lt.almins ) almins = ali
              enddo
c
c             Cycle over *basis* shells here, and try to find the
c              basis shell that best matches the occ shell above:
              irec = irec0
              inbox = inbox0
              iatshl = iatshl0
              keybox = 0
              alminl = zero
              do 58 jshl=1,numshl(ityp)
                inbox = 2*inbox
                if( inrec*(iatshl/inrec).eq.iatshl )then
                  irec = irec + 1
                  inbox = 1
                endif
                iatshl = iatshl + 1
                if( lshel(jshl,ityp).eq.lshl )then
c                 Determine min alpha for this *basis* shell
                  alminj = ala(1,jshl,ityp)
                  do  jal=1,nala(jshl,ityp)
                    alj = ala(jal,jshl,ityp)
                    if( alj.lt.alminj ) alminj = alj
                  enddo
                  if( alminj.eq.almins )then
c                   Perfect match, use this box occupancy
                    keybox = inboxs(irec,ibox) / inbox
                    goto 59
                  elseif( alminj.le.almins .and. alminj.gt.alminl )then
c                   Better match, provisionally load this box occupancy
                    keybox = inboxs(irec,ibox) / inbox
cpas: as of 2.58 this was unfortunately left out ... (inserted 2.60)
                    alminl = alminj
                  endif
                endif
c
   58         continue
   59         continue
c
              if( EVEN(keybox) ) goto 800
c
              nali = nalsh(ishl,ityp)
c
              almins = alsh(1,ishl,ityp)
              do  ial=2,nali
                ali = alsh(ial,ishl,ityp)
                if( ali.lt.almins ) almins = ali
              enddo
              almin = almina
              rsqshl = cutgrid / (almin +almins)
              if( rsqbox.gt.rsqshl ) goto 800
c
              call MKZERO( nrbox, wkr(ia0e+1) )
              if( do_gga )then
                do  nr=1,nrbox
C                 wkr(ia0e+nr) = zero
                  wkr(ia1e+nr) = zero
                  wkr(ia2e+nr) = zero
                enddo
              endif
c
              do 190 ial=1,nali
                ali = alsh(ial,ishl,ityp)
                cali = occshl*calsh(ial,ishl,ityp)
                cali0 = cali
c               Just do off-diags once, for Pete's sake:
                do 180 jal=ial,1,-1
                  alj = alsh(jal,ishl,ityp)
                  alsum = ali + alj
                  als2 = alsum*alsum
                  twoals = alsum + alsum
                  coef = cali*calsh(jal,ishl,ityp)/fourpi
c                 Double the off-diags:
                  cali = cali0 + cali0
c                 Looser cutoff used to get quicker convergence
                  algt = ali
                  if( alj.gt.ali ) algt = alj
                  twoalgt = algt + algt
                  twoalc = algt + almin
                  rsqcut = cutgrid/twoalc
c
c                 Skip if alpha is too large
                  if( rsqbox.gt.rsqcut ) goto 180
c
c                 Sums: occ*Ca*Cb*(a+b)**n*exp(-(a+b)r**2); n=0,1,2
                  if( do_gga )then
                    do  nr=1,nrbox
                      gauss = coef*EXP( -alsum*wkr(irsq+nr) )
                      wkr(ia0e+nr) = wkr(ia0e+nr) + gauss
                      wkr(ia1e+nr) = wkr(ia1e+nr) + gauss*alsum
                      wkr(ia2e+nr) = wkr(ia2e+nr) + gauss*als2
                    enddo
                  else
                    do  nr=1,nrbox
                      gauss = coef*EXP( -alsum*wkr(irsq+nr) )
                      wkr(ia0e+nr) = wkr(ia0e+nr) + gauss
                    enddo
                  endif
c
c                 Close alpha loops:
c
  180           continue
  190         continue
c
              if( do_gga )then
c               Assemble density and derivatives for GGA potential:
c
                if( cntrch )then
c                  ... for global (es) reference density
                  iprho  = icnrh
                  ipdrho = idcrho
                  ipddrh = iddcrh
                else
c                  ... atomic (pseudo-ionic) xc reference density
                  iprho  = iatrh
                  ipdrho = idrho
                  ipddrh = iddrh
                endif
c
                call MKSHL2D( lshl, nrbox, wkr(irngr+1), wkr(irsq+1),
     $           wkr(ia0e+1),wkr(ia1e+1),wkr(ia2e+1),
     $           wkr(iprho+1),wkr(ipdrho+1),wkr(ipddrh+1) )
c
              else
c               Just need to assemble density ...
c
                if( cntrch )then
c                  ... to update compensating countercharge for es
                  iprho = icnrh
                else
c                  ... for LDA calculation (no density derivs needed)
                  iprho = iatrh
                endif
c
                call MKSHL( lshl, nrbox, wkr(irsq+1),
     $           wkr(ia0e+1), wkr(iprho+1) )
c
              endif
c
c             Close shell loop:
c
  800       continue
c
c           Check if we have core density to consider
c
            call MKZERO( nrbox, wkr(icore+1) )
            if( nrsc.ne.1 )then
              rsqcut = radmsh(nrsc,ityp)**2
              if( rsqbox.gt.rsqcut ) goto 195
c
              if( do_gga )then
c               Do core 2nd deriv before 1st deriv, and use (idcor)
c               space for scratch  (do not need core 2nd deriv
c               separate from the total rho 2nd deriv).
c               Interpolate to get 2nd deriv. of core density ...
                call SPLIND2( nrsc, radmsh(1,ityp), corden(1,ityp),y2,
     $                        nrbox,wkr(irngr+1),wkr(idcor+1) )
c                ... and put into total 2nd deriv of density for GGA
                do  nr=1,nrbox
                  wkr(iddrh+nr) = wkr(iddrh+nr) + wkr(idcor+nr)
                enddo
c
c               Interpolate to get deriv. of core density ...
                call SPLIND( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                       nrbox,wkr(irngr+1),wkr(idcor+1) )
c                ... and put into density derivative for GGA potential
                do  nr=1,nrbox
                  wkr(idrho+nr) = wkr(idrho+nr) + wkr(idcor+nr)
                enddo
              endif
c
c             Interpolate to get core density on grid ...
              call SPLINT( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                     nrbox,wkr(irngr+1),wkr(icore+1) )
c
            endif
  195       continue
c
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * *
c             Have total isolated-atom density on this part
c                 of the grid; compute xc-quantities.
c         * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
            do 200 nr=1,nrbox
              rhoxc = wkr(iatrh+nr)
              rhoes = wkr(iatrh+nr) + wkr(icnrh+nr)
              corho = wkr(icore+nr)
              rhopt = rhoxc + corho
c
              if( rhopt.gt.xcrhocut )then
c
                if( do_gga )then
c                 Assemble gradient terms needed for GGA
c
                  rmag = wkr(irngr+nr)
                  if( rmag.gt.1.d-6 )then
c                   For r>0, general Cartesian derivs from radial
c                    ... note dcor already added into drho above
                    drh = wkr(idrho+nr)
c                    ... terms for Cartesian reference derivatives:
                    d1rhdr = drh/rmag
                    d2rh = wkr(iddrh+nr) - d1rhdr
c
c                     ... need gradients for the atomic GGA potential
                    grho = ABS( drh )
                    g2rho = wkr(iddrh+nr) + two*d1rhdr
                    gpggp = grho*wkr(iddrh+nr)
c
                    if( anycntr )then
c                     Need to produce derivatives of (es) reference density
c                     rather than pseudo-ionic reference density used for xc
                      drh = wkr(idrho+nr) + wkr(idcrho+nr)
                      d1rhdr = drh/rmag
                      d2rh = wkr(iddrh+nr) + wkr(iddcrh+nr) - d1rhdr
                    endif
c
                    ux = ( rv(1,nr) - ratom(1) ) / rmag
                    uy = ( rv(2,nr) - ratom(2) ) / rmag
                    uz = ( rv(3,nr) - ratom(3) ) / rmag
                    ir = iv0 + nr
c                    ... first derivs: X,Y,Z
                    drho(ir,1) = drho(ir,1) + ux*drh
                    drho(ir,2) = drho(ir,2) + uy*drh
                    drho(ir,3) = drho(ir,3) + uz*drh
c                    ... second derivs: XX,YY,ZZ, XY,YZ,ZX
                    drho(ir,4) = drho(ir,4) + ux*ux*d2rh + d1rhdr
                    drho(ir,5) = drho(ir,5) + uy*uy*d2rh + d1rhdr
                    drho(ir,6) = drho(ir,6) + uz*uz*d2rh + d1rhdr
                    drho(ir,7) = drho(ir,7) + ux*uy*d2rh
                    drho(ir,8) = drho(ir,8) + uy*uz*d2rh
                    drho(ir,9) = drho(ir,9) + uz*ux*d2rh
                  else
c                   At r=0, dp/du=0, and only d2p/du2 survives limits ...
                    d2rh = wkr(iddrh+nr)
c
                    grho = zero
                    g2rho = three*wkr(iddrh+nr)
                    gpggp = zero
c
                    if( anycntr ) d2rh = wkr(iddrh+nr) + wkr(iddcrh+nr)
c
                    ir = iv0 + nr
                    drho(ir,4) = drho(ir,4) + d2rh
                    drho(ir,5) = drho(ir,5) + d2rh
                    drho(ir,6) = drho(ir,6) + d2rh
                  endif
c
                  call DFTGGA( rhopt, grho,g2rho,gpggp,
     $                         vxcpt,excpt, dxcpt,dxc2pt )
c
c                 End do_gga = .true.
                else
c                 do_gga = .false.
c
                  call DFTLDA( rhopt, vxcpt,excpt )
c
                endif
c
c               Add atom density to sum of atom densities
                rho0(iv0+nr) = rho0(iv0+nr) + rhoes
                rhoc(iv0+nr) = rhoc(iv0+nr) + corho
c
c               Construct "fast" xc-potentials, and update:
                cfac = EXP( -alcutxc*wkr(irsq+nr) )
                fstfcn = cfac / ( (one-xcfac)*cfac + xcfac )
c
                vfast = vxcpt*fstfcn
                vslow = vxcpt - vfast
c
                efast = excpt*fstfcn
                eslow = excpt - efast
c
                vxc0(iv0+nr) = vxc0(iv0+nr) + vfast
                exc0(iv0+nr) = exc0(iv0+nr) + efast
c
c               Integrate atomic electron density and xc-energy on grid
                rhoatm = rhoatm + rhoes
                rhoatxc = rhoatxc + rhoxc
                coratm = coratm + corho
                vxcatm = vxcatm + vxcpt*rhoxc
                vxcfst = vxcfst + vfast*rhoxc
                vxcasl = vxcasl + vslow*rhoxc
                excatm = excatm + excpt*rhoxc
                excfst = excfst + efast*rhoxc
                excasl = excasl + eslow*rhoxc
                exccor = exccor + excpt*corho
                corfst = corfst + efast*corho
c
              endif
c
  200       continue
c
c           Close lattice loop:
c
  900     continue
c
c         Close box loop:
c
  300   continue
c
c       Save isolated atom grid integral quantities:
c
        valatm( 1,iatm) = rhoatm
        valatm( 2,iatm) = rhoatxc
        valatm( 3,iatm) = coratm
        valatm( 4,iatm) = vxcatm
        valatm( 5,iatm) = vxcfst
        valatm( 6,iatm) = vxcasl
        valatm( 7,iatm) = excatm
        valatm( 8,iatm) = excfst
        valatm( 9,iatm) = excasl
        valatm(10,iatm) = exccor
        valatm(11,iatm) = corfst
c
c       Close atom loop:
c
 1000 continue
c
      if( nprocl .gt. 1 )then
c       Reduce all atomic grid integrals to master
c
        nreduc = nvalatm*natm
        call MPREDUC8( iprocl0, nreduc, valatm, wkr, icomm )
C        call MPI_Reduce( valatm, wkr, nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   iprocl0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'grid0xc - Reduce ints' )
        call DCOPY( nreduc, wkr,1, valatm,1 )
c
        call MPBCAST8( iprocl0, nreduc, valatm, icomm )
C        call MPI_Bcast( valatm, nreduc, MPI_DOUBLE_PRECISION,
C     $   iprocl0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'grid0xc - Bcast ints' )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Assemble, analyze, and output final atomic integrals
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do 2000 iatm=1,natm
        ityp = itypa(iatm)
c
        nshli = nocc(ityp)
        noccsh = 0
        occatm = zero
        occatxc = zero
        do 2038 ishl=1,nshli
          lshl = locc(ishl,ityp)
          occshl = occsh(ishl,ityp)
          if( occshl.eq.zero ) goto 2038
c
          occatm = occatm + occshl
          if( occshl.gt.zero ) occatxc = occatxc + occshl
          noccsh = noccsh + 1
 2038   continue
        if( noccsh.eq.0 ) goto 2000
c
c       Retrieve local density/exc-den for weighted nearby sampling
c
        vxcwtnr = xcnrat(1,iatm)
        excwtnr = xcnrat(2,iatm)
        vxcptnr = xcnrat(3,iatm)
        excptnr = xcnrat(4,iatm)
c
c       Atom-resolved grid integral quantities:
c
        rhoatm  = valatm( 1,iatm)
        rhoatxc = valatm( 2,iatm)
        coratm  = valatm( 3,iatm)
        vxcatm  = valatm( 4,iatm)
        vxcfst  = valatm( 5,iatm)
        vxcasl  = valatm( 6,iatm)
        excatm  = valatm( 7,iatm)
        excfst  = valatm( 8,iatm)
        excasl  = valatm( 9,iatm)
        exccor  = valatm(10,iatm)
        corfst  = valatm(11,iatm)
c
        rhoatm = weight*rhoatm
        rhoatxc = weight*rhoatxc
        coratm = weight*coratm
        vxcatm = weight*vxcatm
        vxcfst = weight*vxcfst
        vxcasl = weight*vxcasl
        excatm = weight*excatm
        excfst = weight*excfst
        excasl = weight*excasl
        exccor = weight*exccor
        corfst = weight*corfst
        vxcatnr = (occatxc-rhoatxc)*vxcptnr
        vxcatwt = (occatxc-rhoatxc)*vxcwtnr
        excatnr = (occatxc-rhoatxc)*excptnr
        excatwt = (occatxc-rhoatxc)*excwtnr
c
        rhoatom(iatm,1) = rhoatm
        rhotot = rhotot + rhoatm
        rhoatom(iatm,2) = coratm
        cortot = cortot + coratm
        exc0grd = exc0grd + excatm
        exc0gsl = exc0gsl + excasl
        exccgrd = exccgrd + exccor
        exc0pt = exc0pt + excatnr
        exc0wt = exc0wt + excatwt
c
        grderr = ABS( rhoatm - znuc(ityp) )
        poorho = grderr .gt. gdrho
        if( poorho ) ipoor = ipoor + 1
        if( grderr .gt. badrho ) ibad = ibad + 1
        if( lstout.gt.1 .or. poorho )then
 9021     format(1x,'Atom',i4,'   Z=',f5.1,'   occ=',f6.2,
     $              '   rho=',f14.10,'   core=',f14.10)
 9022     format( 5x, a, 5(1x,f12.8) )
          write(IWR,9021) iatm,znuc(ityp),occatm,rhoatm,coratm
          if( ABS(rhoatxc-rhoatm).gt.cutocc )
     $    write(IWR,9031)                        rhoatxc
 9031     format(1x,'    ',4x,'     ',5x  ,'       ',5x,
     $              'rho(xc)=',f14.10)
        endif
        if( lstout.gt.3 .or. poorho )then
          write(IWR,9022) 'vxc(afsw): ',vxcatm,vxcfst,vxcasl,vxcatwt
          write(IWR,9022) 'exc(afsw): ',excatm,excfst,excasl,excatwt
          if( nrsc.ne.1 )then
            write(IWR,9022) 'core(af): ',exccor,corfst
            write(IWR,9022) 'ectot',excatm+exccor,excfst+corfst,excasl
          endif
        endif
        if( lstout.gt.4 )then
          write(IWR,9022)  'vxcpt(wp):',vxcwtnr,vxcptnr,vxcatwt,vxcatnr
          write(IWR,9022)  'excpt(wp):',excwtnr,excptnr,excatwt,excatnr
        endif
c
c       Close atom loop:
c
 2000 continue
c
c Set maximum lattice index that was used to limit later grid routines:
c
      if( nprocl .gt. 1 )then
c       A hack to communicate max required nlat:
        call MKZERO( nprocl, wkr )
        wkr(1) = DBLE( nlatmx ) + 0.1
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nprocl,1,1, wkr, icomm )
        if( iprocl.eq.iprocl0 )then
          rlatmax = wkr(1)
          do  jpr=2,nprocl
            if( wkr(jpr) .gt. rlatmax ) rlatmax = wkr(jpr)
          enddo
        endif
        call MPBCAST8( iprocl0, 1,rlatmax, icomm )
        nlatmx = rlatmax
      endif
      nlat = nlatmx
c
c Add a grid diagnostic message:
c
      if( ibad.gt.0 )then
        call STOPXERR( 'grid0xc-grid is not fine enough' )
      elseif( ipoor.gt.0 )then
        write(IWR,*) '***** CAUTION: grid fineness questionable'
      endif
c
      if( nprocl .gt. 1 )then
c       Gather the full grid fields to the master/root node:
c
cpas:
        call TIMER( 'grid0xc compute' )
c
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nptr,ir1proc,nrproc, rho0, icomm )
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nptr,ir1proc,nrproc, rhoc, icomm )
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nptr,ir1proc,nrproc, vxc0, icomm )
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nptr,ir1proc,nrproc, exc0, icomm )
        if( do_gga )then
          do  idrh=1,9
            call MPMERGER8( nprocl, iprocl0, iprocl,
     $       nptr,ir1proc,nrproc, drho(1,idrh), icomm)
          enddo
        endif
cpas:
        call TIMER( 'grid0xc Merge  ' )
c
c       For now, broadcast these to all nodes:
        call MPBCAST8( iprocl0, nptr, rho0, icomm )
        call MPBCAST8( iprocl0, nptr, rhoc, icomm )
        call MPBCAST8( iprocl0, nptr, vxc0, icomm )
        call MPBCAST8( iprocl0, nptr, exc0, icomm )
        if( do_gga )then
          do  idrh=1,9
            call MPBCAST8( iprocl0, nptr, drho(1,idrh), icomm )
          enddo
        endif
cpas:
        call TIMER( 'grid0xc Bcast  ' )
c
c       End tp
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
