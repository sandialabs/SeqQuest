c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PESLOFRC
c
c
      subroutine PESLOFRC( IWR, do_gga, ndim,
     $ frcesl,stresl, frcxc0,strxc0, frcxcc,strxcc,
     $ xcalcut,xcfac,xcrhocut,  natm,ntyp,nshld,nald, nlat  ,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,  nrd,nrcor,radmsh,corden, y2, u,  ratm,rlat,
     $ hh, weight, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wkr, frctp,
     $ ves, rv )
c---------------------------------------------------------------
c Purpose: compute grid force term from spherical atom es, and
c          get "phantom" force from grid due to spherical xc
c
c Written: Peter A. Schultz, 17-December-1998, v2.30
c
c Revision history:
c  11May07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  23Jan06-PAS/2.59: align tp/serial
c  24Sep04-PAS:2.59-tp: clean-up task-parallel routines
c   6May03-DMC: Add local copies of the three frc* arrays, an array to
c               hold a local copy of the str* arrays, and work space
c               to hold a second copy of each array.  This space is
c               contiguous, and thus can be distributed across procs
c               with a single send.  All this for parallelization.
c   6Nov04-PAS/2.58a: add dGGA stress analytic-grid correction
c   8Apr02-PAS/2.53: gga stress
c   7Mar02-PAS/2.52: explicit constants to d0.
c  15Nov00-PAS/2.46: GGA
c  18Nov99-PAS/2.42: phantom stress terms for xc
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  23Sep99-PAS/2.39: countercharge bug: fix density derivatives
c  22Sep99-PAS/2.39: simple bugfix (thx Dave Raczkowski)
c---------------------------------------------------------------
c
c Note:
c  This routine modelled off the "grid0xc" routine.
c
c  *** NB *** Changes to this routine should only be done in
c  tandem with changes to "grid0xc" and "potnr0".
c
c  The use of cutoffs in this routine should be compatible
c  with that used in calculating the grid wave functions.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
c Output arrays:
      DIMENSION  frcesl(3,natm), stresl(3,3)
      DIMENSION  frcxc0(3,natm), strxc0(3,3)
      DIMENSION  frcxcc(3,natm), strxcc(3,3)
ctp: parallel scratch:
      DIMENSION  frctp(3,natm,6)
      DIMENSION  strtp(3,3,10)
c Input arrays:
      LOGICAL    do_gga
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
      DIMENSION  ves(*)
c Scratch arrays:
      DIMENSION  y2(*),u(*)
      DIMENSION  rv(3,*), wkr(*)
c
c Local declarations:
      DIMENSION  strdxc0(3,3),strdxcc(3,3)
      DIMENSION  stre(3,3),strx(3,3),strc(3,3), strd0(3,3),strdc(3,3)
      DIMENSION  ratom(3), r0box(3),boxctr(3)
      LOGICAL  poorho, cntrch
      DATA  gdrho / 0.001d0 /, badrho / 0.05d0 /, cutocc / 0.9d-8 /
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  three / 3.d0 /
      DATA  alcutfac / 1.33333333d0 /
c flag for additional output 0-none/1-some/2-full
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
ctp:
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      call MKZERO( 3*natm, frcesl )
      call MKZERO( 9, stresl )
      call MKZERO( 3*natm, frcxc0 )
      call MKZERO( 9, strxc0 )
      call MKZERO( 3*natm, frcxcc )
      call MKZERO( 9, strxcc )
c
      call MKZERO( 9, strdxc0 )
      call MKZERO( 9, strdxcc )
c
      if( lstout.gt.0 )  write(IWR,9010)
 9010 format(/1x,'ESLOFRC: analysis on grid of isolated atom forces')
c
      nlatmx = 1
c
c >>>>>>>>>> Cycle over atoms:
c
c  Initialize bitmap index for shell/box occupancy:
      irec00 = 0
      inbox00 = 0
      iatshl00 = 0
      iorb00 = 0
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
ctp:
c       Each processor will do every nprocs element in this loop.
        if( MOD( iatm, nprocs ) .ne. iproc ) goto 1000
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
        fxat = zero
        fyat = zero
        fzat = zero
        fxc = zero
        fyc = zero
        fzc = zero
        fxes = zero
        fyes = zero
        fzes = zero
        call MKZERO( 9, strx )
        call MKZERO( 9, stre )
        call MKZERO( 9, strc )
        call MKZERO( 9, strd0 )
        call MKZERO( 9, strdc )
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
        irpot = 0
        do 300 ibox=1,nbox
c
c         Get basic box info
          call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $     ira0,irb0,irc0, idra,idrb,idrc )
          nrbox = idra*idrb*idrc
c
          iv0 = irpot
          irpot = irpot + nrbox
c
          nobox = noboxs(ibox)
          if( nobox.le.0 ) goto 300
c
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          call MKRBOX( hh,r0box, idra,idrb,idrc, rv )
c
c         Set up for scratch space below:
          ir0 = 0
          irsq  = ir0
          irngr = ir0 +  1*nrbox
          iatrh = ir0 +  2*nrbox
          icnrh = ir0 +  3*nrbox
          idrho = ir0 +  4*nrbox
          icndr = ir0 +  5*nrbox
          icore = ir0 +  6*nrbox
          idcor = ir0 +  7*nrbox
          ia0e  = icore
          ia1e  = idcor
          ia2e  = ir0 +  8*nrbox
          iddrh = ir0 +  9*nrbox
          inext = iddrh + nrbox
c         I neglect this memory check as I am in BIG workspace
c          if( inext.gt.mxwkrd )
c     $     call STOPXERR( 'ESLOFRC: not enough scratch memory' )
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
              wkr(irngr+nr) = SQRT( wkr(irsq+nr) )
              wkr(iatrh+nr) = zero
              wkr(icnrh+nr) = zero
              wkr(idrho+nr) = zero
              wkr(icndr+nr) = zero
              wkr(iddrh+nr) = zero
            enddo
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
              do  nr=1,nrbox
                wkr(ia0e+nr) = zero
                wkr(ia1e+nr) = zero
                wkr(ia2e+nr) = zero
              enddo
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
c                 Skip if alphas too large
                  if( rsqbox.gt.rsqcut ) goto 180
c
c                 Sums: occ*Ca*Cb*(a+b)**n*exp(-(a+b)r**2); n=0,1,2
                  do  nr=1,nrbox
                    gauss = coef*EXP( -alsum*wkr(irsq+nr) )
                    wkr(ia0e+nr) = wkr(ia0e+nr) + gauss
                    wkr(ia1e+nr) = wkr(ia1e+nr) + gauss*alsum
                    wkr(ia2e+nr) = wkr(ia2e+nr) + gauss*als2
                  enddo
c
c                 Close alpha loops:
c
  180           continue
  190         continue
c
              if( do_gga .and. .not. cntrch )then
c               Assemble density and derivatives for GGA potential:
c
                call MKSHL2D( lshl, nrbox, wkr(irngr+1), wkr(irsq+1),
     $           wkr(ia0e+1),wkr(ia1e+1),wkr(ia2e+1),
     $           wkr(iatrh+1),wkr(idrho+1),wkr(iddrh+1) )
c
              else
c               Just need to assemble density and first derivatives ...
c
                if( cntrch )then
c                  ... to update compensating countercharge for es
                  iprho = icnrh
                  ipdrh = icndr
                else
c                  ... for LDA calculation (no second derivs needed)
                  iprho = iatrh
                  ipdrh = idrho
                endif
c
                call MKSHL1D( lshl, nrbox, wkr(irngr+1), wkr(irsq+1),
     $           wkr(ia0e+1),wkr(ia1e+1),
     $           wkr(iprho+1),wkr(ipdrh+1) )
c
              endif
c
c             Close shell loop:
c
  800       continue
c
c           Now check if we have core density to consider
c
            call MKZERO( nrbox, wkr(icore+1) )
            call MKZERO( nrbox, wkr(idcor+1) )
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
              endif
c
c             Interpolate to get deriv. of core density ...
              call SPLIND( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                     nrbox,wkr(irngr+1),wkr(idcor+1) )
c
c             Interpolate to get core density on grid ...
              call SPLINT( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                     nrbox,wkr(irngr+1),wkr(icore+1) )
c
            endif
  195       continue
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * *
c             Have total isolated-atom density on this part
c                 of the grid; compute xc-quantities.
c         * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
            do 200 nr=1,nrbox
c
              rhoxc = wkr(iatrh+nr)
              rhoes = wkr(iatrh+nr) + wkr(icnrh+nr)
              corho = wkr(icore+nr)
              rhopt = rhoxc + corho
c
              if( rhopt.gt.xcrhocut )then
c
                if( do_gga )then
c
c                 Assemble gradient terms needed for GGA
c
                  rmag = wkr(irngr+nr)
                  if( rmag.gt.1.d-6 )then
c                   For r>0, general Cartesian derivs from radial
c                    ... note dcor already added into drho above
                    drh = wkr(idrho+nr) + wkr(idcor+nr)
c                    ... terms for Cartesian reference derivatives:
                    d1rhdr = drh/rmag
                    d2rh = wkr(iddrh+nr) - d1rhdr
c                    ... need gradients for atomic GGA potential
                    grho = ABS( drh )
                    g2rho = wkr(iddrh+nr) + two*d1rhdr
                    gpggp = grho*wkr(iddrh+nr)
                  else
c                   At r=0, dp/du=0, and only d2p/du2 survives limits ...
                    d2rh = wkr(iddrh+nr)
                    grho = zero
                    g2rho = three*wkr(iddrh+nr)
                    gpggp = zero
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
c               Construct "fast" xc-potentials, and update:
                cfac = EXP( -alcutxc*wkr(irsq+nr) )
                fstfcn = cfac / ( (one-xcfac)*cfac + xcfac )
c
                vfast = vxcpt*fstfcn
                vslow = vxcpt - vfast
                if( do_gga )then
                  dxcfast = fstfcn*dxcpt
                  dxc2fast = fstfcn*dxc2pt
                endif
c
                efast = excpt*fstfcn
                eslow = excpt - efast
c
c               Integrate atomic electron density and xc quantities
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
                if( wkr(irngr+nr).gt.1.d-4 )then
                  dx = rv(1,nr) - ratom(1)
                  dy = rv(2,nr) - ratom(2)
                  dz = rv(3,nr) - ratom(3)
c
c                 Force due to slo grid es potential
                  drhoes = wkr(idrho+nr) + wkr(icndr+nr)
                  vesrpt = drhoes*ves(iv0+nr)/wkr(irngr+nr)
                  fxes = fxes + vesrpt*dx
                  fyes = fyes + vesrpt*dy
                  fzes = fzes + vesrpt*dz
                  stre(1,1) = stre(1,1) + vesrpt*dx*dx
                  stre(2,1) = stre(2,1) + vesrpt*dy*dx
                  stre(3,1) = stre(3,1) + vesrpt*dz*dx
                  stre(1,2) = stre(1,2) + vesrpt*dx*dy
                  stre(2,2) = stre(2,2) + vesrpt*dy*dy
                  stre(3,2) = stre(3,2) + vesrpt*dz*dy
                  stre(1,3) = stre(1,3) + vesrpt*dx*dz
                  stre(2,3) = stre(2,3) + vesrpt*dy*dz
                  stre(3,3) = stre(3,3) + vesrpt*dz*dz
c
c                 Phantom xc force grid term
c                 Note: d(rho)/dx = (d(rho)/dr)*(dr/dx); dr/dx = x/r
                  drhoxc = wkr(idrho+nr)
                  vxcrpt = drhoxc*vfast/wkr(irngr+nr)
                  fxat = fxat + vxcrpt*dx
                  fyat = fyat + vxcrpt*dy
                  fzat = fzat + vxcrpt*dz
                  strx(1,1) = strx(1,1) + vxcrpt*dx*dx
                  strx(2,1) = strx(2,1) + vxcrpt*dy*dx
                  strx(3,1) = strx(3,1) + vxcrpt*dz*dx
                  strx(1,2) = strx(1,2) + vxcrpt*dx*dy
                  strx(2,2) = strx(2,2) + vxcrpt*dy*dy
                  strx(3,2) = strx(3,2) + vxcrpt*dz*dy
                  strx(1,3) = strx(1,3) + vxcrpt*dx*dz
                  strx(2,3) = strx(2,3) + vxcrpt*dy*dz
                  strx(3,3) = strx(3,3) + vxcrpt*dz*dz
                  if( do_gga )then
c                   Phantom GGA stress grid term
                    dxcrpt = ( dxcfast * wkr(idrho+nr)**2 +
     $                         dxc2fast* wkr(iddrh+nr) )
     $                       / wkr(irsq+nr)
                    strd0(1,1) = strd0(1,1) - dxcrpt*dx*dx
                    strd0(2,1) = strd0(2,1) - dxcrpt*dy*dx
                    strd0(3,1) = strd0(3,1) - dxcrpt*dz*dx
                    strd0(1,2) = strd0(1,2) - dxcrpt*dx*dy
                    strd0(2,2) = strd0(2,2) - dxcrpt*dy*dy
                    strd0(3,2) = strd0(3,2) - dxcrpt*dz*dy
                    strd0(1,3) = strd0(1,3) - dxcrpt*dx*dz
                    strd0(2,3) = strd0(2,3) - dxcrpt*dy*dz
                    strd0(3,3) = strd0(3,3) - dxcrpt*dz*dz
                  endif
                  if( nrsc.ne.1 )then
c                   Phantom xc force grid term (core)
                    vxcrpt = wkr(idcor+nr)*vfast/wkr(irngr+nr)
                    fxc = fxc + vxcrpt*dx
                    fyc = fyc + vxcrpt*dy
                    fzc = fzc + vxcrpt*dz
                    strc(1,1) = strc(1,1) + vxcrpt*dx*dx
                    strc(2,1) = strc(2,1) + vxcrpt*dy*dx
                    strc(3,1) = strc(3,1) + vxcrpt*dz*dx
                    strc(1,2) = strc(1,2) + vxcrpt*dx*dy
                    strc(2,2) = strc(2,2) + vxcrpt*dy*dy
                    strc(3,2) = strc(3,2) + vxcrpt*dz*dy
                    strc(1,3) = strc(1,3) + vxcrpt*dx*dz
                    strc(2,3) = strc(2,3) + vxcrpt*dy*dz
                    strc(3,3) = strc(3,3) + vxcrpt*dz*dz
                    if( do_gga )then
c                     Phantom GGA stress grid term (core)
                      dxcrpt = ( dxcfast * wkr(idcor+nr)**2 )
     $                         / wkr(irsq+nr)
c                     Note that d2cor part inside ddrh above
                      strdc(1,1) = strdc(1,1) - dxcrpt*dx*dx
                      strdc(2,1) = strdc(2,1) - dxcrpt*dy*dx
                      strdc(3,1) = strdc(3,1) - dxcrpt*dz*dx
                      strdc(1,2) = strdc(1,2) - dxcrpt*dx*dy
                      strdc(2,2) = strdc(2,2) - dxcrpt*dy*dy
                      strdc(3,2) = strdc(3,2) - dxcrpt*dz*dy
                      strdc(1,3) = strdc(1,3) - dxcrpt*dx*dz
                      strdc(2,3) = strdc(2,3) - dxcrpt*dy*dz
                      strdc(3,3) = strdc(3,3) - dxcrpt*dz*dz
                    endif
                  endif
                endif
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
        fxes = weight*fxes
        fyes = weight*fyes
        fzes = weight*fzes
        frcesl(1,iatm) = -fxes
        frcesl(2,iatm) = -fyes
        frcesl(3,iatm) = -fzes
        fxat = weight*fxat
        fyat = weight*fyat
        fzat = weight*fzat
        frcxc0(1,iatm) = -fxat
        frcxc0(2,iatm) = -fyat
        frcxc0(3,iatm) = -fzat
        fxc = weight*fxc
        fyc = weight*fyc
        fzc = weight*fzc
        frcxcc(1,iatm) = -fxc
        frcxcc(2,iatm) = -fyc
        frcxcc(3,iatm) = -fzc
c
        do  jd=1,ndim
          do  id=1,ndim
            stresl(id,jd) = stresl(id,jd) - weight*stre(id,jd)
            strxc0(id,jd) = strxc0(id,jd) - weight*strx(id,jd)
            strxcc(id,jd) = strxcc(id,jd) - weight*strc(id,jd)
            strdxc0(id,jd) = strdxc0(id,jd) - weight*strd0(id,jd)
            strdxcc(id,jd) = strdxcc(id,jd) - weight*strdc(id,jd)
          enddo
        enddo
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
c
        grderr = ABS( rhoatm - znuc(ityp) )
        poorho = grderr .gt. gdrho
        if( lstout.gt.0 .or. poorho )then
 9021     format(1x,'Atom',i4,'   Z=',f5.1,'   occ=',f6.2,
     $              '   rho=',f14.10,'   core=',f14.10)
 9022     format( 5x, a, 5(1x,f12.8) )
          write(IWR,9021) iatm,znuc(ityp),occatm,rhoatm,coratm
          if( ABS(rhoatxc-rhoatm).gt.cutocc )
     $    write(IWR,9031)                        rhoatxc
 9031     format(1x,'    ',4x,'     ',5x  ,'       ',5x,
     $              'rho(xc)=',f14.10)
        endif
        if( lstout.gt.1 .or. poorho )then
          write(IWR,9022) 'vxc(afs ): ',vxcatm,vxcfst,vxcasl
          write(IWR,9022) 'exc(afs ): ',excatm,excfst,excasl
          if( nrsc.ne.1 )then
            write(IWR,9022) 'core(af): ',exccor,corfst
            write(IWR,9022) 'ectot',excatm+exccor,excfst+corfst,excasl
          endif
        endif
        if( lstout.gt.3 )then
          write(IWR,9023)'at,xc0frc=',iatm,fxat,fyat,fzat
          write(IWR,9023)'at,corfrc=',iatm,fxc,fyc,fzc
          write(IWR,9023)'at,vesfrc=',iatm,fxes,fyes,fzes
 9023     format(1x,a,i3,3x,3f12.8)
        endif
c
c       Close atom loop:
c
 1000 continue
ctp:start
c
      if( nprocs .gt. 1 )then
c       Each processor has forces/stress from its part if the grid.
c       Need to sum these into a single array, on master (node 0).
c       For now, also broadcast total back to all processors.
c
c       ... load all the force stuff into packed scratch space ...
        do  iatm=1,natm
          do  id=1,3
            frctp(id,iatm,4) = frcesl(id,iatm)
            frctp(id,iatm,5) = frcxc0(id,iatm)
            frctp(id,iatm,6) = frcxcc(id,iatm)
          enddo
        enddo
c
c        ... sum it onto master
        nreduc = 3*3*natm
        call MPREDUC8( master, nreduc,frctp(1,1,4),frctp(1,1,1),icomm )
C        call MPI_Reduce( frctp(1,1,4), frctp(1,1,1), nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'peslofrc - Reduce frc' )
c
c        ... broadcast summed total to everybody
        call MPBCAST8( master, nreduc, frctp , icomm )
C        call MPI_Bcast( frctp, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'peslofrc - Bcast frc' )
c
c        ... and put the final total forces where they belong:
        do  iatm=1,natm
          do  id=1,3
            frcesl(id,iatm) = frctp(id,iatm,1)
            frcxc0(id,iatm) = frctp(id,iatm,2)
            frcxcc(id,iatm) = frctp(id,iatm,3)
          enddo
        enddo
c
c       Repeat for process for the stress:
        do  jd=1,3
          do  id=1,3
            strtp(id,jd, 6) = stresl(id,jd)
            strtp(id,jd, 7) = strxc0(id,jd)
            strtp(id,jd, 8) = strxcc(id,jd)
            strtp(id,jd, 9) = strdxc0(id,jd)
            strtp(id,jd,10) = strdxcc(id,jd)
          enddo
        enddo
        nreduc = 5*9
        call MPREDUC8( master, nreduc,strtp(1,1,6),strtp(1,1,1),icomm )
C        call MPI_Reduce( strtp(1,1,4), strtp(1,1,1), nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'peslofrc - Reduce str' )
        call MPBCAST8( master, nreduc, strtp , icomm )
C        call MPI_Bcast( strtp, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'peslofrc - Bcast str' )
        do  jd=1,3
          do  id=1,3
            stresl(id,jd) = strtp(id,jd,1)
            strxc0(id,jd) = strtp(id,jd,2)
            strxcc(id,jd) = strtp(id,jd,3)
            strdxc0(id,jd) = strtp(id,jd,4)
            strdxcc(id,jd) = strtp(id,jd,5)
          enddo
        enddo
c
c       End tp
      endif
ctp:end
c
c Put gga-specific stress terms into total stress terms:
      if( do_gga .and. ndim.gt.0 )then
        do  jd=1,ndim
          do  id=1,ndim
            strxc0(id,jd) = strxc0(id,jd) + strdxc0(id,jd)
            strxcc(id,jd) = strxcc(id,jd) + strdxcc(id,jd)
          enddo
        enddo
      endif
c
      if( do_gga .and. ndim.gt.0 )then
        if( lstout.gt.2 )then
          write(IWR,*) 'eslofrc-258a: add new gga stress terms'
          write(IWR,*) 'gga stress grid reference'
 9040     format(5x,3f12.6)
          do  jd=1,ndim
            write(IWR,9040) (strdxc0(id,jd), id=1,ndim)
          enddo
          write(IWR,*) 'gga stress core grid reference'
          do  jd=1,ndim
            write(IWR,9040) (strdxcc(id,jd), id=1,ndim)
          enddo
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
