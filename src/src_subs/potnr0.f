c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POTNR0
c
c
      subroutine POTNR0( IWR, do_gga, xcnrat,
     $ xcalcut,xcfac,xcrhocut,
     $ natm,ntyp,nshld,nald, nlat  , nptr,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,  nrd,nrcor,radmsh,corden, y2, u,  ratm,rlat, nearest,
     $ hh, wkr, mxwkrd )
c---------------------------------------------------------------
c Purpose: Compute "nearby" potentials for grid reference atoms
c
c Written: Peter A. Schultz,  9-December-2000, for v2.46 (gga)
c
c Revision history:
c  11Jun03-PAS/2.55: change DFTGGA arguments
c   8Apr02-PAS/2.53: gga stress
c   7Mar02-PAS/2.52: explicit constants to d0
c  15Aug01-PAS/2.48b: bugfix!!!! (gga/array not initialized)
c  21Jun01-PAS/2.48: replace STOPs
c   9Dec00-PAS/2.46: split off from grid0xc
c---------------------------------------------------------------
c
c Notes:
c  This routine should exactly reproduce the computation of the
c  potentials done in "grid0xc".  The distinction here is that all
c  we care about is a weighted sampling of the nearby potential,
c  for use in refining a correction term to the xc-energy to fix
c  inaccuracies in the grid integrals.  I separate this code (just
c  as I separated "eslofrc", the near-copy of "grid0xc" to used to
c  compute the phantom force) for clarity in coding, and to ease
c  complexity in memory management in "grid0xc".  I judge it
c  simpler to duplicate the computational task and code, than to
c  code in "grid0xc" the logic of how to manage memory and squeeze
c  the weighting operation into boxes.  This way, this weighting
c  operation can more cleanly be matched to that used in "potnear".
c
c  *** NB ***  Code is kept as close to that used in "grid0xc"
c  as possible, even when such code locally here seems redundant
c  or non-conforming, to ease later changes where common code
c  betweeen "grid0xc", "eslofrc", and "potnr0" is altered.
c  Changes to these routines should only be done in tandem.
c
c  The use of cutoffs in this routine should be compatible
c  with that used in calculating the grid wave functions.
c
c  The code for the nearby average must match that used in "potnear".
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c output array:
      DIMENSION  xcnrat(4,natm)
c input arrays:
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
      DIMENSION  nearest(3,*)
c scratch arrays:
      DIMENSION  y2(*),u(*)
      DIMENSION  wkr(*)
c
c local declarations:
      DIMENSION  rnr(3), r0(3),rr1(3),rr2(3),rr3(3), gg(3,3)
      DIMENSION  nrng(3),nnr(3)
c
c  nearby sampling function parameters:
      DATA  alfsmpl,rsampl / 1.2d0,1.75d0 /
      DATA  alcutfac / 1.33333333d0 /
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  three / 3.d0 /
c flag for additional output 0-none/1-some/2-full
      DATA  lstout  / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( lstout.gt.0 ) write(IWR,*) 'POTNR0: grid reference nearby'
      fourpi = pi + pi + pi + pi
c
c  Set up ranges for nearby stuff:
c
      call CROSS( hh(1,1),hh(1,2), gg(1,3) )
      call CROSS( hh(1,2),hh(1,3), gg(1,1) )
      call CROSS( hh(1,3),hh(1,1), gg(1,2) )
      g1mag = SQRT( gg(1,1)**2 + gg(2,1)**2 + gg(3,1)**2 )
      g2mag = SQRT( gg(1,2)**2 + gg(2,2)**2 + gg(3,2)**2 )
      g3mag = SQRT( gg(1,3)**2 + gg(2,3)**2 + gg(3,3)**2 )
      do  ir=1,3
        gg(ir,1) = gg(ir,1)/g1mag
        gg(ir,2) = gg(ir,2)/g2mag
        gg(ir,3) = gg(ir,3)/g3mag
      enddo
      r1eff = ABS( hh(1,1)*gg(1,1) + hh(2,1)*gg(2,1) + hh(3,1)*gg(3,1) )
      r2eff = ABS( hh(1,2)*gg(1,2) + hh(2,2)*gg(2,2) + hh(3,2)*gg(3,2) )
      r3eff = ABS( hh(1,3)*gg(1,3) + hh(2,3)*gg(2,3) + hh(3,3)*gg(3,3) )
c
c  Clear out output array of nearby potentials:
c
      call MKZERO( 4*natm, xcnrat )
c
c >>>>>>>>>> Cycle over atoms:
c
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
c
        nshli = nocc(ityp)
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
c  Bring out xc slow/fast function parameters
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
c        [NB: code in "grid0nr" and "potnear" must match]
c
c  Compute "nearby" xc energy density if grid contains atom
c
        wts = 1.0
        vxcwtnr = zero
        excwtnr = zero
        vxcptnr = zero
        excptnr = zero
c
        do  jr=1,3
          r0(jr) = (nearest(1,iatm)-1)*hh(jr,1)
     $           + (nearest(2,iatm)-1)*hh(jr,2)
     $           + (nearest(3,iatm)-1)*hh(jr,3)
        enddo
c
c for grins, find "nearest" image closest to atom
        ilatnr = 1
        rsqnr = zero
        do  jr=1,3
          rnr(jr) = r0(jr) - ratm(jr,iatm)
          rsqnr = rsqnr + (rnr(jr)-rlat(jr,1))**2
        enddo
        do  ilat=1,nlat
          rsqnr2 = zero
          do  jr=1,3
            rsqnr2 = rsqnr2 + (rnr(jr)-rlat(jr,ilat))**2
          enddo
          if( rsqnr2.lt.rsqnr )then
            ilatnr = ilat
            rsqnr = rsqnr2
          endif
        enddo
        rnr(1) = r0(1) - ratm(1,iatm) - rlat(1,ilatnr)
        rnr(2) = r0(2) - ratm(2,iatm) - rlat(2,ilatnr)
        rnr(3) = r0(3) - ratm(3,iatm) - rlat(3,ilatnr)
c
c  Determine range in grid point dimensions for nearby sampling:
        rng1 = (rsampl/r1eff + half)
        rng2 = (rsampl/r2eff + half)
        rng3 = (rsampl/r3eff + half)
        nrng(1) = rng1
        nrng(2) = rng2
        nrng(3) = rng3
        nnr(1) = 2*nrng(1) + 1
        nnr(2) = 2*nrng(2) + 1
        nnr(3) = 2*nrng(3) + 1
        nrsmpl = nnr(1)*nnr(2)*nnr(3)
        nrpt = (nrsmpl+1) / 2
c
c  to accentuate similarity with "grid0xc":
        nrbox = nrsmpl
c  Set up temp space
        ir0 = 0
        irsq = ir0
        irngr = ir0 + 1*nrbox
        iatrh = ir0 + 2*nrbox
        icore = ir0 + 3*nrbox
        idcor = ir0 + 4*nrbox
        ia0e  = icore
        iwts  = ir0 + 5*nrbox
        inext = iwts + nrbox
        ia1e  = idcor
        ia2e  = iwts
        idrho = ir0 + 6*nrbox
        iddrh = ir0 + 7*nrbox
        inext = iddrh + nrbox
c I could neglect this memory check as I am in BIG workspace.
        if( inext.ge.mxwkrd ) call STOPXERR( 'nr0-wkrd' )
c
c  Generate r**2 on local grid:
c
        do  jr=1,3
          rr3(jr) = rnr(jr)
          do  iprim=1,3
            rr3(jr) = rr3(jr) - (DBLE(nrng(iprim))+1)*hh(jr,iprim)
          enddo
          rnr(jr) = rr3(jr) + hh(jr,1) + hh(jr,2) + hh(jr,3)
        enddo
        nr = 0
c
        do 1053 i3=1,nnr(3)
          do  jr=1,3
            rr3(jr) = rr3(jr) + hh(jr,3)
            rr2(jr) = rr3(jr)
          enddo
          do 1052 i2=1,nnr(2)
            do  jr=1,3
              rr2(jr) = rr2(jr) + hh(jr,2)
              rr1(jr) = rr2(jr)
            enddo
            do 1051 i1=1,nnr(1)
              rr1(1) = rr1(1) + hh(1,1)
              rr1(2) = rr1(2) + hh(2,1)
              rr1(3) = rr1(3) + hh(3,1)
              nr = nr + 1
c
              wkr(irsq+nr) = rr1(1)**2 + rr1(2)**2 + rr1(3)**2
              wkr(irngr+nr) = SQRT( wkr(irsq+nr) )
c
              wkr(iatrh+nr) = zero
c
 1051       continue
 1052     continue
 1053   continue
c
        if( do_gga )then
          do  nr=1,nrbox
            wkr(idrho+nr) = zero
            wkr(iddrh+nr) = zero
          enddo
        endif
c
        do 1800 ishl=1,nshli
          lshl = locc(ishl,ityp)
          twol = DBLE( 2*lshl )
c
          occshl = occsh(ishl,ityp)
          if( occshl.eq.zero ) goto 1800
c
          call MKZERO( nrbox, wkr(ia0e+1) )
          if( do_gga )then
            do  nr=1,nrbox
C             wkr(ia0e+nr) = zero
              wkr(ia1e+nr) = zero
              wkr(ia2e+nr) = zero
            enddo
          endif
c
          do 1190 ial=1,nalsh(ishl,ityp)
            ali = alsh(ial,ishl,ityp)
            cali = occshl*calsh(ial,ishl,ityp)
            cali0 = cali
            do 1180 jal=ial,1,-1
              alsum = ali + alsh(jal,ishl,ityp)
              als2 = alsum*alsum
              coef = cali*calsh(jal,ishl,ityp)/fourpi
              cali = cali0 + cali0
c
c  Sum over shell: occ*Ca*Cb*(a+b)**n*exp(-(a+b)r**2); n=0,1,2
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
c      Close alpha loops:
c
 1180       continue
 1190     continue
c
          if( do_gga )then
c
c           Assemble density and derivatives for GGA potential:
c
            call MKSHL2D( lshl, nrbox, wkr(irngr+1), wkr(irsq+1),
     $       wkr(ia0e+1),wkr(ia1e+1),wkr(ia2e+1),
     $       wkr(iatrh+1),wkr(idrho+1),wkr(iddrh+1) )
c
          else
c
c           Just need to assemble density for LDA calculation:
c
            call MKSHL( lshl, nrbox, wkr(irsq+1),
     $       wkr(ia0e+1), wkr(iatrh+1) )
c
          endif
c
c      Close shell loop:
c
 1800   continue
c
c  Now check if we have core density to worry about.
c
        call MKZERO( nrbox, wkr(icore+1) )
        if( nrsc.ne.1 )then
c
          if( do_gga )then
c Do core 2nd deriv before first and use (idcor) space for scratch.
c Don't need the core 2nd deriv separate from total rho 2nd deriv.
c Interpolate to get 2nd deriv. of core density ...
            call SPLIND2( nrsc, radmsh(1,ityp), corden(1,ityp),y2,
     $                    nrbox,wkr(irngr+1),wkr(idcor+1) )
c  ... and put into 2nd derivative of density for GGA potential:
            do  nr=1,nrbox
              wkr(iddrh+nr) = wkr(iddrh+nr) + wkr(idcor+nr)
            enddo
c
c Interpolate to get deriv. of core density ...
            call SPLIND( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                   nrbox,wkr(irngr+1),wkr(idcor+1) )
c  ... and put into density derivative for GGA potential calculation
            do  nr=1,nrbox
              wkr(idrho+nr) = wkr(idrho+nr) + wkr(idcor+nr)
            enddo
          endif
c
c Interpolate to get core density on grid from radial core density
          call SPLINT( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                 nrbox,wkr(irngr+1),wkr(icore+1) )
c
        endif
c
c  Local density assembled; now do local exc-dens, and weighted sampling
c
        wtsum = zero
        call MKNRWT( alfsmpl,wtsum,rnr,hh,nnr,nnr(2),nnr(3),wkr(iwts+1))
        wts = zero
        vxcwtnr = zero
        excwtnr = zero
        do 1200 nr=1,nrbox
          atrho = wkr(iatrh+nr)
          corho = wkr(icore+nr)
          rhopt = atrho + corho
c
          if( do_gga )then
c
c           Assemble gradient terms needed for GGA
c
            rmag = wkr(irngr+nr)
            if( rmag.gt.1.d-6 )then
c For r>0, general Cartesian derivs from radial derivs ...
c  ... note that dcor already added into drho above
              drh = wkr(idrho+nr)
c  ... terms for Cartesian reference derivatives:
              d1rhdr = drh/rmag
              d2rh = wkr(iddrh+nr) - d1rhdr
c  ... need gradients for the xc-potential generation:
              grho = ABS( drh )
              g2rho = wkr(iddrh+nr) + two*d1rhdr
              gpggp = grho*wkr(iddrh+nr)
c
            else
c At r=0, dp/du=0, and only d2p/du2 survives limits ...
              d2rh = wkr(iddrh+nr)
              grho = zero
              g2rho = three*wkr(iddrh+nr)
              gpggp = zero
            endif
c
            call DFTGGA( rhopt, grho,g2rho,gpggp,
     $                   vxcpt,excpt, dxcpt,dxc2pt )
c
c           End do_gga = .true.
          else
c           do_gga = .false.
c
            call DFTLDA( rhopt, vxcpt,excpt )
c
          endif
c
          wts = wts + wkr(iwts+nr)
c  Construct "fast" xc-potentials, and update:
          cfac = EXP( -alcutxc*wkr(irsq+nr) )
          fstfcn = cfac / ( (one-xcfac)*cfac + xcfac )
c
          vslow = vxcpt - vxcpt*fstfcn
          eslow = excpt - excpt*fstfcn
c
          vxcwtnr = vxcwtnr + wkr(iwts+nr)*vslow
          excwtnr = excwtnr + wkr(iwts+nr)*eslow
c
          if( nr.eq.nrpt )then
c  Construct "fast" part of xc-potentials, and remove:
            vxcptnr = vxcpt - vxcpt*fstfcn
            excptnr = excpt - excpt*fstfcn
          endif
 1200   continue
        vxcwtnr = vxcwtnr/wtsum
        excwtnr = excwtnr/wtsum
c
c Load the nearby potential values into storage arrays for "grid0xc":
        xcnrat(1,iatm) = vxcwtnr
        xcnrat(2,iatm) = excwtnr
        xcnrat(3,iatm) = vxcptnr
        xcnrat(4,iatm) = excptnr
c
c      Close atom loop:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
