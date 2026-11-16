c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCPOT
c
c
      subroutine VLOCPOT( IWR, do_gga, xcalcut,xcfac,xcrhocut,
     $ exc0rad,exc0slo,exccrad, strxc0,strxcc,
     $ rhoatom, excatom,excaslo,exccore, strvxc0,strxcor,
     $ natm,ntyp,nshld,nald,
     $ itypa, numshl,lshel,nala,ala,cala, znuc,occ,lmxnlp1,
     $ nrad,nrd, radmsh,radwt, nrcor, corden,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad,
     $ maxwkr, wkr, y2,u )
c---------------------------------------------------------------
c Purpose: construct local electrostatic and xc-potentials for
c          spherically averaged isolated atoms, used in vlocmat
c
c Written: Peter A. Schultz, 13-November-2000, for v2.46 (gga)
c          Replacing former routine "atpot"
c
c Revision history:
c  24Sep08-PAS/2.62: r2l initialization (cosmetic)
c  27Feb07-PAS/2.60: adjust output flags
c   6Nov04-PAS/2.58a: add dGGA stress analytic-grid correction
c  11Jun03-PAS/2.55: change DFTGGA arguments
c   8Apr02-PAS/2.53: gga stress
c   7Mar02-PAS.2,52: explicit constants to d0
c  21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c  Notes:
c    Potentials evaluated on same radial mesh as non-local pots
c
c    The electrostatic potential may be written as:
c           2 * occat / r + delta ves(r),
c    where occat is the number of electrons on the atom
c    Only the remainder term is computed.
c
c    This routine coded to handle general shells (l=0,1,2,3,...)
c
c    Arrays vat/vxc/excrad() used for scratch before output
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input arrays:
      LOGICAL    do_gga
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  occ(nshld,ntyp)
      DIMENSION  znuc(ntyp),  lmxnlp1(ntyp)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  vpsrad(nrd,4,ntyp)
      DIMENSION  nrcor(ntyp), corden(nrd,ntyp)
c Output arrays:
      DIMENSION  rhoatom(natm,3), vatrad(nrd,ntyp)
      DIMENSION  vesrad(nrd,ntyp),vxcrad(nrd,ntyp),excrad(nrd,ntyp)
c Scratch arrays:
      DIMENSION  excatom(ntyp),excaslo(ntyp),exccore(ntyp)
      DIMENSION  strvxc0(ntyp), strxcor(ntyp)
      DIMENSION  wkr(maxwkr)
      DIMENSION  y2(*),u(*)
c
c Local declarations:
      LOGICAL  cntrch
      DATA  alcutfac / 1.33333333d0 /
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      fourpi = pi + pi + pi + pi
      rtpid2 = half*rtpi
c determine largest exp argument before underflow:
      argmx = half*cutexp
c
      strxc0 = zero
      strxcc = zero
c
c  Loop over atom types:
c
      do 500 ktyp=1,ntyp
        excatom(ktyp) = zero
        excaslo(ktyp) = zero
        exccore(ktyp) = zero
        strvxc0(ktyp) = zero
        strxcor(ktyp) = zero
        if( znuc(ktyp).eq.zero ) goto 500
c
        almin = 9.d9
        do  ishl=1,numshl(ktyp)
          if( occ(ishl,ktyp).ne.zero )then
            do  ial=1,nala(ishl,ktyp)
              ali = ala(ial,ishl,ktyp)
              if( ali.lt.almin ) almin = ali
            enddo
          endif
        enddo
c
c       Set up xc cutoff function
        alminxc = two*almin
        if( xcalcut.lt.zero ) alminxc = -xcalcut*almin
        alcutxc = alcutfac*alminxc
        if( xcalcut.gt.zero ) alcutxc = xcalcut
c
c       Blank radial potential arrays
        nradk = nrad(ktyp)
        call MKZERO( nradk, vesrad(1,ktyp) )
c
c       Set up scratch space in wkr()
        ir0 = 0
        irsq  = ir0
        irngr = ir0 +  1*nradk
        iatrh = ir0 +  2*nradk
        icnrh = ir0 +  3*nradk
        icore = ir0 +  4*nradk
        idcor = ir0 +  5*nradk
        ia0e  = icore
        ia1e  = idcor
        ia2e  = ir0 +  6*nradk
        idrho = ir0 +  7*nradk
        iddrh = ir0 +  8*nradk
        inext = iddrh + nradk
        if( inext.ge.maxwkr ) call STOPXERR( 'vlp-mem' )
c
c       Clear densities and derivatives, initialize radial terms
        call MKZERO( nradk, wkr(iatrh+1) )
        call MKZERO( nradk, wkr(icnrh+1) )
        call MKZERO( nradk, wkr(idrho+1) )
        call MKZERO( nradk, wkr(iddrh+1) )
        call MKZERO( nradk, wkr(idcor+1) )
        do  nr=1,nradk
          wkr(irngr+nr) = radmsh(nr,ktyp)
          wkr(irsq+nr)  = radmsh(nr,ktyp)**2
        enddo
c
        iorbat = 0
        occat = zero
        do 70 ishl=1,numshl(ktyp)
          lshl = lshel(ishl,ktyp)
c
          occshl = occ(ishl,ktyp)
          occat = occat + occshl
          if( occshl.eq.zero ) goto 70
c
          cntrch = .false.
          if( occshl.lt.zero )  cntrch = .true.
c
c         Shell is occupied - begin alpha loops
c
c         These arrays used here for scratch to store
c          sums of ca*cb*(a+b)**n*exp(-r2*(a+b))
          call MKZERO( nradk, wkr(ia0e+1) )
          call MKZERO( nradk, wkr(ia1e+1) )
          call MKZERO( nradk, wkr(ia2e+1) )
c
          do 55 ial=1,nala(ishl,ktyp)
            ali = ala(ial,ishl,ktyp)
            cali = occshl*cala(ial,ishl,ktyp)
            cali0 = cali
c
c           Do off-diagonals once, using cali and cali0 to double.
            do 50 jal=ial,1,-1
              als = ali + ala(jal,ishl,ktyp)
              twoals = als + als
              cicj = cali*cala(jal,ishl,ktyp)
              cali = cali0 + cali0
c             Could put twoals rather than als for slight speedup
              cicja0 = cicj/fourpi
              cicja1 = cicja0*als
              cicja2 = cicja0*als*als
              rtal = SQRT( als )
              alslp1 = als**(lshl+1)
c
              do 40 nr=1,nradk
                r = radmsh(nr,ktyp)
                r2 = r*r
                rtalr = r*rtal
                ar2  = als*r2
c
                if( ar2.lt.argmx )then
                  expar2 =  EXP( -ar2 )
                  sumlt = zero
                  ar2tol = one
                  sumgt = zero
                  r2ll = one
                  erffac = rtpid2/rtal
c                 Aren't recursion relations beautiful?
                  do  l=0,lshl
                    r2l = r2ll
                    sumlt = ar2tol + l*sumlt
                    ar2tol = ar2*ar2tol
                    sumgt = ((2*l+1)*sumgt-r2ll) / twoals
                    r2ll = r2*r2ll
                    erffac = (2*l+1)*erffac / twoals
                  enddo
                  atsr0 = zero
c                 PRECISION question  ERF[single]<->DERF[double]
                  if( r.gt.1.d-9 ) atsr0 = erffac*DERF( rtalr )/r
                  atsr = expar2*(half*sumlt/alslp1+sumgt) + atsr0
                else
                  r2l = zero
                  expar2 = zero
                  atsr0 = rtpid2/rtalr
                  do  l=0,lshl
                    atsr0 = (2*l+1)*atsr0 / twoals
                  enddo
                  atsr = atsr0
                endif
c
                if( .not. cntrch )then
                  wkr(ia0e+nr) = wkr(ia0e+nr) + cicja0*expar2
                  wkr(ia1e+nr) = wkr(ia1e+nr) + cicja1*expar2
                  wkr(ia2e+nr) = wkr(ia2e+nr) + cicja2*expar2
                else
                  wkr(icnrh+nr) = wkr(icnrh+nr) + cicja0*(r2l)*expar2
                endif
                vesrad(nr,ktyp) = vesrad(nr,ktyp) + cicj*atsr
c
   40         continue
c
c             Close alpha loops:
c
   50       continue
   55     continue
c
          if( .not. cntrch )then
c
c           Assemble density and derivatives for XC potential:
c
            call MKSHL2D( lshl, nradk, radmsh(1,ktyp), wkr(irsq+1),
     $       wkr(ia0e+1),wkr(ia1e+1),wkr(ia2e+1),
     $       wkr(iatrh+1),wkr(idrho+1),wkr(iddrh+1) )
c
          endif
c
c         Close shell loop:
c
   70   continue
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Subtract 1/r term from es potential
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        esrho = zero
        eszat = zero
        do  nr=1,nradk
          rhoes = radwt(nr,ktyp)*(wkr(iatrh+nr)+wkr(icnrh+nr))
          esrho = esrho + rhoes*vesrad(nr,ktyp)
          eszat = eszat - rhoes*occat/radmsh(nr,ktyp)
        enddo
        esrho = fourpi*esrho
        eszat = fourpi*eszat
        esat  = -esrho + eszat + eszat
        if( lstout.gt.4 ) write(IWR,*) 'esrho,z,tot=',esrho,eszat,esat
c
        do  nr=1,nradk
c         Add factor of two to convert to rydbergs
          vesrad(nr,ktyp) = two*( vesrad(nr,ktyp)-occat/radmsh(nr,ktyp))
        enddo
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c            Check if we have core density to consider
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        nradc = nrcor(ktyp)
        nrsc = nradc + 1
        call MKZERO( nradk, wkr(idcor+1) )
        if( nrsc.ne.1 )then
c         Interpolate core density
          ddy0 = 1.d31
          ddy1 = 1.d31
          call SPLINE(nrsc,radmsh(1,ktyp),corden(1,ktyp),ddy0,ddy1,y2,u)
c
c         Do core 2nd deriv before 1st deriv, and use (idcor)
c         space for scratch  (do not need core 2nd deriv
c         separate from the total rho 2nd deriv).
c         Interpolate to get 2nd deriv. of core density ...
          call SPLIND2( nrsc, radmsh(1,ktyp), corden(1,ktyp),y2,
     $                  nrsc, radmsh(1,ktyp),wkr(idcor+1) )
c          ... and put into total 2nd deriv of density for GGA
          do  nr=1,nrsc
            wkr(iddrh+nr) = wkr(iddrh+nr) + wkr(idcor+nr)
          enddo
c
c         Interpolate to get deriv. of core density ...
          call SPLIND( nrsc, radmsh(1,ktyp),corden(1,ktyp),y2,
     $                 nrsc, radmsh(1,ktyp),wkr(idcor+1) )
        endif
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c            Have atom density; compute xc quantities
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        call MKZERO( nradk, vatrad(1,ktyp) )
        call MKZERO( nradk, vxcrad(1,ktyp) )
        call MKZERO( nradk, excrad(1,ktyp) )
        rhoatm = zero
        rhoatx = zero
        coratm = zero
        vxcatm = zero
        vxcfst = zero
        vxcasl = zero
        excatm = zero
        excfst = zero
        excasl = zero
        exccor = zero
        corfst = zero
        strax0 = zero
        straxc = zero
        strd0 = zero
        strdc = zero
c
        do 200 nr=1,nradk
c
          rhoxc = wkr(iatrh+nr)
          rhoes = wkr(iatrh+nr) + wkr(icnrh+nr)
c
          rhopt = rhoxc
          if( nradc.ne.0 ) rhopt = rhopt + corden(nr,ktyp)
c
          if( rhopt.lt.xcrhocut ) goto 200
c
          rmag = radmsh(nr,ktyp)
          if( do_gga )then
c           Assemble derivatives for GGA pot/assume no pt at r=0
c           Remember to add core deriv to total density deriv:
            drh = wkr(idrho+nr) + wkr(idcor+nr)
            d1rhdr = drh/rmag
            d2rh = wkr(iddrh+nr) - d1rhdr
            grho = ABS( drh )
            g2rho = wkr(iddrh+nr) + two*d1rhdr
            gpggp = grho*wkr(iddrh+nr)
c
            call DFTGGA( rhopt, grho,g2rho,gpggp,
     $                   vxcpt,excpt, dxcpt,dxc2pt )
c
          else
c
            call DFTLDA( rhopt, vxcpt,excpt )
c
          endif
c
c         Construct "fast" xc-potentials, and update:
          cfac = EXP( -alcutxc*rmag**2 )
          fstfcn = cfac / ( (one-xcfac)*cfac + xcfac )
c
          vfast = vxcpt*fstfcn
          vslow = vxcpt - vfast
c
          efast = excpt*fstfcn
          eslow = excpt - efast
c
          vxcrad(nr,ktyp) = vfast
          excrad(nr,ktyp) = efast
c
c         Integrate atom electron density and xc quantities
c
c         In the following, "radwt" assumed to contain the factor r**2
c         needed in the radial integral but not the 4pi for the radial
c         integral.  This specification of radwt is more convenient
c         in the more heavily called routines "vlrad" and "vlradX"
          rwt = fourpi*radwt(nr,ktyp)
          rhoatm = rhoatm + rwt*      rhoes
          rhoatx = rhoatx + rwt*      rhoxc
          vxcatm = vxcatm + rwt*vxcpt*rhoxc
          vxcfst = vxcfst + rwt*vfast*rhoxc
          vxcasl = vxcasl + rwt*vslow*rhoxc
          excatm = excatm + rwt*excpt*rhoxc
          excfst = excfst + rwt*efast*rhoxc
          excasl = excasl + rwt*eslow*rhoxc
c         Phantom xc force "analytic" term
          strax0 = strax0 + rwt*vfast*wkr(idrho+nr)*rmag
          if( do_gga )then
c           Phantom GGA stress "analytic" term
            dxc0fast = fstfcn*dxcpt*wkr(idrho+nr)**2
     $               + fstfcn*dxc2pt*wkr(iddrh+nr)
            strd0 = strd0 - rwt*dxc0fast
          endif
          if( nradc.gt.0 .and. nr.le.nradc )then
            coratm = coratm + rwt*      corden(nr,ktyp)
            exccor = exccor + rwt*excpt*corden(nr,ktyp)
            corfst = corfst + rwt*efast*corden(nr,ktyp)
c           Phantom xc force "analytic" term (core)
            straxc = straxc + rwt*vfast*wkr(idcor+nr)*rmag
            if( do_gga )then
c             Phantom GGA stress "analytic" term (core)
c             Note: d2cor term in iddrho above
              dxccfast = fstfcn*dxcpt*wkr(idcor+nr)**2
              strdc = strdc - rwt*dxccfast
            endif
          endif
c
  200   continue
c
        excatom(ktyp) = excatm
        excaslo(ktyp) = excasl
        exccore(ktyp) = exccor
        strvxc0(ktyp) = strax0 + strd0
        strxcor(ktyp) = straxc + strdc
        if( lstout.gt.0 )then
          write(IWR,*) ' '
          write(IWR,9200) ktyp,' radial valence charge =',rhoatm
 9200     format(1x,'Atom type',i4,1x,a,f18.12)
          if( ABS(rhoatx-rhoatm) .gt. 1.d-8 )
     $    write(IWR,9200) ktyp,'     xc valence charge =',rhoatx
          if( lstout.gt.3 )then
            write(IWR,9200) ktyp,'    vxc rad energy =',vxcatm
            write(IWR,9200) ktyp,'    vxc fast enrgy =',vxcfst
            write(IWR,9200) ktyp,'    vxc slo energy =',vxcasl
            write(IWR,9200) ktyp,'    exc rad energy =',excatm
            write(IWR,9200) ktyp,'    exc fast enrgy =',excfst
            write(IWR,9200) ktyp,'    exc slo energy =',excasl
            write(IWR,9200) ktyp,'    coulomb energy =',esrho
            write(IWR,9200) ktyp,'    valence stress =',strax0
            write(IWR,9200) ktyp,'    gga stress     =',strd0
          endif
          if( nradc.gt.0 )then
            write(IWR,9200) ktyp,'   core charge =',coratm
            if( lstout.gt.3 )then
              write(IWR,9200) ktyp,'    core rad energy =',exccor
              write(IWR,9200) ktyp,'    core fast enrgy =',corfst
              write(IWR,9200) ktyp,'    core stress     =',straxc
              write(IWR,9200) ktyp,'    core gga stress =',strdc
            endif
          endif
          rhoatom(ktyp,3) = coratm
        endif
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       use the "non-local" potential of the highest l-value
c                  as the local potential:
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        locpot = lmxnlp1(ktyp) + 1
        if( locpot.le.0 )then
          do  nr=1,nradk
            vatrad(nr,ktyp) = vesrad(nr,ktyp) + vxcrad(nr,ktyp)
          enddo
        else
          do  nr=1,nradk
            vatrad(nr,ktyp) = vesrad(nr,ktyp) + vxcrad(nr,ktyp)
     $                      + vpsrad(nr,locpot,ktyp)
          enddo
        endif
c
c      Close atom loop:
c
  500 continue
c
c Assemble isolated atom correction to xc-energy
c
      exc0rad = zero
      exc0slo = zero
      exccrad = zero
      strxc0  = zero
      strxcc  = zero
      do  iatm=1,natm
        ityp = itypa(iatm)
        exc0rad = exc0rad + excatom(ityp)
        exc0slo = exc0slo + excaslo(ityp)
        exccrad = exccrad + exccore(ityp)
        strxc0  = strxc0  + strvxc0(ityp)
        strxcc  = strxcc  + strxcor(ityp)
      enddo
c
      if( lstout.gt.1 )then
c       Summarize atomic reference xc stress terms:
        write(IWR,9220) 'Total atom local valence xc stress=',strxc0
        if( strxcc .ne. zero )
     $  write(IWR,9220) 'Total atom local core xc stress   =',strxcc
 9220   format(1x,a,f20.8)
      endif
c
c    That's all Folks!
c
      RETURN
      END
