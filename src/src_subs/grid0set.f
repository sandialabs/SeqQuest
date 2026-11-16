c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRID0SET
c
c
      subroutine GRID0SET( IWR, lvlout, igrd0fl,idrhfl0,igrhofl,
     $ do_gga,do_spin, anycore,  xc0val,rhoatom,
     $ xcalcut,xcfac,xcrhocut,  maxwk,maxwkr,mxrb,
     $ natm,ntyp,nshld,nald, nlat, nlat1c, nptr, n1r,n2r,n3r,
     $ itypa, numshl,lshel,nala,ala,cala,
     $        nocc,locc,nalsh,alsh,calsh,  znuc, occsh,
     $ nrd,nrcor,radmsh,corden,
     $ ratm,rlat, nearest, hh, weight,
     $ ndbox,nbox,nrec,inrec,inboxs,boxrad,
     $ wksma,wspl,wkrbx,
     $ wk )
c---------------------------------------------------------------
c Purpose: Do grid part of iteration-independent calculations:
c          reference density, reference xc-potentials, and, if
c          we are doing gga, first and second derivatives of
c          reference density.
c
c Written: Peter A. Schultz, 25-May-2001, for v2.47
c
c Revision history:
c  29Apr14-PAS/2.66: master-only grd0 file
c  19Dec07-PAS/2.61: master-only gga grid files
c  15Oct04-PAS/2.59: pass additional space into grid0xc
c  19Jul01-PAS/2.49: igr0dfl change for spin calculation
c  27May01-PAS/2.48: replace STOPs
c  25May01-PAS/2.47: pushed down from main program
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output arrays:
      DIMENSION  xc0val(*)
      DIMENSION  rhoatom(natm,*)
c input declarations:
      LOGICAL    do_gga, do_spin, anycore
      DIMENSION  itypa(natm)
      DIMENSION  numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $  nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  znuc(ntyp), occsh(nshld,ntyp)
      DIMENSION  nrcor(ntyp),radmsh(nrd,ntyp),corden(nrd,ntyp)
      DIMENSION  ratm(3,natm),rlat(3,*), hh(3,3)
      DIMENSION  nearest(3,natm)
      DIMENSION  ndbox(2,2,3), inboxs(*)
c scratch arrays:
      DIMENSION  wksma(natm,*), wspl(nrd,2), wkrbx(maxwkr)
      DIMENSION  wk(maxwk)
c
c local declarations:
      DATA  zero,half,one / 0.d0,0.5d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Carry over the standard output format:
 9020 format(1x,a,3f20.10)
c
c Identify the master processors:
      call MPNODE( iproc )
      call MPNODE0( master )
c
c Space begins at 1st array space in subroutine:
      i01 = 1
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Do grid part of iteration-independent calculations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Partition space in big work array, and check for sufficiency:
c
c Note that the memory need in "potnr0" cannot be well described
c outside of it, and is not protected in the following statements.
c Normally, there will be plenty of memory, but there is a check
c inside "potnr0" to make sure it stays within bounds, just in case.
      j01 = i01
      j02 = j01 + nptr
      j03 = j02 + nptr
      j04 = j03 + nptr
      j05 = j04 + nptr
      j06 = j05 + nptr
      if( do_gga )then
        j06d = j06
        needwk = j06 - 1 + 9*nptr
      else
c       Set j06d so call to grid0xc (unused) derivs is surely in-bounds:
        j06d = j05
        needwk = j06 - 1
      endif
      if( needwk .gt. maxwk )then
        call STOPXERR( 'gr0set-1/short wk() memory in setup grid' )
      endif
c
c  Compute nearby potential for each reference atom on grid:
c
cspin: need to fix nearby stuff when spin is on ...
      call POTNR0( IWR, do_gga, wksma,
     $ xcalcut,xcfac,xcrhocut,
     $ natm,ntyp,nshld,nald, nlat1c, nptr,
     $ itypa, numshl,lshel,nala,ala,cala,
     $        nocc,locc,nalsh,alsh,calsh,  znuc,occsh,
     $ nrd,nrcor,radmsh,corden, wspl(1,1),wspl(1,2),
     $ ratm,rlat, nearest, hh,
     $ wk(j01), maxwk )
c -->  wkr-8s
c
      if( lvlout.gt.2 ) call TIMER('After potnr0     ')
c
      call MKZERO( nptr, wk(j01) )
      call MKZERO( nptr, wk(j02) )
      call MKZERO( nptr, wk(j03) )
      call MKZERO( nptr, wk(j04) )
c
c  Compute sums of grid densities/xc-potentials for reference atoms:
c
c  Compute space in small scratch workspace:
      mxwkb = maxwkr - 3*mxrb
      if( mxwkb.lt.0 ) call STOPXERR( 'gr0set-2/wkrbx memory' )
c
      call GRID0XC( IWR, do_gga, wksma,wksma(1,5),
     $ xcalcut,xcfac,xcrhocut,
     $ rhotot,cortot,exc0grd,exc0gsl,exc0pt,exc0wt,exccgrd, rhoatom,
     $ natm,ntyp,nshld,nald, nlat1c, nptr,
     $ itypa, numshl,lshel,nala,ala,cala,
     $        nocc,locc,nalsh,alsh,calsh,  znuc,occsh,
     $ nrd,nrcor,radmsh,corden, wspl(1,1),wspl(1,2),
     $ ratm,rlat, hh, weight,
     $ ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $ mxwkb, wkrbx,   wkrbx(1+3*mxrb),
c -->         rv(3B)-s wkr(11B)-s
     $ wk(j01), wk(j02), wk(j03), wk(j04), wk(j06d) )
c -->  rho0-o   vxc0-o   exc0-o   rhoc-o   drho-9o
c As configured, using 5*nptr in lda, and (5+9)=14*nptr in gga
c
c Do some processing (remember all grid fields are in boxes):
c
      write(IWR,*)
      corfast = zero
      if( anycore )then
        cortot = weight*DDOT( nptr, one,0, wk(j04),1 )
        corfast = weight*DDOT( nptr, wk(j03),1, wk(j04),1 )
        if( lvlout.gt.2 )then
          write(IWR,9020)'sum of grid-integrated core densities=',cortot
          write(IWR,9020)'grid-integrated core:Exc0 =  ',corfast
        endif
      endif
      if( lvlout.gt.3 )then
        write(IWR,9020)'grid-integrated isolated atom exc =   ',exc0grd
        write(IWR,9020)'atom exc contribution - nearby correct =',exc0pt
        write(IWR,9020)'atom exc contribution - weighted nearby=',exc0wt
      endif
      xc0val(1) = exc0grd
      xc0val(3) = exc0pt
      xc0val(4) = exc0wt
      xc0val(5) = exccgrd
      xc0val(7) = corfast
      xc0val(9) = exc0gsl
      exc0rad = xc0val(2)
      exccrad = xc0val(6)
c
      cor0fix = exccrad - exccgrd
      exc0fix = exc0rad - exc0grd + cor0fix
      if( lvlout.gt.3 )then
       write(IWR,*)
       write(IWR,9020)'sum radial-integrated isolated atom exc=',exc0rad
       write(IWR,9020)'raw exc correction from isolated atoms =',exc0fix
       if( anycore )
     $ write(IWR,9020)'exc core correction from isolated atoms=',cor0fix
      endif
c
      if( nlat1c.lt.nlat .and.lvlout.gt.3 )
     $ write(IWR,*) 'for 1c terms, nlat reduced from',nlat,' to',nlat1c
c     Reset nlat if system totally contained in primary unit cell
      if( nlat1c.eq.1 .and. nlat.ne.1 ) call STOPXERR('bad nlat gr0set')
c
c Save reference atom densities and xc-potentials to grid ...
c
      call DCOPY( nptr, wk(j01),1, wk(j05),1 )
      call FROMBOX( n1r,n2r,n3r, wk(j01), ndbox, wk(j05) )
      call DCOPY( nptr, wk(j02),1, wk(j05),1 )
      call FROMBOX( n1r,n2r,n3r, wk(j02), ndbox, wk(j05) )
      call DCOPY( nptr, wk(j03),1, wk(j05),1 )
      call FROMBOX( n1r,n2r,n3r, wk(j03), ndbox, wk(j05) )
      call DCOPY( nptr, wk(j04),1, wk(j05),1 )
      call FROMBOX( n1r,n2r,n3r, wk(j04), ndbox, wk(j05) )
c
      if( iproc .eq. master )then
        REWIND( unit=igrd0fl )
c   ... sum of isolated atom densities to igrd0fl:
        call WRITBIG( igrd0fl, nptr,  wk(j01) )
c   ... core charge on grid to igrd0fl:
        if( anycore ) call WRITBIG( igrd0fl, nptr, wk(j04) )
c   ... sum of isolated cut-off atom xc-potentials to igrd0fl:
        call WRITBIG( igrd0fl, nptr,  wk(j02) )
c   ... sum of isolated cut-off atom exc-pots to igrd0fl:
        call WRITBIG( igrd0fl, nptr,  wk(j03) )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After grid0xc    ')
c
      if( do_gga .and. ( iproc.eq.master ) )then
c
c       We need to construct initial density gradients
c       from the derivatives of the reference density
c
c       Unbox all 9 derivative fields and write to disk:
c
        REWIND( unit=idrhfl0 )
        jxx = j06
        do  idrh=1,9
          call DCOPY( nptr, wk(jxx),1, wk(j05),1 )
          call FROMBOX( n1r,n2r,n3r, wk(jxx), ndbox, wk(j05) )
          call WRITBIG( idrhfl0, nptr, wk(jxx) )
          jxx = jxx + nptr
        enddo
c
c       Take reference derivatives, and produce gga gradients:
c
        REWIND( unit=idrhfl0 )
        call MKGRHO0( idrhfl0,   nptr,
     $   wk(j01),wk(j02), wk(j03),  wk(j05),wk(j04) )
c -->    grho-1o g2rho-1o gpggp-1o  drho-3s ddrho-1s
c
        REWIND( igrhofl )
c       Save grho-0, g2rho-0, and gpggp-0 to igrhofl:
        call WRITBIG( igrhofl, nptr, wk(j01) )
        call WRITBIG( igrhofl, nptr, wk(j02) )
        call WRITBIG( igrhofl, nptr, wk(j03) )
        if( do_spin )then
c         Need to save up/dn gradients:
          call DSCAL( nptr, half, wk(j01),1 )
          call DSCAL( nptr, half, wk(j02),1 )
          call DSCAL( nptr, half*half, wk(j03),1 )
c         Save grho-0up/dn, g2rho-0up, and gpggp-0up to igrhofl:
          call WRITBIG( igrhofl, nptr, wk(j01) )
          call WRITBIG( igrhofl, nptr, wk(j01) )
          call WRITBIG( igrhofl, nptr, wk(j02) )
          call WRITBIG( igrhofl, nptr, wk(j03) )
c         Save grho-0dn, g2rho-0dn, and gpggp-0dn to igrhofl:
          call WRITBIG( igrhofl, nptr, wk(j01) )
          call WRITBIG( igrhofl, nptr, wk(j02) )
          call WRITBIG( igrhofl, nptr, wk(j03) )
        endif
c
        if( lvlout.gt.2 ) call TIMER('After mkgrho0    ')
c
      endif
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
