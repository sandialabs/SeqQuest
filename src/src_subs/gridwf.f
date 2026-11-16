c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDWF
c
c
      subroutine GRIDWF( ncplx,idwf,
     $ nk,  natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $ coskr,sinkr,
     $ hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wvfcns,wf,
     $ wkbox )
c---------------------------------------------------------------
c Purpose: generate basis orbitals on regular space grid
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c   4Dec08-PAS/2.63: k-parallel scf matrices
c   4Jan08-PAS/2.61: reduced (dist) memory for grid wvfcns
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: cosmetic modification
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wf
c
c Input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  hh(3,3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
c Output arrays:
      DIMENSION  wvfcns(2,*),wf(*)
c Scratch arrays:
      DIMENSION  wkbox(*)
c
c Local declarations:
      DIMENSION  r0box(3),boxctr(3)
c
c >>>> EXECUTABLE CODE:
c
      if( idwf.eq.0 ) RETURN
c
c  Get Local MP info
      call MPNODES_K( nprocs )
      call MPNODE0_K( master )
      call MPNODE_K( iproc )
      call MPCOMM_K( icomm )
c
c  Set up k-parallel
      call KPMINE( nk,nkloc,nk0 )
c
      ir00 = 0
c
c  Cycle over boxes
c
      do 3000 ibox=1,nbox
c
c       Each processor will do every nprocs'th box in this loop.
        if( MOD( ibox, nprocs ) .ne. iproc ) goto 3000
c
c       Get number of orbitals in this box; skip if none
        nobox = noboxs(ibox)
        if( nobox.le.0 ) goto 3000
c
c       Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
c       Generate mesh point coordinates
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
        call MKRBOX( hh,r0box, idra,idrb,idrc, wkbox )
c
c       Compute orbitals on mesh in this box
c
        iwk1 = 1
        iwk2 = iwk1 + 3*nrbox
        iwk3 = iwk2 + nrbox
        iwk4 = iwk3 + nrbox
        iwk5 = iwk4 + nrbox
        iwk6 = iwk5 + nrbox
        iwk7 = iwk6 + nrbox
        iwk8 = iwk7 + nrbox
c
        call BOXWFN( ibox, ncplx,
     $   nk,nk0,nkloc,  natm,ntyp, nshld,nald, nlat,
     $   itypa,numshl,lshel,nala,ala,cala,alamin,
     $   ratm,rlat,coskr,sinkr,
     $   nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $   wvfcns(1,ir00+1), wf(ir00+1),
     $   wkbox      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $   wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
c       Close box loop:
 2900   continue
        ir00 = ir00 + nkloc*nobox*nrbox
 3000 continue
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
c End:
