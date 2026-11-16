c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRDOVLP
c
c
      subroutine GRDOVLP( ncplx, idwf, near2c, s1atom,orbint,
     $ norb,nk, natm,ntyp, nshld,nald, noad, nlat,
     $ itypa,numshl,lshel,nala,ala,cala,  alamin,norba,
     $ ratm,rlat,coskr,sinkr,  vk, phkdoti,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wvfcns,wf,  wkbox, rmesh,
     $ eikr, otmp )
c---------------------------------------------------------------
c Purpose: compute overlaps on grid of orbitals (for nearby)
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c  29Jan10-ACP/2.63 & 6Jul12-PAS: CLeanup
c  17Nov08-PAS/2.63: k-parallel scf matrices
c   4Jan08-PAS/2.61: reduced (dist) memory for grid wvfcns
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c   4Jan00-PAS/2.42: cosmetic clean-up
c   3May99-PAS/2.35: toggle orbint use from flag near2c
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: on-the-fly grid orbitals
c---------------------------------------------------------------
c
c Note: orbint&orbofk and sinkr&coskr are in full nk
c       Other arrays (all scratch arrays) are in local nkloc
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      REAL  wvfcns,wf
      REAL  ovlpr,ovlpi
c
c Input/Output arrays:
      DIMENSION  s1atom(noad,noad,nk,natm),orbint(2,norb,nk)
c
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  norba(ntyp),  vk(3,nk)
      DIMENSION  hh(3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
c
c Scratch arrays:
      DIMENSION  wkbox(*)
      DIMENSION  phkdoti(2,*)
c     DIMENSION  phkdoti(2,nkloc)
      DIMENSION  rmesh(3,*),eikr(2,*)
c      needs to match the larger of s1atom() and orbint()
      DIMENSION  otmp(*)
c
c Local declarations:
      DIMENSION  r0box(3),boxctr(3)
c  Criterion of rapidly varying orbital:
      DATA  alfast / 1.2d0 /
      LOGICAL    EVEN
c Statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES( nprocl )
      call MPNODE0( masterl )
      call MPNODE( iprocl )
      call MPCOMM( icomml )
c
      call MPNODES_K( nprocs )
      call MPNODE0_K( master )
      call MPNODE_K( iproc )
      call MPCOMM_K( icomm )
c
c  Set up k-parallel
      call KPMINE( nk,nkloc,nk0 )
c
c  Clear out non-master data, so that reduce counts analytics only once
      if( nprocl .gt. 1 )then
        lenvec2 = 2*norb*nk
        lenvec1 = noad*noad*nk*natm
        if( iprocl .ne. masterl )then
c         Multi-proc, and not the master; clear data
          if( near2c .eq. 2 )then
            call MKZERO( lenvec2, orbint )
          endif
          call MKZERO( lenvec1, s1atom )
        endif
Cc       Put a barrier here to make sure master is ready
C        barr = zero
C        lenbarr = 1
C        call MPBCAST8( masterl, lenbarr, barr, icomml )
      endif
c
c     If nkloc=0, bypass loop, because nothing to do ...
      if( nkloc .eq. 0 ) goto 1001
c
      ir00 = 0
      inbox = 0
      jnbox = 0
c
c Loop over boxes of grid points
c
      do 1000 ibox=1,nbox
c
c       Each processor will do every nprocs'th box in this loop
        if( MOD( ibox, nprocs ) .ne. iproc ) goto 1000
c
c       Get basic box data
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
        nrkbox = nkloc*nrbox
c
c       Get number of orbitals in this box; if none, skip
        nobox = noboxs(ibox)
        if( nobox.le.0 ) goto 1000
c
c       Get grid orbitals:
c
        if( idwf.eq.0 )then
c         Compute orbitals on mesh in this box
c
c         Generate mesh point coordinates
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          call MKRBOX( hh,r0box, idra,idrb,idrc, wkbox )
          ir00 = 0
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
     $     nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $     itypa,numshl,lshel,nala,ala,cala,alamin,
     $     ratm,rlat,coskr,sinkr,
     $     nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $     wvfcns(1,ir00+1), wf(ir00+1),
     $     wkbox      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $     wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
        endif
c
        if( ncplx.eq.2 .and. near2c.eq.2 )then
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          call MKRBOX( hh,r0box, idra,idrb,idrc, rmesh )
c
          kr0 = 0
          do  ikloc=1,nkloc
            ik = nk0 + ikloc
            xk = vk(1,ik)
            yk = vk(2,ik)
            zk = vk(3,ik)
c           calculate exp(ik.r) for evaluation of slowpot counter-terms
            do  nr=1,nrbox
              rdotk = xk*rmesh(1,nr) + yk*rmesh(2,nr) + zk*rmesh(3,nr)
              eikr(1,kr0+nr) = COS( rdotk )
              eikr(2,kr0+nr) = SIN( rdotk )
            enddo
            kr0 = kr0 + nrbox
          enddo
        endif
c
c       Loop over atoms/shells
c
        irec = 0
        iatshl = 0
        irshl = ir00
        inext = 0
        do 200 iatm=1,natm
          ityp = itypa(iatm)
c
          if( ncplx.eq.2 .and. near2c.eq.2 )then
            do  ikloc=1,nkloc
              ik = nk0 + ikloc
              xk = vk(1,ik)
              yk = vk(2,ik)
              zk = vk(3,ik)
              ridotk = xk*ratm(1,iatm)+yk*ratm(2,iatm)+zk*ratm(3,iatm)
              phkdoti(1,ikloc) = COS( ridotk )
              phkdoti(2,ikloc) = SIN( ridotk )
            enddo
          endif
c
          knbox = inbox
          krec = irec
          katshl = iatshl
          krshl = irshl
          iorbat00 = 0
          do 182 ishl=1,numshl(ityp)
            norbshi = 2*lshel(ishl,ityp) + 1
            inbox = inbox/2
            if( inrec*(iatshl/inrec).eq.iatshl )then
              irec = irec + 1
              inbox = inboxs(irec,ibox)
            endif
            iatshl = iatshl + 1
            iorbat0 = iorbat00
            iorbat00 = iorbat00 + norbshi
            i0 = inext
            inext = i0 + norbshi
            if( EVEN(inbox) ) goto 182
c
            nali = nala(ishl,ityp)
            almaxi = ala(nali,ishl,ityp)
            do  ial=1,nali-1
              ali = ala(ial,ishl,ityp)
              if( ali.gt.almaxi ) almaxi = ali
            enddo
c
c           Take care of 2-center terms first:
c
            if( near2c.eq.2 .and. almaxi.ge.alfast )then
c             This is a fast orbital, do this term
              iorb = i0
              irorb = irshl
              do 170 iorbsh=1,norbshi
                iorb = iorb + 1
                ir0 = irorb
                if( ncplx.eq.1 )then
                  osumr = zero
                  do  nr=1,nrbox
                    osumr = osumr + wf(ir0+nr)
                  enddo
                  orbint(1,iorb,1) = orbint(1,iorb,1) - weight*osumr
                  orbint(2,iorb,1) = zero
                elseif( ncplx.eq.2 )then
                  kr0 = 0
                  do  ikloc=1,nkloc
                    ik = nk0 + ikloc
                    osumr = zero
                    osumi = zero
c                   Need SUM( wvfcn(r)*exp(-ik.r) )
                    do  nr=1,nrbox
C eikr(2,irbox+nrbox*(ikloc-1))  (irbox=1,nrbox; ikloc=1,nkloc)
c wvfcns(1:2,1:nrbox(ibox),1:nkloc,1:norbshi(ishl,itypa(iatm)),1:ishl(iatm),1:natm,1:nbox)
                      osumr = osumr + wvfcns(1,ir0+nr)*eikr(1,kr0+nr)
     $                              + wvfcns(2,ir0+nr)*eikr(2,kr0+nr)
                      osumi = osumi - wvfcns(1,ir0+nr)*eikr(2,kr0+nr)
     $                              + wvfcns(2,ir0+nr)*eikr(1,kr0+nr)
                    enddo
                    ir0 = ir0 + nrbox
                    kr0 = kr0 + nrbox
c
                    coskri = phkdoti(1,ikloc)
                    sinkri = phkdoti(2,ikloc)
                    orbint(1,iorb,ik) = orbint(1,iorb,ik)
     $               - weight*(osumr*coskri-osumi*sinkri)
                    orbint(2,iorb,ik) = orbint(2,iorb,ik)
     $               - weight*(osumr*sinkri+osumi*coskri)
                  enddo
                endif
                irorb = irorb + nrkbox
  170         continue
            endif
c
            jnbox = knbox
            jrec = krec
            jatshl = katshl
            jrshl = krshl
            jorbat00 = 0
            do 181 jshl=1,ishl
              norbshj = 2*lshel(jshl,ityp) + 1
c
              jnbox = jnbox/2
              if( inrec*(jatshl/inrec).eq.jatshl )then
                jrec = jrec + 1
                jnbox = inboxs(jrec,ibox)
              endif
              jatshl = jatshl + 1
              jorbat0 = jorbat00
              jorbat00 = jorbat00 + norbshj
              if( EVEN(jnbox) ) goto 181
c
              joff = 0
              if( jshl.eq.ishl )then
                norbshj = 0
                joff = 1
              endif
c
              nalj = nala(jshl,ityp)
              almaxj = ala(nalj,jshl,ityp)
              do  jal=1,nalj-1
                alj = ala(jal,jshl,ityp)
                if( alj.gt.almaxj ) almaxj = alj
              enddo
c
              if( almaxi.lt.alfast .and. almaxj.lt.alfast )then
c               Both slow shells, skip
                jrshl = jrshl + norbshj*nrkbox
                goto 181
              endif
c
              irorb = irshl
              iorbat = iorbat0
              do 172 iorbsh=1,norbshi
                iorbat = iorbat + 1
c
                norbshj = norbshj + joff
c
                jrorb = jrshl
                jorbat = jorbat0
                do 171 jorbsh=1,norbshj
                  jorbat = jorbat + 1
c
                  ir0 = irorb
                  jr0 = jrorb
c
                  do 60 ikloc=1,nkloc
                    ik = nk0 + ikloc
                    ovlpr = zero
                    ovlpi = zero
                    if( ncplx.eq.1 )then
                      do  nr=1,nrbox
                        ovlpr= ovlpr + wf(ir0+nr)*wf(jr0+nr)
                      enddo
                    elseif( ncplx.eq.2 )then
                      do  nr=1,nrbox
c                       ovlp =+ wvfcns(ir0+nr)*CONJG(wvfcns(jr0+nr))
                        ovlpr= ovlpr + wvfcns(1,ir0+nr)*wvfcns(1,jr0+nr)
     $                               + wvfcns(2,ir0+nr)*wvfcns(2,jr0+nr)
                        ovlpi= ovlpi + wvfcns(2,ir0+nr)*wvfcns(1,jr0+nr)
     $                               - wvfcns(1,ir0+nr)*wvfcns(2,jr0+nr)
                      enddo
                    endif
                    ovlpr = weight*ovlpr
                    ovlpi = weight*ovlpi
c
c                   Subtract grid overlap to form "defect" overlap:
c
                    s1atom(iorbat,jorbat,ik,iatm) =
     $               s1atom(iorbat,jorbat,ik,iatm) - ovlpr
                    if( iorbat.ne.jorbat )then
                      s1atom(jorbat,iorbat,ik,iatm) =
     $                 s1atom(jorbat,iorbat,ik,iatm) - ovlpi
                    endif
c
                    ir0 = ir0 + nrbox
                    jr0 = jr0 + nrbox
   60             continue
c
                  jrorb = jrorb + nrkbox
  171           continue
c
                irorb = irorb + nrkbox
  172         continue
c
              jrshl = jrshl + norbshj*nrkbox
  181       continue
c
            irshl = irshl + norbshi*nrkbox
  182     continue
c
c         Close atom loop:
  200   continue
c
c       Close box loop:
  900   continue
        ir00 = ir00 + nobox*nrkbox
 1000 continue
c     Bypass from nkloc=0:
 1001 continue
c
      if( nprocl .gt. 1 )then
c       Need to sum results over all processors
        if( near2c .eq. 2 )then
          call MPREDUC8( masterl, lenvec2, orbint, otmp, icomml )
          if( iprocl.eq.masterl ) call DCOPY( lenvec2, otmp,1, orbint,1)
          call MPBCAST8( masterl, lenvec2, orbint, icomml )
        endif
        call MPREDUC8( masterl, lenvec1, s1atom, otmp, icomml )
        if( iprocl.eq.masterl ) call DCOPY( lenvec1, otmp,1, s1atom,1 )
        call MPBCAST8( masterl, lenvec1, s1atom, icomml )
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Kill those terms where the orbitals are slowly varying
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      i0 = 0
      do 1200 iatm=1,natm
        ityp = itypa(iatm)
c
        iorbat00 = 0
        do 1185 ishl=1,numshl(ityp)
          norbshi = 2*lshel(ishl,ityp) + 1
          iorbat0 = iorbat00
          iorbat00 = iorbat00 + norbshi
c
          nali = nala(ishl,ityp)
          almaxi = ala(nali,ishl,ityp)
          do  ial=1,nali-1
            ali = ala(ial,ishl,ityp)
            if( ali.gt.almaxi ) almaxi = ali
          enddo
c
          if( near2c.eq.2 .and. almaxi.lt.alfast )then
c           Slow orbital, blank this term
            do  ik=1,nk
              do  iorb=i0+1,i0+norbshi
                orbint(1,iorb,ik) = zero
                orbint(2,iorb,ik) = zero
              enddo
            enddo
          endif
          i0 = i0 + norbshi
c
          iorbat = iorbat0
          do 1175 iorbsh=1,norbshi
            iorbat = iorbat + 1
c
            jorbat00 = 0
            do 1180 jshl=1,ishl
              norbshj = 2*lshel(jshl,ityp) + 1
              jorbat0 = jorbat00
              jorbat00 = jorbat00 + norbshj
c
              nalj = nala(jshl,ityp)
              almaxj = ala(nalj,jshl,ityp)
              do  jal=1,nalj-1
                alj = ala(jal,jshl,ityp)
                if( alj.gt.almaxj ) almaxj = alj
              enddo
              if( jshl.eq.ishl ) norbshj = iorbsh
c
              if( almaxi.lt.alfast .and. almaxj.lt.alfast )then
c               Both slow orbitals, blank term
                jorbat = jorbat0
                do 1070 jorbsh=1,norbshj
                  jorbat = jorbat + 1
                  if( iorbat.eq.jorbat )then
                    do  ik=1,nk
                      s1atom(iorbat,jorbat,ik,iatm) = zero
                    enddo
                  else
                    do  ik=1,nk
                      s1atom(iorbat,jorbat,ik,iatm) = zero
                      s1atom(jorbat,iorbat,ik,iatm) = zero
                    enddo
                   endif
 1070           continue
                goto 1180
              endif
c
 1180       continue
c
 1175     continue
 1185   continue
c
c    Close atom loop:
c
 1200 continue
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
