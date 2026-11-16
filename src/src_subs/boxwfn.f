c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BOXWFN
c
c
      subroutine BOXWFN( ibox, ncplx,
     $ nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,ratm,rlat,coskr,sinkr,
     $ nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $ wvfcns, wf,
     $ r, sumg,rsq, rv1,rv2,rv3,rv4,rv5 )
c---------------------------------------------------------------
c Purpose: generate orbitals on regular space grid for a box
c
c Written: Peter A. Schultz, 28-May-1997, for v2.20 (boxes)
c
c Revision history:
c  18Jan09-PAS/2.63: k-parallel
c   7Mar02-PAS/2.52: ang constants extracted
c  21Jun01-PAS/2.48: replace STOPs
c  21Dec99-PAS/2.42: cosmetic clean-up
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns,wf
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  inboxs(nrec,ibox), boxctr(3)
      DIMENSION  r(3,*)
c output arrays:
      DIMENSION  wvfcns(2,*),wf(*)
c scratch arrays:
      DIMENSION  sumg(*),rsq(*)
      DIMENSION  rv1(*),rv2(*),rv3(*),rv4(*),rv5(*)
c local declarations:
      DIMENSION  ratom(3)
c
      DIMENSION  vnorm(5)
      DATA  zero, three / 0.d0, 3.d0 /
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
      ang4 = vnorm(3) / vnorm(4)
      ang5 = vnorm(3) / vnorm(5)
c
      nrkbox = nkloc*nrbox
c
      ir00 = 0
      inbox = 0
c
c    Clear out wave function array:
c
      nwf = ncplx*nkloc*nobox*nrbox
      if( ncplx.eq.1 )then
        call MKSZERO( nwf, wf(ir00+1) )
      else
        call MKSZERO( nwf, wvfcns(1,ir00+1) )
      endif
c
c    Center loop:
c
      irec = 0
      iatshl = 0
      irshl = ir00
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
        almina = alamin(ityp)
c
        do 800 ishl=1,numshl(ityp)
          norbshi = 2*lshel(ishl,ityp) + 1
          inbox = inbox/2
          if( inrec*(iatshl/inrec).eq.iatshl )then
            irec = irec + 1
            inbox = inboxs(irec,ibox)
          endif
          iatshl = iatshl + 1
          if( EVEN(inbox) ) goto 800
c
          lshl = lshel(ishl,ityp)
          nali = nala(ishl,ityp)
c
          almins = ala(1,ishl,ityp)
          do  ial=2,nali
            ali = ala(ial,ishl,ityp)
            if( ali.lt.almins ) almins = ali
          enddo
          almin = almina
          rngmax = SQRT( cutgrid / (almin +almins) )
c
c  Cycle through lattice vectors:
c
          do 900 ilat=1,nlat
c
            rngbox = zero
            do  jr=1,3
              ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
              rngbox = rngbox + (ratom(jr)-boxctr(jr))**2
            enddo
            rngbox = SQRT( rngbox ) - boxrad
            if( rngbox.gt.rngmax ) goto 900
c
            rsqbox = rngbox**2
            if( rngbox.lt.zero ) rsqbox = zero
c
            do  nr=1,nrbox
              rv1(nr) = r(1,nr) - ratom(1)
              rv2(nr) = r(2,nr) - ratom(2)
              rv3(nr) = r(3,nr) - ratom(3)
              rsq(nr) = rv1(nr)**2 + rv2(nr)**2 + rv3(nr)**2
              sumg(nr) = zero
            enddo
c
            ndo = 0
            do 70 ial=1,nali
              alf = ala(ial,ishl,ityp)
              twoalc = alf + almin
              rsqcut = cutgrid/twoalc
              if( rsqbox.gt.rsqcut ) goto 70
c
              coef = cala(ial,ishl,ityp) / vnorm(lshl+1)
              ndo = ndo + 1
              do  nr=1,nrbox
                sumg(nr) = sumg(nr) + coef*EXP( -alf*rsq(nr) )
              enddo
   70       continue
            if( ndo.eq.0 ) goto 900
c
            if( ncplx.eq.1 )then
c
c >>>>> Real orbitals/molecular systems
c
              ir0 = irshl
              if( lshl.eq.0 )then
c
                do  nr=1,nrbox
                  wf(ir0+nr) = wf(ir0+nr) + sumg(nr)
                enddo
c
              elseif( lshl.eq.1 )then
c
                nr1 = ir0
                nr2 = nr1 + nrkbox
                nr3 = nr2 + nrkbox
                do  nr=1,nrbox
                  wf(nr1+nr) = wf(nr1+nr) + sumg(nr)*rv1(nr)
                  wf(nr2+nr) = wf(nr2+nr) + sumg(nr)*rv2(nr)
                  wf(nr3+nr) = wf(nr3+nr) + sumg(nr)*rv3(nr)
                enddo
c
              elseif( lshl.eq.2 )then
c
                nr1 = ir0
                nr2 = nr1 + nrkbox
                nr3 = nr2 + nrkbox
                nr4 = nr3 + nrkbox
                nr5 = nr4 + nrkbox
                do  nr=1,nrbox
                  wf(nr1+nr) = wf(nr1+nr) + sumg(nr) * rv1(nr)*rv2(nr)
                  wf(nr2+nr) = wf(nr2+nr) + sumg(nr) * rv2(nr)*rv3(nr)
                  wf(nr3+nr) = wf(nr3+nr) + sumg(nr) * rv3(nr)*rv1(nr)
                  wf(nr4+nr) = wf(nr4+nr) + sumg(nr) *
     $             (rv1(nr)**2-rv2(nr)**2)*ang4
                  wf(nr5+nr) = wf(nr5+nr) + sumg(nr) *
     $             (three*rv3(nr)**2-rsq(nr))*ang5
                enddo
c
              else
c
                call STOPXERR( 'boxwfn-l/lshl is too high' )
c
              endif
c
            else
c
c >>>>> Complex orbitals/periodic systems
c
c start out with some arithmetic to accelerate things
              if( lshl.eq.1 )then
                do  nr=1,nrbox
                  rv1(nr) = sumg(nr)*rv1(nr)
                  rv2(nr) = sumg(nr)*rv2(nr)
                  rv3(nr) = sumg(nr)*rv3(nr)
                enddo
              elseif( lshl.eq.2 )then
c the order in which I do this critical for memory mgmt
                do  nr=1,nrbox
                  rsq(nr) = sumg(nr)*(three*rv3(nr)**2-rsq(nr))*ang5
                  rv4(nr) = sumg(nr)*(rv1(nr)**2-rv2(nr)**2)*ang4
                  rv5(nr) = sumg(nr)*rv1(nr)*rv2(nr)
                  rv2(nr) = sumg(nr)*rv2(nr)*rv3(nr)
                  rv3(nr) = sumg(nr)*rv3(nr)*rv1(nr)
                enddo
              endif
c
              ir0 = irshl
              do 400 ikloc=1,nkloc
                ik = nk0 + ikloc
                facr = coskr(ik,ilat)
                faci = sinkr(ik,ilat)
c
                if( lshl.eq.0 )then
c
                  do  nr=1,nrbox
                    wvfcns(1,ir0+nr) = wvfcns(1,ir0+nr) + facr*sumg(nr)
                    wvfcns(2,ir0+nr) = wvfcns(2,ir0+nr) + faci*sumg(nr)
                  enddo
c
                elseif( lshl.eq.1 )then
c
                  nr1 = ir0
                  nr2 = nr1 + nrkbox
                  nr3 = nr2 + nrkbox
                  do  nr=1,nrbox
                    wvfcns(1,nr1+nr) = wvfcns(1,nr1+nr) + facr*rv1(nr)
                    wvfcns(2,nr1+nr) = wvfcns(2,nr1+nr) + faci*rv1(nr)
                    wvfcns(1,nr2+nr) = wvfcns(1,nr2+nr) + facr*rv2(nr)
                    wvfcns(2,nr2+nr) = wvfcns(2,nr2+nr) + faci*rv2(nr)
                    wvfcns(1,nr3+nr) = wvfcns(1,nr3+nr) + facr*rv3(nr)
                    wvfcns(2,nr3+nr) = wvfcns(2,nr3+nr) + faci*rv3(nr)
                  enddo
c
                elseif( lshl.eq.2 )then
c
                  nr1 = ir0
                  nr2 = nr1 + nrkbox
                  nr3 = nr2 + nrkbox
                  nr4 = nr3 + nrkbox
                  nr5 = nr4 + nrkbox
                  do  nr=1,nrbox
                    wvfcns(1,nr1+nr) = wvfcns(1,nr1+nr) + facr*rv5(nr)
                    wvfcns(2,nr1+nr) = wvfcns(2,nr1+nr) + faci*rv5(nr)
                    wvfcns(1,nr2+nr) = wvfcns(1,nr2+nr) + facr*rv2(nr)
                    wvfcns(2,nr2+nr) = wvfcns(2,nr2+nr) + faci*rv2(nr)
                    wvfcns(1,nr3+nr) = wvfcns(1,nr3+nr) + facr*rv3(nr)
                    wvfcns(2,nr3+nr) = wvfcns(2,nr3+nr) + faci*rv3(nr)
                    wvfcns(1,nr4+nr) = wvfcns(1,nr4+nr) + facr*rv4(nr)
                    wvfcns(2,nr4+nr) = wvfcns(2,nr4+nr) + faci*rv4(nr)
                    wvfcns(1,nr5+nr) = wvfcns(1,nr5+nr) + facr*rsq(nr)
                    wvfcns(2,nr5+nr) = wvfcns(2,nr5+nr) + faci*rsq(nr)
                  enddo
c
                else
c
                  call STOPXERR( 'boxwfn-l/lshl is too high' )
c
                endif
c
c  Close bloch vector loop:
c
                ir0 = ir0 + nrbox
  400         continue
c
            endif
c
c  Close lattice vector loop:
c
  900     continue
c
c  Close shell/atom loops:
c
          irshl = irshl + norbshi*nrkbox
  800   continue
 1000 continue
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
