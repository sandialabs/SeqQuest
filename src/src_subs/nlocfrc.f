c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLOCFRC
c
c
      subroutine NLOCFRC( nafrc,iforce, ndim, frcnl,strnl,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,  nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl,  ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,vpsrad, wkrad,
     $ dmat,    gntcomb, angsav, frctp )
c---------------------------------------------------------------
c Purpose: calculate the force and stress due to the non-local
c          part of the electron-effective nuclear potential
c          interaction
c
c Written: P.J. Feibelman, extensively revised by P.A. Schultz
c
c Revision history:
c  15Oct08-PAS/2.62: stripe-parallel
c  11May07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c   6May03-DMC: Add frctp to argument list to provide a local copy
c               of the frcnl array, for parallelization.
c   7Mar02-PAS/2.52: ang constants extracted
c  23Jul01-PAS/2.49: compact dmat
c  21Jun01-PAS/2.48: replace STOPs
c  22May01-PAS/2.47: condense radial scratch space
c  15Mar00-PAS/2.43: rename routines to "nlXXXX"
c  15Feb00-PAS/2.43: opt single atom star angular calc.
c  15Oct99-PAS/2.41: Raczkowski optimizations integrated/cleaned up
c   5Oct99-DBR/2.41: optimization: angular terms
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  31Aug99-PAS/2.39: consistent cutoffs - split stars
c  17Dec98-PAS/2.30: clean force/stress arrays
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  29Aug98-PAS/2.25: full rework of bottom loops: finally a full
c    separation of angular/radial terms, gaussian terms condensed;
c    algebra simplified; lower triangle only; renamed NLOCFRC
c  11Aug98-PAS/2.24: fix cutoff bug, full proper cutoff scaling
c  21Jul98-PAS/2.23: clean up cutoffs, extend I:N terms
c  22May98-PAS/2.22: bugfix
c  28Apr97-PAS/2.17: naming conventions
c  11Nov94-PAS/2.12: setup radial integral routine installed
c  27May94-PAS/2.11: big cleanup/ optimize
c   9Mar94-PAS/2.10: fix nlat, j-index bugs; match to setup routine;
c               purge complex declarations; partial cleanup
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c  Allows for non-locality up to l=lmxnlp1(ktyp)-1
c
c  Pseudopotential of atom type ktyp is assumed to fall to zero
c  at large radius as exp(-r**2*almnnl(ktyp)).
c  almnnl(ktyp) is taken from the fits of Bachelet, et al.
c
c  Because of the nature of the non-local potential, the scaling of
c  the matrix elements is not the same as for the local potential
c  integrals (as the norm of the 3-center overlap gaussian). Getting
c  an effective and efficient cutoff scheme is a bit tricky, and the
c  code implemented in the following reflects this.  The scheme
c  implemented here is overly inclusive. It does incorporate scaling
c  of the magnitude of the matrix element properly, but, at best,
c  takes only partial advantage of the "importance" scaling of the
c  matrix elements, i.e., the density matrix or overlap.  Whether an
c  efficient and effective scheme can be fashioned to incorporate
c  this importance evaluation is left as an exercise for the reader.
c    - 11Aug98-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  frcnl(3,natm), strnl(3,3)
c Input arrays:
      DIMENSION  iforce(natm)
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  gntcomb(*)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),vpsrad(nrd,4,ntyp)
      DIMENSION  dmat(nk,*)
c scratch arrays: (missing wkrad dim=nald**2)
      DIMENSION  wkrad(nrd,*)
      DIMENSION  angsav(5*5*5*5*3*81)
c parallel scratch:
      DIMENSION  frctp(3,natm), strtp(3,3)
c
c Local declarations:
      LOGICAL    iatmi, diagshl
      PARAMETER  (lintd=6*6*18)
      DIMENSION  vnlint(lintd)
      DIMENSION  vnln(lintd),vnli(lintd),vnlj(lintd)
      DIMENSION  engri(3),engrj(3)
      DIMENSION  rnuc(3), ri(3),rj(3), risav(3)
c  need neqlatd to be pretty big for bulk
      PARAMETER  (neqlatd=72)
      DIMENSION  ristar(3,neqlatd)
      DIMENSION  aj(9,9),ai(9,9,neqlatd), ilateq(neqlatd)
      DIMENSION  ylmj(2,25),ylmi(2,25,neqlatd)
      DIMENSION  ylmj_t(2,5,5),ylmi_t(2,5,5,neqlatd)
      DIMENSION  gradaj(3,9,9),gradai(3,9,9,neqlatd),
     $           grylmj(2,3,25),grylmi(2,3,25,neqlatd)
      DIMENSION  grylmj_t(3,2,5,5),grylmi_t(3,2,5,5,neqlatd)
c
      DIMENSION  npprl(5)
      DATA       npprl / 1,2,2,2,2 /
c
      DIMENSION  vnorm(9)
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'FDEV/NLOCFRC: stripe-parallel'
      call FLGETIWR( IWR )
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxf = cutgwfn
      emaxsl = cutslow
      emaxii = cutii
      emxslo = MAX( emaxsl, emaxii )
c
      rt4pi = two*rtpi
c
      mxeqil = 0
      call MKZERO( 3*natm, frcnl )
      call MKZERO( 9, strnl )
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c     
c  Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c     
c  Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,   
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
c    Nucleus loop:
c
      do 3000 nuc=1,natm
c
Cc     Each processor will do every nprocs element in this loop.
C        if( MOD( nuc, nprocs ) .ne. iproc ) goto 3000
c
      ktyp = itypa(nuc)
      if( znuc(ktyp).eq.zero ) goto 3000
c
      lmxp1 = lmxnlp1(ktyp)
      if( lmxp1.lt.1 ) goto 3000
c
      alaminn = almnnl(ktyp)
      lmxnl = lmxp1 - 1
c
c    Nucleus lattice loop:
c
      do 2900 klat=1,nlat
c
      do  ir=1,3
        rnuc(ir) = ratm(ir,nuc) + rlat(ir,klat)
      enddo
c
c     Atom j loop:
c
      jjnext = 0
      jnext = myj0
      do 2000 jatm=myatm1,myatm2
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        do  ir=1,3
          rj(ir) = ratm(ir,jatm) - rnuc(ir)
        enddo
        rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
        rjmag = SQRT( rjsquar )
c
        alaminj = alamin(jtyp)
        emax = emxslo
        rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c       Test overlap of Atom j and short-range pot due to "nuc"
        if( rjsquar.gt.rjnmax ) goto 2000
c
        ljmax = 0
        do  jshl=1,nshlj
          if( lshel(jshl,jtyp).gt.ljmax ) ljmax = lshel(jshl,jtyp)
        enddo
        call AINITF( aj,gradaj )
        if( rjmag .lt. 1.d-12 )then
          lpowj = 0
          ljmxp1 = 1
          ylmj(1,1) = one / rt4pi
          call MKZERO( 6, grylmj )
c
          ylmj_t(1,1,1) = ylmj(1,1)
          call MKZERO( 6, grylmj_t )
c
        else
          lpowj = 1
          ljmxp1 = ljmax + lmxp1
          call AOFRF( aj,gradaj, rj )
          call DYLM( ylmj,grylmj, ljmax+lmxnl, rj )
c
          jlm = 0
          do  mmj=1,ljmxp1
            do  llj=mmj,ljmxp1
              jlm = jlm + 1
              do  jpprl=1,npprl(mmj)
                ylmj_t(jpprl,mmj,llj) = ylmj(jpprl,jlm)
                grylmj_t(1,jpprl,mmj,llj) = grylmj(jpprl,1,jlm)
                grylmj_t(2,jpprl,mmj,llj) = grylmj(jpprl,2,jlm)
                grylmj_t(3,jpprl,mmj,llj) = grylmj(jpprl,3,jlm)
              enddo
            enddo
          enddo
c
        endif
c
c       Atom j shell loop:
c
        jorbatnxt = 0
        j = jstart
        do 1800 jshl=1,nshlj
          nalj = nala(jshl,jtyp)
          lshlj = lshel(jshl,jtyp)
          norbshj = 2*lshlj + 1
          lj0 = lshlj*lshlj
c
          j0 = j
          jorbat0 = jorbatnxt
          jorbatnxt = jorbat0 + norbshj
c
          alminj = ala(1,jshl,jtyp)
          almaxj = ala(nalj,jshl,jtyp)
          do  jal=1,nalj
            alj = ala(jal,jshl,jtyp)
            if( alj.lt.alminj ) alminj = alj
            if( alj.gt.almaxj ) almaxj = alj
          enddo
          alaminj = alminj
c
          emax = emxslo
          rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c         Test whether this radial function will contribute:
          if( rjsquar.gt.rjnmax ) goto 1799
c
c         Atom i loop:
c
          ijnext = jjstart
          inext = jstart
          do 1000 iatm=jatm,natm
            ityp = itypa(iatm)
            nshli = numshl(ityp)
            norbati = norba(ityp)
            istart = inext
            inext = istart + norbati
c
            matij = norbati*norbatj
            ijstart = ijnext
            ijnext = ijstart + matij*ncplx
c
c           Following helps ensure only i.ge.j done:
            iorbatn0 = 0
            ishl1 = 1
            if( iatm.eq.jatm )then
              ishl1 = jshl
              do  ishl=1,ishl1-1
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                istart = istart + norbshi
                iorbat0 = iorbatn0
                iorbatn0 = iorbat0 + norbshi
              enddo
            endif
c
            alamini = alamin(ityp)
c           Criterion for including contribution of Atom i:
            ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
            rijnmax = (emaxsl-ex)*(alaminj+alamini+alaminn) /
     $                          ( (alaminj+alaminn)*alamini )
            riinmax = (emaxii-ex)*(alaminj+alamini+alaminn) /
     $                          ( (alaminj+alaminn)*alamini )
c
            if( iatm.eq.jatm )then
              emaxij = emaxii
              rijnmx0 = riinmax
            else
              emaxij = emaxsl
              rijnmx0 = rijnmax
            endif
            rinmax = emaxij*(alamini+alaminn) / (alamini*alaminn)
            rinmx1 = rinmax
            rjovlp = rjsquar * ( alaminj / (alaminj+alaminn) )**2
c
            limax = 0
            do  ishl=ishl1,nshli
              if( lshel(ishl,ityp).gt.limax ) limax = lshel(ishl,ityp)
            enddo
c
c           Find nearest iatm image:
            risqnxt = 9.d9
            do 90 ilat=1,nlat2c
              do  ir=1,3
                ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
              enddo
              risquar = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
              if( risquar.lt.risqnxt )then
                risqnxt = risquar
                ilatnxt = ilat
              endif
   90       continue
c
c           Lattice vector loop for Atom i:
c
c           Need to preset initial ilat now.
c           I:I term, if it appears in star, appears in ilat=1.
c           If it is found, exit "star search", do the integrals,
c           and return with ilat1=2 to create non-I:I part of star
            ilat1 = 1
c
            do 900 ilatx=1,nlat2c
              risquar = risqnxt
              risqnxt = 9.d9
c
c             First test that iatm overlaps nuc ...
              if( risquar.gt.rinmx1 ) goto 1000
c
c               ... and test for 3-center overlap as well:
c                   (I align the ri and rj vectors for this test)
              rimag = SQRT( risquar )
              if( risquar.gt.rjovlp )then
                rinmx2 = rijnmx0 * ( (alaminj+alaminn) /
     $           (alaminj+alaminn-alaminj*rjmag/rimag) )**2
                if( risquar.gt.rinmx2 ) goto 1000
              endif
c
c             Align the iatm test vector along the jatm vector
              if( rjmag .gt. 1.d-6 )then
                scale = rimag/rjmag
                do  ir=1,3
                  risav(ir) = scale*rj(ir)
                enddo
              else
                do  ir=1,3
                  risav(ir) = ratm(ir,iatm)+rlat(ir,ilatnxt) - rnuc(ir)
                enddo
              endif
c
              f1 = alaminj / (alaminj+alaminn)
              rijn0 = zero
              do  ir=1,3
                rijn0 = rijn0 + (f1*rj(ir)-risav(ir))**2
              enddo
c
              iatmi = .false.
              ieqlat = 0
              do 190 ilat=ilat1,nlat2c
c
                do  ir=1,3
                  ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
                enddo
                risqtst = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
c               Compile orientation info of atoms such that risq=risquar
c
                if( risqtst-risquar .gt. 1.d-6
     $           .or. ( iatmi .and. risqtst.ge.risquar ) )then
c
c                 This risq is bigger, check if it is next biggest ...
                  if( risqtst.lt.risqnxt )then
                    risqnxt = risqtst
                    ilatnxt = ilat
                  endif
c
                elseif( ABS( risqtst-risquar ) .le. 1.d-6 )then
c
c                 This risq is within current star, process it
c
                  if( ilat.eq.1 .and. iatm.eq.jatm )then
                    if( rijn0.gt.riinmax ) goto 190
                    iatmi = .true.
                  else
c                   Check tighter cutoff ...
                    if( rijn0.gt.rijnmax ) goto 190
                  endif
c
                  ieqlat = ieqlat + 1
                  if( ieqlat.gt.neqlatd ) goto 190
c
                  ilateq(ieqlat) = ilat
                  do  ir=1,3
                    ristar(ir,ieqlat) = ri(ir)
                  enddo
c
                  call AINITF( ai(1,1,ieqlat),gradai(1,1,1,ieqlat) )
                  if( risqtst .lt. 1.d-12 )then
                    limxp1 = 1
                    ylmi(1,1,ieqlat) = one / rt4pi
                    call MKZERO( 6, grylmi(1,1,1,ieqlat) )
c
                    ylmi_t(1,1,1,ieqlat) = ylmi(1,1,ieqlat)
                    call MKZERO( 6, grylmi_t(1,1,1,1,ieqlat) )
c
                  else
                    limxp1 = limax + lmxp1
                    call AOFRF( ai(1,1,ieqlat),gradai(1,1,1,ieqlat), ri)
                    call DYLM( ylmi(1,1,ieqlat),grylmi(1,1,1,ieqlat),
     $                         limax+lmxnl, ri )
c
                    ilm = 0
                    do  mmi=1,limxp1
                      do  lli=mmi,limxp1
                        ilm = ilm + 1
                        do  ipprl=1,npprl(mmi)
                          ylmi_t(ipprl,mmi,lli,ieqlat) =
     $                     ylmi(ipprl,ilm,ieqlat)
                          grylmi_t(1,ipprl,mmi,lli,ieqlat) =
     $                     grylmi(ipprl,1,ilm,ieqlat)
                          grylmi_t(2,ipprl,mmi,lli,ieqlat) =
     $                     grylmi(ipprl,2,ilm,ieqlat)
                          grylmi_t(3,ipprl,mmi,lli,ieqlat) =
     $                     grylmi(ipprl,3,ilm,ieqlat)
                        enddo
                      enddo
                    enddo
c
                  endif
c
                endif
c
  190         continue
c
              if( iatmi )then
                ilat1 = 2
                if( rijn0.gt.riinmax ) goto 900
              else
                ilat1 = 1
                if( rijn0.gt.rijnmax ) goto 900
              endif
c
c             Check how many atoms in star we picked up
              neqilat = ieqlat
              if( neqilat.gt.mxeqil ) mxeqil = neqilat
              if( neqilat.gt.neqlatd )then
                write(IWR,*) '>>>>> neqilat,neqlatd=',neqilat,neqlatd
                call STOPXERR( 'neqlatd - exceeded in NLOCFRC' )
              endif
              if( neqilat.le.0 ) goto 900
c
c             Find largest l for Atom i with alpha that will contribute:
              ilmax = -1
              do 280 ishl=ishl1,nshli
                nali = nala(ishl,ityp)
                lshli = lshel(ishl,ityp)
c
                almini = ala(1,ishl,ityp)
                almaxi = ala(nali,ishl,ityp)
                do  ial=1,nali
                  ali = ala(ial,ishl,ityp)
                  if( ali.lt.almini ) almini = ali
                  if( ali.gt.almaxi ) almaxi = ali
                enddo
c
                emax = emaxf
                if( almaxi.lt.alfast .and. almaxj.lt.alfast )emax=emaxsl
                if( iatmi ) emax = emaxii
c
                ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
                rijnmx = (emax-ex)*(alaminj+almini+alaminn) /
     $                           ( (alaminj+alaminn)*almini )
                if( rijn0.le.rijnmx )then
                  if( lshli.gt.ilmax ) ilmax = lshli
                endif
  280         continue
c
              if( ilmax.lt.0 ) goto 900
c
              rimag = SQRT( risquar )
c
              lisq = limxp1*limxp1
              ljsq = ljmxp1*ljmxp1
c
              maxlj = (lshlj+1)*(lshlj+1)
              maxli = (ilmax+1)*(ilmax+1)
c
c             Carve up space in angsav() for arrays in angsetf ...
c
              iangs = 1
c              ... angsav() needs:
              nangs = limxp1*ljmxp1*7*lmxp1*9*9
              iangt = iangs + nangs
c              ... angtmp() needs:
              nangt = lisq*ljsq
              igntj = iangt + nangt
c              ... gntcj() (complex, use 4* real space for safety) needs:
              ngntj = 4*(ljsq*maxlj*lmxp1*lmxp1)
              ignti = igntj + ngntj
c             ...  cgntci() (also complex) needs:
              ngnti = 4*(lisq*maxli*lmxp1*lmxp1)
c
              if( neqilat.eq.1 )then
c
c             If I have only one atom in a star, then compute the angular
c             stuff up here only once, rather than repeat for every shell
c             If more than one atom per star, I either have to duplicate
c             angsav for each element in the star, or, as done here, compute
c             it below for each star. 15Feb00-PAS/2.43
c
                ieqlat = 1
                call NLANGF(
     $           maxlj,maxli,ljmxp1,limxp1,lmxp1,ljsq,lisq,
     $           ylmj_t,ylmi_t(1,1,1,ieqlat),
     $           grylmj_t,grylmi_t(1,1,1,1,ieqlat),gntcomb,
     $           angsav,angsav(iangt),angsav(igntj),angsav(ignti) )
              endif
c
c             Atom i shell loop:
c
              iorbatnxt = iorbatn0
              i = istart
              do 800 ishl=ishl1,nshli
                nali = nala(ishl,ityp)
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                li0 = lshli*lshli
c
                i0 = i
                iorbat0 = iorbatnxt
                iorbatnxt = iorbat0 + norbshi
c
                diagshl = .false.
                if( i0.eq.j0 ) diagshl = .true.
c
                almini = ala(1,ishl,ityp)
                almaxi = ala(nali,ishl,ityp)
                do  ial=1,nali
                  ali = ala(ial,ishl,ityp)
                  if( ali.lt.almini ) almini = ali
                  if( ali.gt.almaxi ) almaxi = ali
                enddo
c
                emax = emaxf
                if( almaxi.lt.alfast .and. almaxj.lt.alfast )emax=emaxsl
                if( iatmi ) emax = emaxii
c
                lpowmx = lpowj
                lpowi = 0
                if( rimag .ge. 1.d-12 )then
                  lpowi = 1
                  lpowmx = 1
                endif
                limxpd = limxp1 + lpowi
                ljmxpd = ljmxp1 + lpowj
c
                lpowmx = lshli + lshlj + 1 + lpowmx
                lstmx = lpowmx*lmxp1
                lijmx = limxpd*ljmxpd*lstmx
c
                do  lij=1,lijmx
                  vnln(lij) = zero
                  vnli(lij) = zero
                  vnlj(lij) = zero
                enddo
c
c               Alpha loops:
c
                nx = 0
                do 320 jal=1,nalj
                  alj = ala(jal,jshl,jtyp)
                  calj = cala(jal,jshl,jtyp)
                  twoalj = alj + alj
c
c                 Test whether this alpha contributes:
                  rjnmx = emax*(alj+alaminn) / (alj*alaminn)
                  if( rjsquar.gt.rjnmx ) goto 320
c
                  fj1 = alj / (alj+alaminn)
                  rijn = zero
                  do  ir=1,3
                    rijn = rijn + (fj1*rj(ir)-risav(ir))**2
                  enddo
c
                  do 310 ial=1,nali
                    ali = ala(ial,ishl,ityp)
                    cali = cala(ial,ishl,ityp)
                    twoali = ali + ali
c
c                   Test whether this pair of alphas contributes:
                    ex = rjsquar*alj*alaminn / (alj+alaminn)
                    rijnmx = (emax-ex)*(alj+ali+alaminn) /
     $                               ( (alj+alaminn)*ali )
                    if( rijn.gt.rijnmx ) goto 310
c
                    nx = nx + 1
c
                    alsum = ali + alj
                    exparg = ali*risquar + alj*rjsquar
                    bargi = twoali*rimag
                    bargj = twoalj*rjmag
                    bargs = bargi + bargj
c
c                   Evaluate radial integrals
c
                    call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $               exparg,bargs,alsum,radmsh(1,ktyp), wkrad )
                    nr0 = nrmn - 1
                    nrlen = nrmx - nrmn + 1
c
                    call NLRAD( vnlint, alsum,bargj,bargi,bargs,
     $               exparg, ljmxpd,limxpd,lmxp1,lstmx,lpowmx,
     $               nrlen,nrd,nr0, radmsh(1,ktyp),vpsrad(1,1,ktyp),
     $               wkrad(1,1), wkrad(1,8), wkrad(1,15),
     $               wkrad(1,22),wkrad(1,23),wkrad(1,24),
     $               wkrad(1,25),wkrad(1,26),wkrad(1,27) )
c
c                   Sum into running totals:
                    cnln = cali*calj
                    cnli = cnln*twoali
                    cnlj = cnln*twoalj
                    do  lij=1,lijmx
                      vnln(lij) = vnln(lij) + cnln*vnlint(lij)
                      vnli(lij) = vnli(lij) + cnli*vnlint(lij)
                      vnlj(lij) = vnlj(lij) + cnlj*vnlint(lij)
                    enddo
c
c                   Close alpha loops:
  310             continue
  320           continue
c
                if( nx.eq.0 ) goto 799
c
c               Loop over star of Atom j:
c
                do 490 ieqlat=1,neqilat
c
                  ilat = ilateq(ieqlat)
                  do  ir=1,3
                    ri(ir) = ristar(ir,ieqlat)
                  enddo
c
                  if( neqilat.gt.1 )then
                    call NLANGF(
     $               maxlj,maxli,ljmxp1,limxp1,lmxp1,ljsq,lisq,
     $               ylmj_t,ylmi_t(1,1,1,ieqlat),
     $               grylmj_t,grylmi_t(1,1,1,1,ieqlat),gntcomb,
     $               angsav,angsav(iangt),angsav(igntj),angsav(ignti) )
                  endif
c
c                 Angular loops:
c
                  ij0 = ijstart + jorbat0*norbati + iorbat0
                  ij0 = ij0 - norbati
                  mli1 = 1
                  do 420 mlj=1,norbshj
                    lj = lj0 + mlj
                    jorbat = jorbat0 + mlj
                    j = j0 + mlj
c
                    ij0 = ij0 + norbati
                    if( diagshl ) mli1 = mlj
                    do 410 mli=mli1,norbshi
                      li = li0 + mli
                      iorbat = iorbat0 + mli
                      i = i0 + mli
                      ij = ij0 + mli
c
                      denlilj = vnorm(li)*vnorm(lj)
c
c                     Do sum over density matrix:
                      dmsum = zero
                      if( i.eq.j )then
                        do  k=1,nk
C                          dmsum = dmsum + coskr(k,ilat)*dmat(i,i,k)
                          dmsum = dmsum + coskr(k,ilat)*dmat(k,ij)
                        enddo
                      else
c                       i.ne.j
                        do  k=1,nk
C                          dmsum = dmsum + coskr(k,ilat)*dmat(i,j,k)
                          dmsum = dmsum + coskr(k,ilat)*dmat(k,ij)
                        enddo
                        if( ncplx.eq.2 )then
                          ijc = ij + matij
                          do  k=1,nk
C                            dmsum = dmsum + sinkr(k,ilat)*dmat(j,i,k)
                            dmsum = dmsum + sinkr(k,ilat)*dmat(k,ijc)
                          enddo
                        endif
c                       Double because only doing one triangle
                        dmsum = dmsum + dmsum
                      endif
                      dmsum = dmsum / denlilj
c
                      call NLINTF( enprj,enpri, engrj,engri,
     $                 vnln,vnlj,vnli, rjmag,rimag, lj,li,
     $                 ljmxp1,limxp1,lmxp1,lpowmx,lpowj,lpowi,lstmx,
     $                 aj(1,lj),    gradaj(1,1,lj),
     $                 ai(1,li,ieqlat),gradai(1,1,li,ieqlat), angsav )
c
                      fxi = dmsum*( engri(1) + enpri*ri(1) )
                      fyi = dmsum*( engri(2) + enpri*ri(2) )
                      fzi = dmsum*( engri(3) + enpri*ri(3) )
                      fxj = dmsum*( engrj(1) + enprj*rj(1) )
                      fyj = dmsum*( engrj(2) + enprj*rj(2) )
                      fzj = dmsum*( engrj(3) + enprj*rj(3) )
c
                      frcnl(1,iatm) = frcnl(1,iatm) + fxi
                      frcnl(2,iatm) = frcnl(2,iatm) + fyi
                      frcnl(3,iatm) = frcnl(3,iatm) + fzi
c
                      frcnl(1,jatm) = frcnl(1,jatm) + fxj
                      frcnl(2,jatm) = frcnl(2,jatm) + fyj
                      frcnl(3,jatm) = frcnl(3,jatm) + fzj
c
                      frcnl(1,nuc) = frcnl(1,nuc) - fxi - fxj
                      frcnl(2,nuc) = frcnl(2,nuc) - fyi - fyj
                      frcnl(3,nuc) = frcnl(3,nuc) - fzi - fzj
c
                      strnl(1,1) = strnl(1,1) - fxi*ri(1) - fxj*rj(1)
                      strnl(2,1) = strnl(2,1) - fyi*ri(1) - fyj*rj(1)
                      strnl(3,1) = strnl(3,1) - fzi*ri(1) - fzj*rj(1)
                      strnl(1,2) = strnl(1,2) - fxi*ri(2) - fxj*rj(2)
                      strnl(2,2) = strnl(2,2) - fyi*ri(2) - fyj*rj(2)
                      strnl(3,2) = strnl(3,2) - fzi*ri(2) - fzj*rj(2)
                      strnl(1,3) = strnl(1,3) - fxi*ri(3) - fxj*rj(3)
                      strnl(2,3) = strnl(2,3) - fyi*ri(3) - fyj*rj(3)
                      strnl(3,3) = strnl(3,3) - fzi*ri(3) - fzj*rj(3)
c
c                     Close angular loops:
  410               continue
  420             continue
c
c                 Close star of iatm loop:
  490           continue
c
  799           continue
c               Close Atom i shell loop:
                i = i0 + norbshi
  800         continue
c
c              Close Atom i lattice loop:
  900       continue
c
c           Close Atom i loop:
 1000     continue
c
 1799     continue
c         Close Atom j shell loop:
          j = j0 + norbshj
 1800   continue
c
c       Close Atom j loop:
 2000 continue
c
c       Close Nucleus lattice loop:
c
 2900 continue
c
c       Close Nucleus loop:
c
 3000 continue
c
      if( nprocs .gt. 1 )then
c       Each local processor now has forces over its subset.
c       Need to sum these to a single value, on the master proc (0).
c       But for now, we will Broadcast total to all processors.
c
        nreduc = 3*natm
c
        call MPREDUC8( master, nreduc, frcnl , frctp , icomm )
C        call MPI_Reduce( frcnl , frctp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocfrc - Reduce frcnl' )
        call DCOPY( nreduc, frctp,1, frcnl,1 )
c
        call MPBCAST8( master, nreduc, frcnl , icomm )
C        call MPI_Bcast( frcnl , nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocfrc - Bcast frcnl' )
c
        nreduc = 9
        call MPREDUC8( master, nreduc, strnl , strtp , icomm )
C        call MPI_Reduce( strnl , strtp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocfrc - Reduce strnl' )
        call DCOPY( nreduc, strtp,1, strnl,1 )
c
        call MPBCAST8( master, nreduc, strnl , icomm )
C        call MPI_Bcast( strnl , nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocfrc - Bcast strnl' )
c
c       End tp
      endif
c
c    That's all Folks!
c
      RETURN
      END
