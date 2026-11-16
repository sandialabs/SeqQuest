c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PNLOCFRC
c
c
      subroutine PNLOCFRC( nafrc,iforce, ndim, frcnl,strnl,
     $ norb,nk, natm,ntyp,nshld,nald, nlat,  nlat2c,
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
c output arrays:
      DIMENSION  frcnl(3,natm), strnl(3,3)
c input arrays:
      DIMENSION  iforce(natm)
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  gntcomb(*)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),vpsrad(nrd,4,ntyp)
      DIMENSION  dmat(norb,norb,nk)
c scratch arrays: (missing dim=nald**2)
      DIMENSION  wkrad(nrd,*)
      DIMENSION  angsav(5*5*5*5*3*81)
ctp: task-parallel scratch:
      DIMENSION  frctp(3,natm), strtp(3,3)
c local declarations:
      PARAMETER  (lintd=6*6*18)
      DIMENSION  vnlint(lintd)
      DIMENSION  vnln(lintd),vnli(lintd),vnlj(lintd)
      DIMENSION  engri(3),engrj(3)
      LOGICAL    iatmi
      DIMENSION  rnuc(3), ri(3),rj(3), rjsav(3)
c  need neqlatd to be pretty big for bulk
      PARAMETER  (neqlatd=72)
      DIMENSION  rjstar(3,neqlatd)
      DIMENSION  ai(9,9),aj(9,9,neqlatd), lateqj(neqlatd)
      DIMENSION  ylmi(2,25),ylmj(2,25,neqlatd)
      DIMENSION  ylmi_t(2,5,5),ylmj_t(2,5,5,neqlatd)
      DIMENSION  gradai(3,9,9),gradaj(3,9,9,neqlatd),
     $           grylmi(2,3,25),grylmj(2,3,25,neqlatd)
      DIMENSION  grylmi_t(3,2,5,5),grylmj_t(3,2,5,5,neqlatd)
c
      DIMENSION  ime(5)
      DATA       ime / 1,2,2,2,2 /
c
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
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
      mxeqjl = 0
      call MKZERO( 3*natm, frcnl )
      call MKZERO( 9, strnl )
c
ctp:
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c    Nucleus loop:
c
      do 3000 nuc=1,natm
c
ctp:
c     Each processor will do every nprocs element in this loop.
        if( MOD( nuc, nprocs ) .ne. iproc ) goto 3000
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
c    Atom i loop:
c
      inext = 0
      do 2000 iatm=1,natm
        ityp = itypa(iatm)
        nshli = numshl(ityp)
        istart = inext
        inext = istart + norba(ityp)
c
        do  ir=1,3
          ri(ir) = ratm(ir,iatm) - rnuc(ir)
        enddo
        risquar = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
        rimag = SQRT( risquar )
c
        alamini = alamin(ityp)
        emax = emxslo
        rinmax = emax*(alamini+alaminn) / (alamini*alaminn)
c
c  Test overlap of atom i and short-range pot due to "nuc"
        if( risquar.gt.rinmax ) goto 2000
c
        limax = 0
        do  ishl=1,nshli
          if( lshel(ishl,ityp).gt.limax ) limax = lshel(ishl,ityp)
        enddo
        call AINITF( ai,gradai )
        if( rimag .lt. 1.d-12 )then
          lpowi = 0
          limxp1 = 1
          ylmi(1,1) = one/rt4pi
          call MKZERO( 6, grylmi )
c
          ylmi_t(1,1,1) = ylmi(1,1)
          call MKZERO( 6, grylmi_t )
c
        else
          lpowi = 1
          limxp1 = limax + lmxp1
          call AOFRF( ai,gradai, ri )
          call DYLM( ylmi,grylmi, limax+lmxnl, ri )
c
          ilm = 0
          do  mmi=1,limxp1
            do  lli=mmi,limxp1
              ilm = ilm + 1
              do  ipprl=1,ime(mmi)
                ylmi_t(ipprl,mmi,lli) = ylmi(ipprl,ilm)
                grylmi_t(1,ipprl,mmi,lli) = grylmi(ipprl,1,ilm)
                grylmi_t(2,ipprl,mmi,lli) = grylmi(ipprl,2,ilm)
                grylmi_t(3,ipprl,mmi,lli) = grylmi(ipprl,3,ilm)
              enddo
            enddo
          enddo
c
        endif
c
        i = istart
        do 1800 ishl=1,nshli
          nali = nala(ishl,ityp)
          lshli = lshel(ishl,ityp)
          li0 = lshli*lshli
          norbshi = 2*lshli + 1
          i0 = i
c
          almini = ala(1,ishl,ityp)
          almaxi = ala(nali,ishl,ityp)
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.lt.almini ) almini = ali
            if( ali.gt.almaxi ) almaxi = ali
          enddo
          alamini = almini
c
          emax = emxslo
          rinmax = emax*(alamini+alaminn) / (alamini*alaminn)
c
c  Test whether this radial function will contribute:
          if( risquar.gt.rinmax ) goto 1800
c
c    Atom j loop:
c
          jnext = 0
          do 1000 jatm=1,iatm
            jtyp = itypa(jatm)
            jstart = jnext
            jnext = jstart + norba(jtyp)
            nshlj = numshl(jtyp)
c
c  Following helps ensure only j.le.i done:
            if( jatm.eq.iatm ) nshlj = ishl
c
            alaminj = alamin(jtyp)
c  Criterion for including contribution of atom j:
            ex = risquar*alamini*alaminn / (alamini+alaminn)
            rijnmax = (emaxsl-ex)*(alamini+alaminj+alaminn) /
     $                          ( (alamini+alaminn)*alaminj )
            riinmax = (emaxii-ex)*(alamini+alaminj+alaminn) /
     $                          ( (alamini+alaminn)*alaminj )
c
            if( jatm.eq.iatm )then
              emaxij = emaxii
              rijnmx0 = riinmax
            else
              emaxij = emaxsl
              rijnmx0 = rijnmax
            endif
            rjnmax = emaxij*(alaminj+alaminn) / (alaminj*alaminn)
            rjnmx1 = rjnmax
            riovlp = risquar * ( alamini / (alamini+alaminn) )**2
c
            ljmax = 0
            do  jshl=1,nshlj
              if( lshel(jshl,jtyp).gt.ljmax ) ljmax = lshel(jshl,jtyp)
            enddo
c
c  Find nearest jatm image:
            rjsqnxt = 9.d9
            do 90 jlat=1,nlat2c
              do  jr=1,3
                rj(jr) = ratm(jr,jatm) + rlat(jr,jlat) - rnuc(jr)
              enddo
              rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
              if( rjsquar.lt.rjsqnxt )then
                rjsqnxt = rjsquar
                jlatnxt = jlat
              endif
   90       continue
c
c     Lattice vector loop:
c
c Need to preset initial jlat now.  I:I term, if it appears in star,
c appears in jlat=1.  If it is found, exit "star search", do the
c integrals, and return with jlat1=2 to create non-I:I part of star
            jlat1 = 1
c
            do 900 jlatx=1,nlat2c
              rjsquar = rjsqnxt
              rjsqnxt = 9.d9
c
c  First test that jatm overlaps nuc ...
              if( rjsquar.gt.rjnmx1 ) goto 1000
c
c   ... and test for 3-center overlap as well:
c       (I align the ri and rj vectors for this test)
              rjmag = SQRT( rjsquar )
              if( rjsquar.gt.riovlp )then
                rjnmx2 = rijnmx0 * ( (alamini+alaminn) /
     $           (alamini+alaminn-alamini*rimag/rjmag) )**2
                if( rjsquar.gt.rjnmx2 ) goto 1000
              endif
c
c  Align the jatm test vector along the iatm vector
              if( rimag .gt. 1.d-6 )then
                scale = rjmag/rimag
                do  jr=1,3
                  rjsav(jr) = scale*ri(jr)
                enddo
              else
                do  jr=1,3
                  rjsav(jr) = ratm(jr,jatm)+rlat(jr,jlatnxt) - rnuc(jr)
                enddo
              endif
c
              f1 = alamini / (alamini+alaminn)
              rijn0 = zero
              do  jr=1,3
                rijn0 = rijn0 + (f1*ri(jr)-rjsav(jr))**2
              enddo
c
              iatmi = .false.
              jeq = 0
              do 190 jlat=jlat1,nlat2c
c
                do  jr=1,3
                  rj(jr) = ratm(jr,jatm) + rlat(jr,jlat) - rnuc(jr)
                enddo
                rjsqtst = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
c  Compile orientation info for atoms such that rjsq=rjsquar
c
                if( rjsqtst-rjsquar .gt. 1.d-6
     $           .or. ( iatmi .and. rjsqtst.ge.rjsquar ) )then
c
c  This rjsq is bigger, check if it is next biggest ...
                  if( rjsqtst.lt.rjsqnxt )then
                    rjsqnxt = rjsqtst
                    jlatnxt = jlat
                  endif
c
                elseif( ABS(rjsqtst-rjsquar) .le. 1.d-6 )then
c
c  This rjsq is within current star, process it
c
                  if( jlat.eq.1 .and. jatm.eq.iatm )then
                    if( rijn0.gt.riinmax ) goto 190
                    iatmi = .true.
                  else
c                   Check tighter cutoff ...
                    if( rijn0.gt.rijnmax ) goto 190
                  endif
c
                  jeq = jeq + 1
                  if( jeq.gt.neqlatd ) goto 190
c
                  lateqj(jeq) = jlat
                  do  jr=1,3
                    rjstar(jr,jeq) = rj(jr)
                  enddo
c
                  call AINITF( aj(1,1,jeq),gradaj(1,1,1,jeq) )
                  if( rjsqtst .lt. 1.d-12 )then
                    ljmxp1 = 1
                    ylmj(1,1,jeq) = one/rt4pi
                    call MKZERO( 6, grylmj(1,1,1,jeq) )
c
                    ylmj_t(1,1,1,jeq) = ylmj(1,1,jeq)
                    call MKZERO( 6, grylmj_t(1,1,1,1,jeq) )
c
                  else
                    ljmxp1 = ljmax + lmxp1
                    call AOFRF( aj(1,1,jeq),gradaj(1,1,1,jeq), rj )
                    call DYLM( ylmj(1,1,jeq),grylmj(1,1,1,jeq),
     $                         ljmax+lmxnl, rj )
c
                    jlm = 0
                    do  mmj=1,ljmxp1
                      do  llj=mmj,ljmxp1
                        jlm = jlm + 1
                        do  jpprl=1,ime(mmj)
                          ylmj_t(jpprl,mmj,llj,jeq) =
     $                     ylmj(jpprl,jlm,jeq)
                          grylmj_t(1,jpprl,mmj,llj,jeq) =
     $                     grylmj(jpprl,1,jlm,jeq)
                          grylmj_t(2,jpprl,mmj,llj,jeq) =
     $                     grylmj(jpprl,2,jlm,jeq)
                          grylmj_t(3,jpprl,mmj,llj,jeq) =
     $                     grylmj(jpprl,3,jlm,jeq)
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
                jlat1 = 2
                if( rijn0.gt.riinmax ) goto 900
              else
                jlat1 = 1
                if( rijn0.gt.rijnmax ) goto 900
              endif
c
c  Check how many atoms in star we picked up
              neqjl = jeq
              if( neqjl.gt.mxeqjl ) mxeqjl = neqjl
              if( neqjl.gt.neqlatd )then
                write(IWR,*) '>>>>> neqjl,neqlatd=',neqjl,neqlatd
                call STOPXERR( 'neqlatd - exceeded in NLOCFRC' )
              endif
              if( neqjl.le.0 ) goto 900
c
c  Find largest l for atom j for which an alpha will contribute:
              jlmax = -1
              do 280 jshl=1,nshlj
                nalj = nala(jshl,jtyp)
                lshlj = lshel(jshl,jtyp)
c
                alminj = ala(1,jshl,jtyp)
                almaxj = ala(nalj,jshl,jtyp)
                do  jal=1,nalj
                  alj = ala(jal,jshl,jtyp)
                  if( alj.lt.alminj ) alminj = alj
                  if( alj.gt.almaxj ) almaxj = alj
                enddo
c
                emax = emaxf
                if( almaxi.lt.alfast .and. almaxj.lt.alfast )emax=emaxsl
                if( iatmi ) emax = emaxii
c
                ex = risquar*alamini*alaminn / (alamini+alaminn)
                rijnmx = (emax-ex)*(alamini+alminj+alaminn) /
     $                           ( (alamini+alaminn)*alminj )
                if( rijn0.le.rijnmx )then
                  if( lshlj.gt.jlmax ) jlmax = lshlj
                endif
  280         continue
c
              if( jlmax.lt.0 ) goto 900
c
              rjmag = SQRT( rjsquar )
c
              lisq = limxp1*limxp1
              ljsq = ljmxp1*ljmxp1
c
              maxli = (lshli+1)*(lshli+1)
              maxlj = (jlmax+1)*(jlmax+1)
c
c carve up space in angsav() for arrays in angsetf:
              iangs = 1
c angsav() needs:
              nangs = limxp1*ljmxp1*7*lmxp1*9*9
              iangt = iangs + nangs
c angtmp() needs:
              nangt = lisq*ljsq
              ignti = iangt + nangt
c gntci() (complex, use 4* real space for safety) needs:
              ngnti = 4*(lisq*maxli*lmxp1*lmxp1)
              igntj = ignti + ngnti
c cgntcj() (also complex) needs:
              ngntj = 4*(ljsq*maxlj*lmxp1*lmxp1)
c
              if( neqjl.eq.1 )then
c
c  If I have only one atom in a star, then compute the angular
c  stuff up here only once, rather than repeat for every shell
c  If more than one atom per star, I either have to duplicate
c  angsav for each element in the star, or, as done here, compute
c  it below for each star. 15Feb00-PAS/2.43
c
                jeq = 1
                call NLANGF(
     $           maxli,maxlj,limxp1,ljmxp1,lmxp1,lisq,ljsq,
     $           ylmi_t,ylmj_t(1,1,1,jeq),
     $           grylmi_t,grylmj_t(1,1,1,1,jeq),gntcomb,
     $           angsav,angsav(iangt),angsav(ignti),angsav(igntj) )
              endif
c
c    Loop over shells on j-atm
c
              j = jstart
              do 800 jshl=1,nshlj
                nalj = nala(jshl,jtyp)
                lshlj = lshel(jshl,jtyp)
                lj0 = lshlj*lshlj
                norbshj = 2*lshlj + 1
                j0 = j
c
                norbshj0 = norbshj
                joff = 0
                if( j0.eq.i0 )then
                  norbshj0 = 0
                  joff = 1
                endif
c
                alminj = ala(1,jshl,jtyp)
                almaxj = ala(nalj,jshl,jtyp)
                do  jal=1,nalj
                  alj = ala(jal,jshl,jtyp)
                  if( alj.lt.alminj ) alminj = alj
                  if( alj.gt.almaxj ) almaxj = alj
                enddo
c
                emax = emaxf
                if( almaxi.lt.alfast .and. almaxj.lt.alfast )emax=emaxsl
                if( iatmi ) emax = emaxii
c
                lpowmx = lpowi
                lpowj = 0
                if( rjmag .ge. 1.d-12 )then
                  lpowj = 1
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
c    Alpha loops:
c
                nx = 0
                do 320 ial=1,nali
                  ali = ala(ial,ishl,ityp)
                  cali = cala(ial,ishl,ityp)
                  twoali = ali + ali
c
c  Test whether this alpha contributes:
                  rinmx = emax*(ali+alaminn) / (ali*alaminn)
                  if( risquar.gt.rinmx ) goto 320
c
                  fi1 = ali / (ali+alaminn)
                  rijn = zero
                  do  jr=1,3
                    rijn = rijn + (fi1*ri(jr)-rjsav(jr))**2
                  enddo
c
                  do 310 jal=1,nalj
                    alj = ala(jal,jshl,jtyp)
                    calj = cala(jal,jshl,jtyp)
                    twoalj = alj + alj
c
c  Test whether this pair of alphas contributes:
                    ex = risquar*ali*alaminn / (ali+alaminn)
                    rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                               ( (ali+alaminn)*alj )
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
c  Evaluate radial integrals
c
                    call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $               exparg,bargs,alsum,radmsh(1,ktyp), wkrad )
                    nr0 = nrmn - 1
                    nrlen = nrmx - nrmn + 1
c
                    call NLRAD( vnlint,  alsum,bargi,bargj,bargs,
     $               exparg, limxpd,ljmxpd,lmxp1,lstmx,lpowmx,
     $               nrlen,nrd,nr0, radmsh(1,ktyp),vpsrad(1,1,ktyp),
     $               wkrad(1,1),wkrad(1,8),wkrad(1,15),
     $               wkrad(1,22),wkrad(1,23),wkrad(1,24),
     $               wkrad(1,25),wkrad(1,26),wkrad(1,27) )
c
c  Sum into running totals:
c
                    cnln = cali*calj
                    cnli = cnln*twoali
                    cnlj = cnln*twoalj
                    do  lij=1,lijmx
                      vnln(lij) = vnln(lij) + cnln*vnlint(lij)
                      vnli(lij) = vnli(lij) + cnli*vnlint(lij)
                      vnlj(lij) = vnlj(lij) + cnlj*vnlint(lij)
                    enddo
c
c    Close alpha loops:
c
  310             continue
  320           continue
c
                if( nx.eq.0 ) goto 800
c
                do 490 jeq=1,neqjl
c
                  jlat = lateqj(jeq)
                  do  jr=1,3
                    rj(jr) = rjstar(jr,jeq)
                  enddo
c
                  if( neqjl.gt.1 )then
                    call NLANGF(
     $               maxli,maxlj,limxp1,ljmxp1,lmxp1,lisq,ljsq,
     $               ylmi_t,ylmj_t(1,1,1,jeq),
     $               grylmi_t,grylmj_t(1,1,1,1,jeq),gntcomb,
     $               angsav,angsav(iangt),angsav(ignti),angsav(igntj) )
                  endif
c
c    Angular loops:
c
                  norbshj = norbshj0
                  i = i0
                  do 420 li=li0+1,li0+norbshi
                    i = i + 1
                    norbshj = norbshj + joff
                    j = j0
                    do 410 lj=lj0+1,lj0+norbshj
                      j = j + 1
                      denlilj = vnorm(li)*vnorm(lj)
c
c  Do sum over density matrix:
                      if( i.eq.j )then
                        dmsum = zero
                        do  k=1,nk
                          dmsum = dmsum + coskr(k,jlat)*dmat(i,i,k)
                        enddo
                      else
                        dmsum = zero
                        do  k=1,nk
                          dmsum = dmsum + coskr(k,jlat)*dmat(i,j,k)
     $                                  - sinkr(k,jlat)*dmat(j,i,k)
                        enddo
                        dmsum = dmsum + dmsum
                      endif
                      dmsum = dmsum/denlilj
c
                      call NLINTF( enpri,enprj, engri,engrj,
     $                 vnln,vnli,vnlj, rimag,rjmag, li,lj,
     $                 limxp1,ljmxp1,lmxp1,lpowmx,lpowi,lpowj,lstmx,
     $                 ai(1,li),    gradai(1,1,li),
     $                 aj(1,lj,jeq),gradaj(1,1,lj,jeq), angsav )
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
c    Close angular loops:
c
  410               continue
  420             continue
c
c    Close star of jatm loop:
c
  490           continue
c
c    Close atom j shell loop:
c
  800         j = j0 + norbshj
c
c    Close atom j lattice loop:
c
  900       continue
c
c    Close atom j loop:
c
 1000     continue
c
c    Close atom i shell loop:
c
 1800   i = i0 + norbshi
c
c    Close atom i loop:
c
 2000 continue
c
c    Close nucleus lattice loop:
c
 2900 continue
c
c    Close nucleus loop:
c
 3000 continue
ctp:start
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
ctp:end
c
c    That's all Folks!
c
      RETURN
      END
