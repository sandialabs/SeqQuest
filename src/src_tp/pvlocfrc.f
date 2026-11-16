c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PVLOCFRC
c
c
      subroutine PVLOCFRC( nafrc,iforce, ndim ,frcloc,strloc,
     $ dolocxc, norb,nk, natm,ntyp,nshld,nald, nlat,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ lmxnlp1, ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vpsrad, vatrad,vesrad, wkrad,
     $ naldsq, vlocez,vloces,vlocxc, twoali,twoalj,twoals, p,
     $ wksm,
     $ dmat, frctp )
c---------------------------------------------------------------
c Purpose: calculate forces&stress from spherical local potentials
c
c Written: Peter A. Schultz, based on P.J. Feibelman's original
c
c Revision history:
c  11May07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  15Jun05-PAS/2.59: File units
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines.
c   6May03-DMC: Add frclocl to parameter list, to provide a local
c               copy of the frcloc array, for parallelization
c   7Mar02-PAS/2.52: ang constants extracted
c  23Jul01-PAS/2.49: compact dmat
c  22May01-PAS/2.47: condense radial scratch space
c  14Mar00-PAS/2.43: rename "vl" routines
c  21Dec99-PAS/2.42: cosmetic clean-up
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  17Dec98-PAS/2.30: clean force/stress/arrays
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  30Sep98-PAS/2.27: rework angular integral section (including
c                    interesting legacy bug), rename
c  27Aug98-PAS/2.25: purge (N**7) equivs and condense arithmetic;
c                 use "engausf' only for reference occupied terms
c  11Aug98-PAS/2.24: match 3-center routines again
c  23Jun98-PAS/2.23: output controls + occ/slow cutoffs done right
c                    split basis
c  28Apr97-PAS/2.17: naming conventions
c   3Jan96-PAS/2.13: clean cutoffs/opt
c   6Jan95-PAS/2.11: purge explicit complex
c  17Oct94-PAS/2.10: cleanup and optimization
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c output arrays:
      DIMENSION  frcloc(3,natm), strloc(3,3)
c input arrays:
      LOGICAL    dolocxc
      DIMENSION  iforce(natm)
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $  nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  znuc(ntyp), occsh(nshld,ntyp),occij(nshld,nshld,ntyp)
      DIMENSION  lmxnlp1(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp),
     $           vpsrad(nrd,4,ntyp),vatrad(nrd,ntyp),vesrad(nrd,ntyp)
      DIMENSION  dmat(norb,norb,nk)
c scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vlocez(6,naldsq),vloces(6,naldsq),vlocxc(6,naldsq)
      DIMENSION  twoali(naldsq),twoalj(*),twoals(*)
      DIMENSION  p(3,naldsq)
      DIMENSION  wksm(naldsq,*)
ctp: parallel scratch:
      DIMENSION  frctp(3,natm)
      DIMENSION  strtp1(11),strtp2(11)
c local declarations:
      LOGICAL    iatmi
      LOGICAL    keepsh, do_es,do_ps,do_xc
      DIMENSION  vesint(6),vpsint(6),vlocint(6)
      DIMENSION  fez(3),fesi(3),fesj(3),fxci(3),fxcj(3)
c
      DIMENSION  rnuc(3),ri(3),rj(3), rij(3)
      DIMENSION  ci(10,9),cj(10,9),gradci(3*10,9),gradcj(3*10,9)
      DIMENSION  vzero(270)
c
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxf = cutgwfn
      emaxsl = cutslow
      emaxii = cutii
      emxslo = MAX( emaxsl, emaxii )
c
      fourpi = pi + pi + pi + pi
c
ctp:
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c  Initialize stress and force
c
      strszz = zero
      strsnn = zero
      call MKZERO( 3*natm, frcloc )
      call MKZERO( 9, strloc )
c
      call MKZERO( 270, vzero )
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
      alaminn = two*alamin(ktyp)
      almin = 1.d9
      do  ishl=1,nocc(ktyp)
        lshl = locc(ishl,ktyp)
        occshl = occsh(ishl,ktyp)
        if( occshl.ne.zero )then
          do  ial=1,nalsh(ishl,ktyp)
            if( alsh(ial,ishl,ktyp).lt.almin ) almin=alsh(ial,ishl,ktyp)
          enddo
        endif
      enddo
      alaminn = two*almin
c
      nmat = nk*norb*norb
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
c
      alamini = alamin(ityp)
      emax = emxslo
      rinmax = emax*(alamini+alaminn) / (alamini*alaminn)
c
c  Test overlap of atom i and short-range pot due to "nuc"
      if( risquar.gt.rinmax ) goto 2000
c
      if( klat.ne.1 .or. iatm.ne.nuc )then
c  Add nuc-iatm core repulsion to forces.  Note: skipped IF
c  no nuc:iatm overlap (electronic and nuclear contributions cancel).
        zi = znuc(ityp)
        do 28 ishl=1,nshli
c
          norbshi = 2*lshel(ishl,ityp) + 1
          almini = ala(1,ishl,ityp)
          do  ial=1,nala(ishl,ityp)
            ali = ala(ial,ishl,ityp)
            if( ali.lt.almini ) almini = ali
          enddo
c
          do 26 jshl=1,ishl
c
c Note that ishl and jshl must have same angular momentum
c and that z contribution is in jshl.le.ishl part
            occshl = norbshi*occij(jshl,ishl,ityp)
            if( occshl.eq.zero ) goto 26
c
            alminj = ala(1,jshl,ityp)
            do  jal=1,nala(jshl,ityp)
              alj = ala(jal,jshl,ityp)
              if( alj.lt.alminj ) alminj = alj
            enddo
            alsum = almini + alminj
            emax = emaxii
            rismax = emax*(alsum+alaminn) / (alsum*alaminn)
            if( risquar.gt.rismax )  zi = zi - occshl
   26     continue
c
   28   continue
        zztot = znuc(ktyp)*zi
c
        engyadd = zztot / SQRT( risquar )
        facn = engyadd/risquar
        frcloc(1,nuc) = frcloc(1,nuc) - two*facn*ri(1)
        frcloc(2,nuc) = frcloc(2,nuc) - two*facn*ri(2)
        frcloc(3,nuc) = frcloc(3,nuc) - two*facn*ri(3)
        risqd = zero
        do  jd=1,ndim
          risqd = risqd + ri(jd)*ri(jd)
          do  id=1,ndim
            strloc(id,jd) = strloc(id,jd) - facn*ri(id)*ri(jd)
          enddo
        enddo
        strszz = strszz - facn*risqd
      endif
c
      call CINITF( ci,gradci )
      call COFRF( ci,gradci,ri )
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
c  criterion for including contribution of atom j:
        ex = risquar*alamini*alaminn / (alamini+alaminn)
        rijnmax = (emaxsl-ex)*(alamini+alaminj+alaminn) /
     $                      ( (alamini+alaminn)*alaminj )
        riinmax = (emaxii-ex)*(alamini+alaminj+alaminn) /
     $                      ( (alamini+alaminn)*alaminj )
c
c    Lattice vector loop:
c
        do 900 jlat=1,nlat2c
c
          do  jr=1,3
            rij(jr) = ratm(jr,iatm) - ratm(jr,jatm) - rlat(jr,jlat)
            rj(jr) = ratm(jr,jatm) + rlat(jr,jlat) - rnuc(jr)
          enddo
          rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
          f1 = alamini / (alamini+alaminn)
          rijn0 = zero
          do  jr=1,3
            rijn0 = rijn0 + (f1*ri(jr)-rj(jr))**2
          enddo
c
c  Test for significant overlap of atom j with i and nuc
          if( jatm.eq.iatm .and. jlat.eq.1 )then
            iatmi = .true.
            if( rijn0.gt.riinmax ) goto 900
          else
            iatmi = .false.
            if( rijn0.gt.rijnmax ) goto 900
          endif
c
          call CINITF( cj,gradcj )
          call COFRF( cj,gradcj,rj )
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
            if( almaxi.lt.alfast .and. almaxj.lt.alfast ) emax = emaxsl
            if( iatmi ) emax = emaxii
c
            ex = risquar*alamini*alaminn / (alamini+alaminn)
            rijnmx = (emax-ex)*(alamini+alminj+alaminn) /
     $                       ( (alamini+alaminn)*alminj )
            if( rijn0.gt.rijnmx ) goto 800
c
            keepsh = .false.
            occ = zero
            if( iatmi )then
              if( ishl.gt.jshl )then
                occ = two*occij(ishl,jshl,ityp)
              else
                occ = occij(jshl,ishl,ityp)
              endif
              if( occ.ne.zero ) keepsh = .true.
            endif
c
            lmx = lshli + lshlj + 2
c
c    Alpha loops:
c
            ix = 0
            ixx = 0
            do 320 ial=1,nali
              ali = ala(ial,ishl,ityp)
c
c  Test whether this alpha contributes:
              rinmx = emax*(ali+alaminn) / (ali*alaminn)
              if( risquar.gt.rinmx .and. .not. keepsh ) goto 320
c
              fi1 = ali / (ali+alaminn)
              rijn = zero
              do  jr=1,3
                rijn = rijn + (fi1*ri(jr)-rj(jr))**2
              enddo
c
              do 310 jal=1,nalj
                alj = ala(jal,jshl,jtyp)
c
c  Test whether this pair of alphas contributes:
                ex = risquar*ali*alaminn / (ali+alaminn)
                rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                           ( (ali+alaminn)*alj )
                if( rijn.gt.rijnmx .and. .not. keepsh ) goto 310
c
                ix = ix + 1
c
                alsum = ali + alj
                exparg = ali*risquar + alj*rjsquar
                twoals(ix) = alsum + alsum
                tooali = ali + ali
                tooalj = alj + alj
                twoali(ix) = tooali
                twoalj(ix) = tooalj
                p(1,ix) = tooali*ri(1) + tooalj*rj(1)
                p(2,ix) = tooali*ri(2) + tooalj*rj(2)
                p(3,ix) = tooali*ri(3) + tooalj*rj(3)
                pmag = SQRT( p(1,ix)**2 + p(2,ix)**2 + p(3,ix)**2 )
c
                enfac = fourpi*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
c
                if( keepsh )then
c
c  Compute radial integrals, put matrix elements of -2/r in vlocez:
c
c expfac is for underflow control
                  expfac = zero
                  if( exparg.lt.argmx ) expfac =  EXP( -exparg )
                  call VLRADZ( vlocez(1,ix),lmx,alsum,pmag,-one,
     $                         exparg,expfac )
c
c      compute contribution to force integral from hellmann-feynmann
c      term involving the spherical part of the density matrix
c
                  do  ladd=1,lmx
                    vlocez(ladd,ix) = vlocez(ladd,ix)*znuc(ktyp)*enfac
                  enddo
                else
c  blank it, since it is too hard to stop its use later ...
                  do  ladd=1,lmx
                    vlocez(ladd,ix) = zero
                  enddo
                endif
c
                do_es = .false.
                do_ps = .false.
                do_xc = .false.
                do  ladd=1,lmx
                  vesint(ladd) = zero
                  vpsint(ladd) = zero
                  vlocint(ladd) = zero
                enddo
                if( rijn.gt.rijnmx ) goto 307
c
c install underflow cutoffs:
                call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $           exparg,pmag,alsum, radmsh(1,ktyp), wkrad )
                if( nrmx.le.0 ) goto 307
c
                ixx = ixx + 1
c Flag to compute matrix elements of ...
c    ... short range spherical electrostatic potential
                do_es = .true.
c    ... hellman-feynmann contribution from local pspot
                do_ps = .false.
                if( lmxnlp1(ktyp).ge.0 ) do_ps = .true.
c    ... all of the short range potential (to get xc piece)
                do_xc = dolocxc
c
                nrlen = nrmx - nrmn + 1
                lp = lmxnlp1(ktyp) + 1
                call VLRADF( vesint,vpsint,vlocint, do_es,do_ps,do_xc,
     $           lmx,alsum,pmag,exparg,    nrlen,
     $           radmsh(nrmn,ktyp),radwt(nrmn,ktyp), vesrad(nrmn,ktyp),
     $           vpsrad(nrmn,lp,ktyp),vatrad(nrmn,ktyp),
     $           wkrad, wkrad(1, 8),wkrad(1, 9),wkrad(1,10),
     $                  wkrad(1,11),wkrad(1,12),wkrad(1,13) )
c
  307           continue
c
                if( dolocxc )then
c                Atom reference xc potential forces computed here
                 do  ladd=1,lmx
                  vlocez(ladd,ix) = vlocez(ladd,ix) - enfac*vesint(ladd)
                  vloces(ladd,ix) = enfac*(vesint(ladd)+vpsint(ladd))
                  vlocxc(ladd,ix) = enfac*vlocint(ladd)
     $                            - vloces(ladd,ix)
                 enddo
                else
c                Atom reference xc potential forces computed on grid
                 do  ladd=1,lmx
                  vlocez(ladd,ix) = vlocez(ladd,ix) - enfac*vesint(ladd)
                  vloces(ladd,ix) = enfac*(vesint(ladd)+vpsint(ladd))
                  vlocxc(ladd,ix) = zero
                 enddo
                endif
c
  310         continue
  320       continue
c
            nx = ix
            if( nx.eq.0 ) goto 800
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
                if( keepsh .and. li.eq.lj )then
c
                  call VLOCF0( fez(1),fez(2),fez(3),
     $             li,lj, ci(1,li),cj(1,lj), gradci(1,li),gradcj(1,lj),
     $             nx, vlocez, twoals, p,
     $             wksm(1,1), wksm(1,2),wksm(1,5),
     $             wksm(1,8),wksm(1,9),wksm(1,10),
     $             wksm(1,11),wksm(1,12),wksm(1,13),wksm(1,14) )
c
                  occfac = occ/denlilj
                  fez(1) = occfac*fez(1)
                  fez(2) = occfac*fez(2)
                  fez(3) = occfac*fez(3)
c
                  frcloc(1,nuc) = frcloc(1,nuc) - fez(1)
                  frcloc(2,nuc) = frcloc(2,nuc) - fez(2)
                  frcloc(3,nuc) = frcloc(3,nuc) - fez(3)
c
                  do  jd=1,ndim
                    strsnn = strsnn - fez(jd)*ri(jd)/two
                    do  id=1,ndim
                      strloc(id,jd) = strloc(id,jd) - fez(id)*ri(jd)/two
                    enddo
                  enddo
c
                endif
c
                if( ixx.eq.0 ) goto 410
c
c  Compute density matrix k-sum:
c
                if( i.eq.j )then
                  dmsum = zero
                  do  k=1,nk
                    dmsum = dmsum + coskr(k,jlat)*dmat(i,i,k)
                  enddo
                else
                  dmsum = zero
                  do  k=1,nk
                    dmsum = dmsum + coskr(k,jlat)*dmat(i,j,k)
     $                            - sinkr(k,jlat)*dmat(j,i,k)
                  enddo
                  dmsum = dmsum + dmsum
                endif
                dmsum = dmsum/denlilj
c
                fesi0 = zero
                fxci0 = zero
                fesj0 = zero
                fxcj0 = zero
                do  ir=1,3
                  fesi(ir) = zero
                  fxci(ir) = zero
                  fesj(ir) = zero
                  fxcj(ir) = zero
                enddo
c
                if( dolocxc )then
c
                  call VLOCF2(
     $             fesi0,fesi, fesj0,fesj, fxci0,fxci, fxcj0,fxcj,
     $             li,lj, ci(1,li),cj(1,lj),gradci(1,li),gradcj(1,lj),
     $             nx, vloces,vlocxc, twoali,twoalj, p,
     $             wksm(1, 1), wksm(1, 2),wksm(1, 3),wksm(1, 6),
     $                        wksm(1, 9),wksm(1,10),wksm(1,13),
     $             wksm(1,16),wksm(1,17),wksm(1,18),
     $             wksm(1,19),wksm(1,20),wksm(1,21),wksm(1,22) )
c
                else
c
                  call VLOCF1( fesi0,fesi, fesj0,fesj,
     $             li,lj, ci(1,li),cj(1,lj),gradci(1,li),gradcj(1,lj),
     $             nx, vloces, twoali,twoalj, p,
     $             wksm(1, 1), wksm(1, 2),wksm(1, 3),wksm(1, 6),
     $                        wksm(1, 9),wksm(1,10),wksm(1,13),
     $             wksm(1,16),wksm(1,17),wksm(1,18),
     $             wksm(1,19),wksm(1,20),wksm(1,21),wksm(1,22) )
c
                endif
c
                do 406 ir=1,3
                  fesi(ir) = dmsum*( fesi(ir) + fesi0*ri(ir) )
                  fxci(ir) = dmsum*( fxci(ir) + fxci0*ri(ir) )
                  fesj(ir) = dmsum*( fesj(ir) + fesj0*rj(ir) )
                  fxcj(ir) = dmsum*( fxcj(ir) + fxcj0*rj(ir) )
  406           continue
c
                frcloc(1,iatm) = frcloc(1,iatm) + fesi(1) + fxci(1)
                frcloc(2,iatm) = frcloc(2,iatm) + fesi(2) + fxci(2)
                frcloc(3,iatm) = frcloc(3,iatm) + fesi(3) + fxci(3)
c
                frcloc(1,jatm) = frcloc(1,jatm) + fesj(1) + fxcj(1)
                frcloc(2,jatm) = frcloc(2,jatm) + fesj(2) + fxcj(2)
                frcloc(3,jatm) = frcloc(3,jatm) + fesj(3) + fxcj(3)
c
                frcloc(1,nuc)  = frcloc(1,nuc)  - fesi(1) - fesj(1)
                frcloc(2,nuc)  = frcloc(2,nuc)  - fesi(2) - fesj(2)
                frcloc(3,nuc)  = frcloc(3,nuc)  - fesi(3) - fesj(3)
c
c The factors for the xc terms in the stress calculations below
c are suspect.  Original code was bogus, and I'm guessing about
c the correct coordinate factors and coefficients.  Fortunately
c this code is inactive with the xc flags set as they are.
c At some point I'll have to clear this up. 30Sep98-PAS
c
                do  jd=1,ndim
                  do  id=1,ndim
                    strloc(id,jd) = strloc(id,jd)
     $               - fesi(id)*ri(jd) - fesj(id)*rj(jd)
     $               - fxci(id)*rij(jd)
     $               + fxcj(id)*rij(jd)
                  enddo
                enddo
c
c    Close angular loops:
c
  410         continue
  420       continue
c
c    Close atom j shell loop:
c
  800     j = j0 + norbshj
c
c    Close atom j lattice loop:
c
  900   continue
c
c    Close atom j loop:
c
 1000 continue
c
c    Close atom i shell loop:
c
 1800 i = i0 + norbshi
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
c       Each local processor has forces over its share of work.
c       Need to sum these to master (node 0).
c       For now, also broadcast total to all processors.
c
        nreduc = 3*natm
        call MPREDUC8( master, nreduc, frcloc, frctp , icomm )
C        call MPI_Reduce( frcloc, frctp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocfrc - Reduce frcloc' )
        call DCOPY( 3*natm, frctp,1, frcloc,1 )
c
        call MPBCAST8( master, nreduc, frcloc, icomm )
C        call MPI_Bcast( frcloc, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocfrc - Bcast frcloc' )
c
c       Want to add both strszz and strsnn to what is accumulated.
c       Will put it in as part of the strtp array for comms.
c
        call DCOPY( 9, strloc,1, strtp1,1 )
        strtp1(10) = strszz
        strtp1(11) = strsnn
        nreduc = 11
c
        call MPREDUC8( master, nreduc, strtp1, strtp2, icomm )
C        call MPI_Reduce( strtp1, strtp2, nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocfrc - Reduce strloc' )
c
        call MPBCAST8( master, nreduc, strtp2, icomm )
C        call MPI_Bcast( strtp2, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocfrc - Bcast strlll' )
c
        call DCOPY( 9, strtp2,1, strloc,1 )
        strszz = strtp2(10)
        strsnn = strtp2(11)
c
c       End tp
      endif
ctp:end
c
c  Print some results:
c
      if( lstout.gt.2 )then
        call FLGETIWR( IWR )
        write(IWR,'(4x,a,f20.8)') 'stress/zz=',strszz
        write(IWR,'(4x,a,f20.8)') 'stress/nn=',strsnn
      endif
c
c    That's all Folks!
c
      RETURN
      END
