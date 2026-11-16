c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLOCMAT
c
c
      subroutine NLOCMAT( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl, ratm, rlat,coskr,sinkr, gntcomb,
     $ nrd,nrad, radmsh,vpsrad, wkrad,
     $ vnlr,vnli,
     $ vnlij,  angsav, vnltp )
c---------------------------------------------------------------
c Purpose: calculate matrix elements due to non-local part of
c          the electron:effective-nuclear-potential interaction
c
c Revision history:
c  17Sep08-PAS/2.62: stripe-distributed-parallel
c  20Jun07-PAS/2.60: merge tp and serial (2oct07: small cosmetic)
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  16Apr03-DMC: Add vnltp to parameter list, to provide a local copy
c               of the vnlij array, for parallelism.
c   2Oct02-PAS/2.54: optimize NLINT, and clean out consequences
c   7Mar02-PAS/2.52: ang constants extracted
c  21Jun01-PAS/2.48: replace STOPs
c  22May01-PAS/2.47: condense radial scratch space
c   6Mar00-PAS/2.43: rename "nl" routines, "ennloc" to "nlocmat"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Sep98-PAS/2.27: cosmetic - naming convention changed
c   5Aug98-PAS/2.24: fix cutoff bug, full proper cutoff scaling
c  21Jul98-PAS/2.23: clean up cutoffs, extend I:N terms
c  22May98-PAS/2.22: bugfix
c  28Apr97-PAS/2.17: naming conventions
c  18Jan94-PAS/2.04: fix nlat, j-indexing bugs; rationalize
c         arguments; full cleanup, match sij,tij,enij,ennloc
c  21Jun89-PAS: extract ANGSET and some of LLOOPS out of under
c               orbital loops, pack angsav
c    May89-PAS: extract RADNLOC out of under angular loops
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c  Allows for non-locality up to l=lmxnlp1(ktyp)-1
c
c  The pseudopotential for atom type ktyp is assumed to fall to
c  zero as exp(-r**2*almnnl(ktyp)).  almnnl(ktyp) here is
c  taken from the fits of Bachelet, et al.
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
c Output array:
      DIMENSION  vnlij(nk,*)
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),vpsrad(nrd,4,ntyp)
c
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vnlr(nk),vnli(nk)
      DIMENSION  angsav(*)
c     Angsav scratch space packed below needs space 
c       5x5x5x3x9x9 + 4*(5x9x3x3) + 4x(5x5x9x3x3)=30375+1620+8100=40095
c Task-parallel scratch:
      DIMENSION  vnltp(*)
c
c Local declarations:
      LOGICAL    Lmatreduce
      LOGICAL    iatmi, diagshl
      LOGICAL    izero(9,9)
      DIMENSION  vnlint(5*5*15),vnltmp(5*5*15)
      DIMENSION  rnuc(3),ri(3),rj(3),risav(3)
c neqlatd in main routine should match this too.
      PARAMETER  ( neqlatd = 72 )
      DIMENSION  aj(9,9),ai(9,9,neqlatd), ilateq(neqlatd)
      DIMENSION  ylmj(2,25),ylmi(2,25,neqlatd)
c
      DIMENSION  vnorm(9)
c lstout sets level of diagnostic output, least:0-5:most
      DATA  lstout / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
      Lmatreduce = .false.
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
      nmat = nk*norb*norb
c
      mxeqil = 0
      if( lstout.gt.0 ) write(IWR,*) 'NLOCMAT: v-nonlocal diagnostics'
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
      nmatstr = ncplx*nk*matstr
c
c  Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
      nmatpr = mymatpr*ncplx*nk
c
      if( nprocs .gt. 0 ) call MKZERO( nmatpr, vnlij )
c
c    Nucleus loop:
c
      do 3000 nuc=1,natm
c
c     Each processor will do every nprocs element in this loop.
      if( Lmatreduce .and. ( MOD( nuc, nprocs ) .ne. iproc ) ) goto 3000
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
        call AINIT( aj )
        if( rjmag .lt. 1.d-12 )then
          ljmxp1 = 1
          ylmj(1,1) = one / rt4pi
c
        else
          ljmxp1 = ljmax + lmxp1
          call AOFR( aj, rj )
          call YLMREAL( ylmj, ljmax+lmxnl, rj )
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
c
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
            do 900 ilatx=1,nlat2c
              risquar = risqnxt
              risqnxt = 9.d9
c
c             First test that iatm overlaps nuc ...
              if( risquar.gt.rinmx1 ) goto 1000
c
c              ... and test for 3-center overlap as well:
c               (I align the ri and rj vectors for this test)
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
              do 190 ilat=1,nlat2c
c
                do  ir=1,3
                  ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
                enddo
                risqtst = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
c               Compile orientation info of atoms such that risq=risquar
c
                if( risqtst-risquar .gt.  1.d-6 )then
c
c                 This risq is bigger, check if it is next biggest ...
c
                  if( risqtst.lt.risqnxt )then
c                    ... it indeed is, and record it:
                    risqnxt = risqtst
                    ilatnxt = ilat
                  endif
c
                elseif( ABS( risqtst-risquar ) .le. 1.d-6 )then
c
c                 This risq is within current star, process it
c
                  if( ilat.eq.1 .and. iatm.eq.jatm )then
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
c
                  call AINIT( ai(1,1,ieqlat) )
                  if( risqtst .lt. 1.d-12 )then
                    limxp1 = 1
                    ylmi(1,1,ieqlat) = one / rt4pi
c
                  else
                    limxp1 = limax + lmxp1
                    call AOFR( ai(1,1,ieqlat), ri )
                    call YLMREAL( ylmi(1,1,ieqlat), limax+lmxnl, ri )
c
                  endif
c
                endif
c
  190         continue
c
              if( iatmi )then
                if( rijn0.gt.riinmax ) goto 900
              else
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
              maxlj = (lshlj+1)*(lshlj+1)
              maxli = (ilmax+1)*(ilmax+1)
c
c             Carve up space in angsav() for arrays in angset ...
              iangs = 1
c              ... angsav():
              nangs = ljmxp1*limxp1*limxp1*lmxp1*9*9
              igntj = iangs + nangs
c              ... gntci() (complex, allow 4* real space for safety):
              ngntj = 4*(ljmxp1*maxlj*lmxp1*lmxp1)
              ignti = igntj + ngntj
c              ... cgntcj() (also complex):
              ngnti = 4*(limxp1*limxp1*maxli*lmxp1*lmxp1)
c
c             Compute sums of three-j symbols:
              call NLANG( maxlj,maxli,ljmxp1,limxp1,lmxp1,rjmag,rimag,
     $         ylmj,gntcomb,angsav,angsav(igntj),angsav(ignti) )
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
                lpowmx = lshli + lshlj + 1
                lstmx = lpowmx*lmxp1
                lijmx = limxp1*ljmxp1*lstmx
                do  lij=1,lijmx
                  vnlint(lij) = zero
                enddo
c
c               Alpha loops:
c
                nx = 0
                do 320 jal=1,nalj
                  alj = ala(jal,jshl,jtyp)
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
c
c                   Test whether this pair of alphas contributes:
                    ex = rjsquar*alj*alaminn / (alj+alaminn)
                    rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                               ( (alj+alaminn)*ali )
                    if( rijn.gt.rijnmx ) goto 310
c
                    alsum = ali + alj
                    exparg = ali*risquar + alj*rjsquar
                    bargi = two*ali*rimag
                    bargj = two*alj*rjmag
                    bargs = bargi + bargj
c
c                   Evaluate radial integrals
c
                    call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $               exparg,bargs,alsum,radmsh(1,ktyp), wkrad )
                    nr0 = nrmn - 1
                    nrlen = nrmx - nrmn + 1
c
                    call NLRAD( vnltmp, alsum,bargj,bargi,bargs,
     $               exparg, ljmxp1,limxp1,lmxp1,lstmx,lpowmx,
     $               nrlen,nrd,nr0, radmsh(1,ktyp),vpsrad(1,1,ktyp),
     $               wkrad(1,1), wkrad(1,8), wkrad(1,15),
     $               wkrad(1,22),wkrad(1,23),wkrad(1,24),
     $               wkrad(1,25),wkrad(1,26),wkrad(1,27) )
c
                    fac = cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
                    do  lij=1,lijmx
                      vnlint(lij) = vnlint(lij) + fac*vnltmp(lij)
                    enddo
                    nx = nx + 1
c
c                   Close alpha loops:
  310             continue
  320           continue
c
                if( lstout.gt.3 )then
                  write(IWR,*) '   nuc,iatm,jatm,klat,ilatx,neqilat,nx'
                  write(IWR,9350)  nuc,iatm,jatm,klat,ilatx,neqilat,nx
 9350             format(1x,10i4)
                endif
c
                if( nx.eq.0 ) goto 799
c
c               Angular loops:
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
c                   Evaluate contribution to vnlij(i,j)
c
                    if( neqilat.eq.1 )then
                      call NLINT( vnlr,vnli, vnlint,
     $                 ljmxp1,limxp1,lmxp1,lstmx,lpowmx, rjmag,rimag,
     $                 lj,li, aj(1,lj),ai(1,li,1), ylmi,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 ilateq(1) )
                    else
c
                      do 405 lmpj=1,lj
                        if( aj(lmpj,lj).eq.zero )then
                          do  lmpi=1,li
                            izero(lmpj,lmpi) = .true.
                          enddo
                        else
                          do 404 lmpi=1,li
                            do  lateqi=1,neqilat
                              if( ai(lmpi,li,lateqi).ne.zero ) goto 403
                            enddo
                            izero(lmpj,lmpi) = .true.
                            goto 404
  403                       continue
                            izero(lmpj,lmpi) = .false.
  404                     continue
                        endif
  405                 continue
c
                      call NLINTEQ( vnlr,vnli, vnlint,
     $                 ljmxp1,limxp1,lmxp1,lstmx,lpowmx, rjmag,rimag,
     $                 lj,li,izero,aj(1,lj),ai,ylmi,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 neqilat,ilateq, wkrad )
c -->                                  vnltmp(neqilat)-s
                    endif
c
                    fac = one / (vnorm(li)*vnorm(lj))
c
                    do  k=1,nk
C                      vnlij(i,j,k) = vnlij(i,j,k) + fac*vnlr(k)
                      vnlij(k,ij) = vnlij(k,ij) + fac*vnlr(k)
                    enddo
                    if( ncplx.eq.2 .and. i.ne.j )then
                      ijc = ij + matij
                      do  k=1,nk
C                       vnlij(j,i,k) = vnlij(j,i,k) + fac*vnli(k)
                        vnlij(k,ijc) = vnlij(k,ijc) + fac*vnli(k)
                      enddo
                    endif
c
                    if( lstout.gt.4 )then
                      write(IWR,*) '      i,j,vnl=',i,j,vnlr(1),vnli(1)
                    endif
c
c                   Close angular loops:
  410             continue
  420           continue
c
  799           continue
c               Close Atom i shell loop:
                i = i0 + norbshi
  800         continue
c
c             Close Atom i lattice loop:
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
c >>>> Reassamble final matrix
c
      if( nprocs.gt.1 ) call TIMER('  nlocmat compute')
      if( Lmatreduce .and. nprocs .gt. 1 )then
c
c       Each local processor has full non-local Hamiltonian over its
c       share of atomic potentials.
c       Reduce (sum) these to master, to get total non-local Ham matrix.
c
        call MPREDUC8( master, nmat, vnlij, vnltp, icomm )
C        call MPI_Reduce( vnlij, vnltp, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocmat - Reduce vnlij' )
c
c       Put the matrix into the right place
        call DCOPY( nmat, vnltp,1, vnlij,1 )
c
        call TIMER('  nlocmat Reduce ')
c
c       The final output matrix is only on the master processor
c
      endif
c
c  Dump diagnostics for debugging if wanted:
c
      write(IWR,*) 'NLOCMAT: Max neqilat=',mxeqil
C      if( iproc .eq. master )then
C        if( lstout.gt.0 )then
C          do  k=1,nk
C            sumii = zero
C            sumsqii = zero
C            sumij = zero
C            sumsqij = zero
C            do  i=1,norb
C              sumii = sumii + vnlij(i,i,k)
C              sumsqii = sumsqii + vnlij(i,i,k)**2
C              do  j=1,norb
C                sumij = sumij + vnlij(i,j,k)
C                sumsqij = sumsqij + vnlij(i,j,k)**2
C              enddo
C            enddo
C            write(IWR,*) '  k,sumii,sqii=',k,sumii,sumsqii
C            write(IWR,*) '  k,sumij,sqij=',k,sumij,sumsqij
C          enddo
C        endif
C        if( lstout.ge.2 )then
C          write(IWR,*)     '  diagonal elements (i,i,k=1):'
C          write(IWR,9400)  (vnlij(i,i,1),i=1,norb)
C 9400     format(1x,1p6d13.4/2x,6d13.4)
C        endif
C        if( lstout.ge.3 )then
C          write(IWR,*)     '  all elements (i,j,k=1):'
C          write(IWR,9400)  ((vnlij(i,j,1),i=1,norb),j=1,norb)
C        endif
Cc       End master-only diagnostic code
C      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLOCMAT1
c
c
      subroutine NLOCMAT1( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl, ratm, rlat,coskr,sinkr, gntcomb,
     $ nrd,nrad, radmsh,vpsrad, wkrad,
     $ vnlr,vnli,
     $ vnlij,  angsav, vnltp )
c---------------------------------------------------------------
c Purpose: calculate matrix elements due to non-local part of
c          the electron:effective-nuclear-potential interaction
c
c Revision history:
c   7Sep08-PAS/2.62: switch ij-ji to natural order
c  20Jun07-PAS/2.60: merge tp and serial (2oct07: small cosmetic)
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  16Apr03-DMC: Add vnltp to parameter list, to provide a local copy
c               of the vnlij array, for parallelism.
c   2Oct02-PAS/2.54: optimize NLINT, and clean out consequences
c   7Mar02-PAS/2.52: ang constants extracted
c  21Jun01-PAS/2.48: replace STOPs
c  22May01-PAS/2.47: condense radial scratch space
c   6Mar00-PAS/2.43: rename "nl" routines, "ennloc" to "nlocmat"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Sep98-PAS/2.27: cosmetic - naming convention changed
c   5Aug98-PAS/2.24: fix cutoff bug, full proper cutoff scaling
c  21Jul98-PAS/2.23: clean up cutoffs, extend I:N terms
c  22May98-PAS/2.22: bugfix
c  28Apr97-PAS/2.17: naming conventions
c  18Jan94-PAS/2.04: fix nlat, j-indexing bugs; rationalize
c         arguments; full cleanup, match sij,tij,enij,ennloc
c  21Jun89-PAS: extract ANGSET and some of LLOOPS out of under
c               orbital loops, pack angsav
c    May89-PAS: extract RADNLOC out of under angular loops
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c  Allows for non-locality up to l=lmxnlp1(ktyp)-1
c
c  The pseudopotential for atom type ktyp is assumed to fall to
c  zero as exp(-r**2*almnnl(ktyp)).  almnnl(ktyp) here is
c  taken from the fits of Bachelet, et al.
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
c Output array:
      DIMENSION  vnlij(norb,norb,nk)
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),vpsrad(nrd,4,ntyp)
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vnlr(nk),vnli(nk)
      DIMENSION  angsav(5*5*5*5*3*81)
c Task-parallel scratch:
      DIMENSION  vnltp(*)
c
c Local declarations:
      LOGICAL    iatmi, diagshl
      LOGICAL    izero(9,9)
      DIMENSION  vnlint(5*5*15),vnltmp(5*5*15)
      DIMENSION  rnuc(3),ri(3),rj(3),risav(3)
c neqlatd in main routine should match this too.
      PARAMETER  ( neqlatd = 72 )
      DIMENSION  aj(9,9),ai(9,9,neqlatd), ilateq(neqlatd)
      DIMENSION  ylmj(2,25),ylmi(2,25,neqlatd)
c
      DIMENSION  vnorm(9)
c lstout sets level of diagnostic output, least:0-5:most
      DATA  lstout / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/NLOCMAT1: ijkmat'
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
      nmat = nk*norb*norb
c
      mxeqil = 0
      if( lstout.gt.0 ) write(IWR,*) 'NLOCMAT: v-nonlocal diagnostics'
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      if( nprocs .gt. 0 ) call MKZERO( nmat, vnlij )
c
c    Nucleus loop:
c
      do 3000 nuc=1,natm
c
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
c     Atom j loop:
c
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
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
        call AINIT( aj )
        if( rjmag .lt. 1.d-12 )then
          ljmxp1 = 1
          ylmj(1,1) = one / rt4pi
c
        else
          ljmxp1 = ljmax + lmxp1
          call AOFR( aj, rj )
          call YLMREAL( ylmj, ljmax+lmxnl, rj )
c
        endif
c
c       Atom j shell loop:
c
        j = jstart
        do 1800 jshl=1,nshlj
          nalj = nala(jshl,jtyp)
          lshlj = lshel(jshl,jtyp)
          norbshj = 2*lshlj + 1
          lj0 = lshlj*lshlj
c
          j0 = j
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
          inext = jstart
          do 1000 iatm=jatm,natm
            ityp = itypa(iatm)
            nshli = numshl(ityp)
            norbati = norba(ityp)
            istart = inext
            inext = istart + norbati
c
c           Following helps ensure only i.ge.j done:
            ishl1 = 1
            if( iatm.eq.jatm )then
              ishl1 = jshl
              do  ishl=1,ishl1-1
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                istart = istart + norbshi
              enddo
            endif
c
            alamini = alamin(ityp)
c
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
            do 900 ilatx=1,nlat2c
              risquar = risqnxt
              risqnxt = 9.d9
c
c             First test that iatm overlaps nuc ...
              if( risquar.gt.rinmx1 ) goto 1000
c
c              ... and test for 3-center overlap as well:
c               (I align the ri and rj vectors for this test)
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
              do 190 ilat=1,nlat2c
c
                do  ir=1,3
                  ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
                enddo
                risqtst = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
c               Compile orientation info of atoms such that risq=risquar
c
                if( risqtst-risquar .gt.  1.d-6 )then
c
c                 This risq is bigger, check if it is next biggest ...
c
                  if( risqtst.lt.risqnxt )then
c                    ... it indeed is, and record it:
                    risqnxt = risqtst
                    ilatnxt = ilat
                  endif
c
                elseif( ABS( risqtst-risquar ) .le. 1.d-6 )then
c
c                 This risq is within current star, process it
c
                  if( ilat.eq.1 .and. iatm.eq.jatm )then
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
c
                  call AINIT( ai(1,1,ieqlat) )
                  if( risqtst .lt. 1.d-12 )then
                    limxp1 = 1
                    ylmi(1,1,ieqlat) = one / rt4pi
c
                  else
                    limxp1 = limax + lmxp1
                    call AOFR( ai(1,1,ieqlat), ri )
                    call YLMREAL( ylmi(1,1,ieqlat), limax+lmxnl, ri )
c
                  endif
c
                endif
c
  190         continue
c
              if( iatmi )then
                if( rijn0.gt.riinmax ) goto 900
              else
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
              maxlj = (lshlj+1)*(lshlj+1)
              maxli = (ilmax+1)*(ilmax+1)
c
c             Carve up space in angsav() for arrays in angset ...
              iangs = 1
c              ... angsav():
              nangs = ljmxp1*limxp1*limxp1*lmxp1*9*9
              igntj = iangs + nangs
c              ... gntci() (complex, allow 4* real space for safety):
              ngntj = 4*(ljmxp1*maxlj*lmxp1*lmxp1)
              ignti = igntj + ngntj
c              ... cgntcj() (also complex):
              ngnti = 4*(limxp1*limxp1*maxli*lmxp1*lmxp1)
c
c             Compute sums of three-j symbols:
              call NLANG( maxlj,maxli,ljmxp1,limxp1,lmxp1,rjmag,rimag,
     $         ylmj,gntcomb,angsav,angsav(igntj),angsav(ignti) )
c
c             Atom i shell loop:
c
              i = istart
              do 800 ishl=ishl1,nshli
                nali = nala(ishl,ityp)
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                li0 = lshli*lshli
c
                i0 = i
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
                lpowmx = lshli + lshlj + 1
                lstmx = lpowmx*lmxp1
                lijmx = limxp1*ljmxp1*lstmx
                do  lij=1,lijmx
                  vnlint(lij) = zero
                enddo
c
c               Alpha loops:
c
                nx = 0
                do 320 jal=1,nalj
                  alj = ala(jal,jshl,jtyp)
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
c
c                   Test whether this pair of alphas contributes:
                    ex = rjsquar*alj*alaminn / (alj+alaminn)
                    rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                               ( (alj+alaminn)*ali )
                    if( rijn.gt.rijnmx ) goto 310
c
                    alsum = ali + alj
                    exparg = ali*risquar + alj*rjsquar
                    bargi = two*ali*rimag
                    bargj = two*alj*rjmag
                    bargs = bargi + bargj
c
c                   Evaluate radial integrals
c
                    call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $               exparg,bargs,alsum,radmsh(1,ktyp), wkrad )
                    nr0 = nrmn - 1
                    nrlen = nrmx - nrmn + 1
c
                    call NLRAD( vnltmp, alsum,bargj,bargi,bargs,
     $               exparg, ljmxp1,limxp1,lmxp1,lstmx,lpowmx,
     $               nrlen,nrd,nr0, radmsh(1,ktyp),vpsrad(1,1,ktyp),
     $               wkrad(1,1), wkrad(1,8), wkrad(1,15),
     $               wkrad(1,22),wkrad(1,23),wkrad(1,24),
     $               wkrad(1,25),wkrad(1,26),wkrad(1,27) )
c
                    fac = cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
                    do  lij=1,lijmx
                      vnlint(lij) = vnlint(lij) + fac*vnltmp(lij)
                    enddo
                    nx = nx + 1
c
c                   Close alpha loops:
  310             continue
  320           continue
c
                if( lstout.gt.3 )then
                  write(IWR,*) '   nuc,iatm,jatm,klat,ilatx,neqilat,nx'
                  write(IWR,9350)  nuc,iatm,jatm,klat,ilatx,neqilat,nx
 9350             format(1x,10i4)
                endif
c
                if( nx.eq.0 ) goto 799
c
c               Angular loops:
c
                mli1 = 1
                do 420 mlj=1,norbshj
                  lj = lj0 + mlj
                  j = j0 + mlj
c
                  if( diagshl ) mli1 = mlj
                  do 410 mli=mli1,norbshi
                    li = li0 + mli
                    i = i0 + mli
c
c                   Evaluate contribution to vnlij(i,j)
c
                    if( neqilat.eq.1 )then
                      call NLINT( vnlr,vnli, vnlint,
     $                 ljmxp1,limxp1,lmxp1,lstmx,lpowmx, rjmag,rimag,
     $                 lj,li, aj(1,lj),ai(1,li,1), ylmi,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 ilateq(1) )
                    else
c
                      do 405 lmpj=1,lj
                        if( aj(lmpj,lj).eq.zero )then
                          do  lmpi=1,li
                            izero(lmpj,lmpi) = .true.
                          enddo
                        else
                          do 404 lmpi=1,li
                            do  lateqi=1,neqilat
                              if( ai(lmpi,li,lateqi).ne.zero ) goto 403
                            enddo
                            izero(lmpj,lmpi) = .true.
                            goto 404
  403                       continue
                            izero(lmpj,lmpi) = .false.
  404                     continue
                        endif
  405                 continue
c
                      call NLINTEQ( vnlr,vnli, vnlint,
     $                 ljmxp1,limxp1,lmxp1,lstmx,lpowmx, rjmag,rimag,
     $                 lj,li,izero,aj(1,lj),ai,ylmi,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 neqilat,ilateq, wkrad )
c -->                                  vnltmp(neqilat)-s
                    endif
c
                    fac = one / (vnorm(li)*vnorm(lj))
c
                    do  k=1,nk
                      vnlij(i,j,k) = vnlij(i,j,k) + fac*vnlr(k)
                    enddo
                    if( ncplx.eq.2 .and. i.ne.j )then
                      do  k=1,nk
                        vnlij(j,i,k) = vnlij(j,i,k) + fac*vnli(k)
                      enddo
                    endif
c
                    if( lstout.gt.4 )then
                      write(IWR,*) '      i,j,vnl=',i,j,vnlr(1),vnli(1)
                    endif
c
c                   Close angular loops:
  410             continue
  420           continue
c
  799           continue
c               Close Atom i shell loop:
                i = i0 + norbshi
  800         continue
c
c             Close Atom i lattice loop:
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
c >>>> Reassamble final matrix
c
      if( nprocs .gt. 1 )then
        call TIMER('  nlocmat compute')
c
c       Each local processor has full non-local Hamiltonian over its
c       share of atomic potentials.
c       Reduce (sum) these to master, to get total non-local Ham matrix.
c
        call MPREDUC8( master, nmat, vnlij, vnltp, icomm )
C        call MPI_Reduce( vnlij, vnltp, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocmat - Reduce vnlij' )
c
c       Put the matrix into the right place
        call DCOPY( nmat, vnltp,1, vnlij,1 )
c
        call TIMER('  nlocmat Reduce ')
c
c       The final output matrix is only on the master processor
c
      endif
c
c  Dump diagnostics for debugging if wanted:
c
      write(IWR,*) 'NLOCMAT: Max neqilat=',mxeqil
C      if( iproc .eq. master )then
C        if( lstout.gt.0 )then
C          do  k=1,nk
C            sumii = zero
C            sumsqii = zero
C            sumij = zero
C            sumsqij = zero
C            do  i=1,norb
C              sumii = sumii + vnlij(i,i,k)
C              sumsqii = sumsqii + vnlij(i,i,k)**2
C              do  j=1,norb
C                sumij = sumij + vnlij(i,j,k)
C                sumsqij = sumsqij + vnlij(i,j,k)**2
C              enddo
C            enddo
C            write(IWR,*) '  k,sumii,sqii=',k,sumii,sumsqii
C            write(IWR,*) '  k,sumij,sqij=',k,sumij,sumsqij
C          enddo
C        endif
C        if( lstout.ge.2 )then
C          write(IWR,*)     '  diagonal elements (i,i,k=1):'
C          write(IWR,9400)  (vnlij(i,i,1),i=1,norb)
C 9400     format(1x,1p6d13.4/2x,6d13.4)
C        endif
C        if( lstout.ge.3 )then
C          write(IWR,*)     '  all elements (i,j,k=1):'
C          write(IWR,9400)  ((vnlij(i,j,1),i=1,norb),j=1,norb)
C        endif
Cc       End master-only diagnostic code
C      endif
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
