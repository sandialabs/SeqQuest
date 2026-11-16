c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLOCMAT
c
c
      subroutine PNLOCMAT( IWR,
     $ norb,nk, natm,ntyp,nshld,nald, nlat  ,nlat2c,
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
c  20Jun07-PAS/2.60: merge tp and serial
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
c output array:
      DIMENSION  vnlij(norb,norb,nk)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),vpsrad(nrd,4,ntyp)
c scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vnlr(nk),vnli(nk)
      DIMENSION  angsav(5*5*5*5*3*81)
c task-parallel scratch:
      DIMENSION  vnltp(norb,norb,*)
c
c local declarations:
      LOGICAL    iatmi
      LOGICAL    izero(9,9)
      DIMENSION  vnlint(5*5*15),vnltmp(5*5*15)
      DIMENSION  rnuc(3),ri(3),rj(3),rjsav(3)
c neqlatd in main routine should match this too.
      PARAMETER  (neqlatd=72)
      DIMENSION  ai(9,9),aj(9,9,neqlatd), jlateq(neqlatd)
      DIMENSION  ylmi(2,25),ylmj(2,25,neqlatd)
c
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c lstout sets level of diagnostic output, least:0-5:most
      DATA  lstout / 0 /
c---------------------------------------------------------------
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
      rt4pi = two*rtpi
c
      nmat = nk*norb*norb
c
      mxeqjl = 0
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
c       Test overlap of atom i and short-range pot due to "nuc"
        if( risquar.gt.rinmax ) goto 2000
c
        limax = 0
        do  ishl=1,nshli
          if( lshel(ishl,ityp).gt.limax ) limax = lshel(ishl,ityp)
        enddo
        call AINIT( ai )
        if( rimag.lt.1.e-12 )then
          limxp1 = 1
          ylmi(1,1) = one/rt4pi
c
        else
          limxp1 = limax + lmxp1
          call AOFR( ai, ri )
          call YLMREAL( ylmi, limax+lmxnl, ri )
c
        endif
c
        i = istart
        do 1800 ishl=1,nshli
          nali = nala(ishl,ityp)
          lshli = lshel(ishl,ityp)
          norbshi = 2*lshli + 1
          li0 = lshli*lshli
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
c         Test whether this radial function will contribute:
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
c           Following helps ensure only j.le.i done:
            if( jatm.eq.iatm ) nshlj = ishl
c
            alaminj = alamin(jtyp)
c
c           Criterion for including contribution of atom j:
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
c           Find nearest jatm image:
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
            do 900 jlatx=1,nlat2c
              rjsquar = rjsqnxt
              rjsqnxt = 9.d9
c
c             First test that jatm overlaps nuc ...
              if( rjsquar.gt.rjnmx1 ) goto 1000
c
c              ... and test for 3-center overlap as well:
c               (I align the ri and rj vectors for this test)
              rjmag = SQRT( rjsquar )
              if( rjsquar.gt.riovlp )then
                rjnmx2 = rijnmx0 * ( (alamini+alaminn) /
     $           (alamini+alaminn-alamini*rimag/rjmag) )**2
                if( rjsquar.gt.rjnmx2 ) goto 1000
              endif
c
c             Align the jatm test vector along the iatm vector
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
              do 190 jlat=1,nlat2c
c
                do  jr=1,3
                  rj(jr) = ratm(jr,jatm) + rlat(jr,jlat) - rnuc(jr)
                enddo
                rjsqtst = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
c               Compile orientation info of atoms such that rjsq=rjsquar
c
                if( rjsqtst-rjsquar .gt.  1.d-6 )then
c
c                 This rjsq is bigger, check if it is next biggest ...
c
                  if( rjsqtst.lt.rjsqnxt )then
c                    ... it indeed is, and record it:
                    rjsqnxt = rjsqtst
                    jlatnxt = jlat
                  endif
c
                elseif( ABS( rjsqtst-rjsquar ) .le. 1.d-6 )then
c
c                 This rjsq is within current star, process it
c
                  if( jlat.eq.1 .and. jatm.eq.iatm )then
                    iatmi = .true.
                  else
c                   Check tighter cutoff ...
                    if( rijn0.gt.rijnmax ) goto 190
                  endif
c
                  jeq = jeq + 1
                  if( jeq.gt.neqlatd ) goto 190
c
                  jlateq(jeq) = jlat
c
                  call AINIT( aj(1,1,jeq) )
                  if( rjsqtst .lt. 1.d-12 )then
                    ljmxp1 = 1
                    ylmj(1,1,jeq) = one/rt4pi
c
                  else
                    ljmxp1 = ljmax + lmxp1
                    call AOFR( aj(1,1,jeq), rj )
                    call YLMREAL( ylmj(1,1,jeq), ljmax+lmxnl, rj )
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
              neqjlat = jeq
              if( neqjlat.gt.mxeqjl ) mxeqjl = neqjlat
              if( neqjlat.gt.neqlatd )then
                write(IWR,*) '>>>>> neqjlat,neqlatd=',neqjlat,neqlatd
                call STOPXERR( 'neqlatd - exceeded in NLOCFRC' )
              endif
              if( neqjlat.le.0 ) goto 900
c
c             Find largest l for atom j with alpha that will contribute:
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
              maxli = (lshli+1)*(lshli+1)
              maxlj = (jlmax+1)*(jlmax+1)
c
c             Carve up space in angsav() for arrays in angset ...
              iangs = 1
c              ... angsav():
              nangs = limxp1*ljmxp1*ljmxp1*lmxp1*9*9
              ignti = iangs + nangs
c              ... gntci() (complex, allow 4* real space for safety):
              ngnti = 4*(limxp1*maxli*lmxp1*lmxp1)
              igntj = ignti + ngnti
c              ... cgntcj() (also complex):
              ngntj = 4*(ljmxp1*ljmxp1*maxlj*lmxp1*lmxp1)
c
c             Compute sums of three-j symbols:
              call NLANG( maxli,maxlj,limxp1,ljmxp1,lmxp1,rimag,rjmag,
     $         ylmi,gntcomb,angsav,angsav(ignti),angsav(igntj) )
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
                lpowmx = lshli + lshlj + 1
                lstmx = lpowmx*lmxp1
                lijmx = limxp1*ljmxp1*lstmx
                do  lij=1,lijmx
                  vnlint(lij) = zero
                enddo
c
c    Alpha loops:
c
                nx = 0
                do 320 ial=1,nali
                  ali = ala(ial,ishl,ityp)
c
c                 Test whether this alpha contributes:
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
c
c                   Test whether this pair of alphas contributes:
                    ex = risquar*ali*alaminn / (ali+alaminn)
                    rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                               ( (ali+alaminn)*alj )
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
                    call NLRAD( vnltmp, alsum,bargi,bargj,bargs,
     $               exparg, limxp1,ljmxp1,lmxp1,lstmx,lpowmx,
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
c    Close alpha loops:
c
  310             continue
  320           continue
c
                if( lstout.gt.3 )then
                  write(IWR,*) '   nuc,iatm,jatm,klat,jlatx,neqjlat,nx'
                  write(IWR,9350)  nuc,iatm,jatm,klat,jlatx,neqjlat,nx
 9350             format(1x,10i4)
                endif
c
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
c
c                   Evaluate contribution to vnlij(i,j)
c
                    if( neqjlat.eq.1 )then
                      call NLINT( vnlr,vnli, vnlint,
     $                 limxp1,ljmxp1,lmxp1,lstmx,lpowmx, rimag,rjmag,
     $                 li,lj, ai(1,li),aj(1,lj,1), ylmj,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 jlateq(1) )
                    else
c
                      do 405 lmpi=1,li
                        if( ai(lmpi,li).eq.zero )then
                          do  lmpj=1,lj
                            izero(lmpi,lmpj) = .true.
                          enddo
                        else
                          do 404 lmpj=1,lj
                            do  jeq=1,neqjlat
                              if( aj(lmpj,lj,jeq).ne.zero ) goto 403
                            enddo
                            izero(lmpi,lmpj) = .true.
                            goto 404
  403                       izero(lmpi,lmpj) = .false.
  404                     continue
                        endif
  405                continue
c
                      call NLINTEQ( vnlr,vnli, vnlint,
     $                 limxp1,ljmxp1,lmxp1,lstmx,lpowmx, rimag,rjmag,
     $                 li,lj,izero,ai(1,li),aj,ylmj,angsav,
     $                 nk,nlat, coskr,sinkr,
     $                 neqjlat,jlateq, wkrad )
c -->                                  vnltmp(neqjlat)-s
                    endif
c
                    fac = one / (vnorm(li)*vnorm(lj))
                    do  k=1,nk
                      vnlij(i,j,k) = vnlij(i,j,k) + fac*vnlr(k)
                      if( i.ne.j )
     $                 vnlij(j,i,k) = vnlij(j,i,k) - fac*vnli(k)
                    enddo
                    if( lstout.gt.4 )then
                      write(IWR,*) '      i,j,vnl=',i,j,vnlr(1),vnli(1)
                    endif
c
c    Close angular loops:
c
  410             continue
  420           continue
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
c
      if( nprocs .gt. 0 )then
        call TIMER( 'nlocmat compute' )
c
c       Each local processor has full non-local Hamiltonian over its
c       share of atomic potentials.
c       Need to sum these to master, to get total non-local Ham matrix.
c
        call MPREDUC8( master, nmat, vnlij, vnltp, icomm )
C        call MPI_Reduce( vnlij, vnltp, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocmat - Reduce vnlij' )
c
c       Put the matrix into the right place
        call DCOPY( nmat, vnltp,1, vnlij,1 )
Cc
Cc       For now, we will also broadcast the total to all processors.
C        call MPBCAST8( master, nmat, vnlij, icomm )
C        call MPI_Bcast( vnlij, nmat, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pnlocmat - Bcast vnlij' )
      endif
c
c  Dump diagnostics for debugging if wanted:
c
      write(IWR,*) 'NLOCMAT: Max neqjlat=',mxeqjl
      if( iproc .eq. master )then
        if( lstout.gt.0 )then
          do  k=1,nk
            sumii = zero
            sumsqii = zero
            sumij = zero
            sumsqij = zero
            do  i=1,norb
              sumii = sumii + vnlij(i,i,k)
                sumsqii = sumsqii + vnlij(i,i,k)**2
              do  j=1,norb
                sumij = sumij + vnlij(i,j,k)
                sumsqij = sumsqij + vnlij(i,j,k)**2
              enddo
            enddo
            write(IWR,*) '  k,sumii,sqii=',k,sumii,sumsqii
            write(IWR,*) '  k,sumij,sqij=',k,sumij,sumsqij
          enddo
        endif
        if( lstout.ge.2 )then
          write(IWR,*)     '  diagonal elements (i,i,k=1):'
          write(IWR,9400)  (vnlij(i,i,1),i=1,norb)
 9400     format(1x,1p6d13.4/2x,6d13.4)
        endif
        if( lstout.ge.3 )then
          write(IWR,*)     '  all elements (i,j,k=1):'
          write(IWR,9400)  ((vnlij(i,j,1),i=1,norb),j=1,norb)
        endif
c       End master-only diagnostic code
      endif
c
c    That's all Folks!
c
      RETURN
      END
