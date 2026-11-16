c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCMAT
c
c
ctp:
      subroutine PVLOCMAT( IWR,
     $ norb,nk, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vij, vloc1,vloc2, vloc3 )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all but the 2-center I:N terms (done in "vlocmii")
c
c Written: P.A. Schultz, for v2.43, based on the routine "enij"
c          originally written by P.J. Feibelman
c
c Revision history:
c  18Jun07-PAS/2.60: merge tp and serial
c  29Jan07-PAS/2.60: clean up tp; put barrier between Reduces
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines.
c   DMC:  task parallel
c   7Mar02-PAS/2.52: ang constants extracted
c  22May01-PAS/2.47: condense radial scratch space
c   6Mar00-PAS/2.43: rename "vl" routines, "enij" to "vlocmat"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  15Sep99-PAS/2.39: pack big workspace down + simple opt
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Sep98-PAS/2.27: compact temp orbital matrix
c  27Aug98-PAS/2.25: purge (N**7) triplet equivalences
c  11Aug98-PAS/2.24: match 3-center routines again
c  16Jul98-PAS/2.23: slow+occ cutoffs finally dealt with properly,
c                    split off I:N into separate routine (ENII),
c                    accomodate split basis
c  28Apr97-PAS/2.17: naming conventions
c  13NOv96-PAS/2.16: opt k-loops, k-first temp Ham arrays
c  29Sep95-PAS/2.13: fix F(N-II) shell-pruning, tighten cutoffs,
c                    extract k-dependent loops from bottom for opt
c   3Jun94-PAS/2.10: exc matrix install
c  18Jan94-PAS/2.04: full cleanup, match sij,tij,enij,ennloc
c  21Jun89-PAS: reorganized to extract SHTRNG out from under
c               angular loops
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c
c  Note that vij() matrix is incremented while the vxc/exc arrays
c  are effectively initialized here and are put out in vloc1 and vloc2.
c  They are output in vlocN(i,j,k) order, but are used as scratch
c  space in kij order, requiring some fast stepping.  If nk=1,
c  this order is identical, and can do some shortcuts.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c output arrays:
      DIMENSION  vij(norb,norb,nk)
      DIMENSION  vloc1(nk,norb,norb),vloc2(nk,norb,norb)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $ nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  znuc(ntyp), occsh(nshld,ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  vatrad(nrd,ntyp)
      DIMENSION  vesrad(nrd,ntyp),vxcrad(nrd,ntyp),excrad(nrd,ntyp)
c scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vloc3(nk,norb,norb)
c
c local declarations:
c
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      LOGICAL    keepsh
      DIMENSION  rnuc(3),ri(3),rj(3), p(3), ci(10,9),cj(10,9)
      DIMENSION  vatint(5),vesint(5),vxcint(5),excint(5)
c
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
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
c
      emxslo = emaxsl
c
      fourpi = pi + pi + pi + pi
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      nmat = nk*norb*norb
c
      if( lstout.ne.0 )  write(IWR,9000)
 9000 format(/3x,'nuc  klat iatm',9x,'E(Z-Z)',12x,'E(E-Z)',12x,'TOT')
 9001 format(1x,'ENRI',3i3,2x,3f20.10)
 9002 format(5x,a,3i3,2x,f20.10)
 9003 format(5x,a,4i3,2x,f20.10)
c
      call MKZERO( nmat, vloc1 )
      call MKZERO( nmat, vloc2 )
      call MKZERO( nmat, vloc3 )
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
      alaminn = two*alamin(ktyp)
      almin = 1.d9
      do  ishl=1,nocc(ktyp)
        occshl = occsh(ishl,ktyp)
        if( occshl.ne.zero )then
          do  ial=1,nalsh(ishl,ktyp)
            if( alsh(ial,ishl,ktyp).lt.almin ) almin=alsh(ial,ishl,ktyp)
          enddo
        endif
      enddo
      alaminn = two*almin
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
      call CINIT( ci )
      call COFR( ci,ri )
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
      emax = emaxf
      if( almaxi.lt.alfast ) emax = emaxsl
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
        if( jatm.eq.iatm )then
          nshlj = ishl
          jlat1 = 2
        else
          jlat1 = 1
        endif
c
        alaminj = alamin(jtyp)
c  criterion for including contribution of atom j:
        ex = risquar*alamini*alaminn / (alamini+alaminn)
        rijnmax = (emaxsl-ex)*(alamini+alaminj+alaminn) /
     $                      ( (alamini+alaminn)*alaminj )
c
c    Lattice vector loop:
c
        do 900 jlat=jlat1,nlat2c
c
          do  jr=1,3
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
          if( rijn0.gt.rijnmax ) goto 900
c
          call CINIT( cj )
          call COFR( cj,rj )
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
c
            ex = risquar*alamini*alaminn / (alamini+alaminn)
            rijnmx = (emax-ex)*(alamini+alminj+alaminn) /
     $                       ( (alamini+alaminn)*alminj )
            if( rijn0.gt.rijnmx ) goto 800
c
            lmx = lshli + lshlj + 1
            call MKZERO( 25, vessh )
            call MKZERO( 25, vijsh )
            call MKZERO( 25, vxcsh )
            call MKZERO( 25, excsh )
c
c    Alpha loops:
c
            enijsum = zero
            coulsum = zero
            do 320 ial=1,nali
              ali = ala(ial,ishl,ityp)
c
c  Test whether this alpha contributes:
              rinmx = emax*(ali+alaminn) / (ali*alaminn)
              if( risquar.gt.rinmx ) goto 320
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
                if( rijn.gt.rijnmx ) goto 310
c
                alsum = ali + alj
                exparg = ali*risquar + alj*rjsquar
                do  ir=1,3
                  p(ir) = two*(ali*ri(ir)+alj*rj(ir))
                enddo
                pmag = SQRT( p(1)*p(1)+p(2)*p(2)+p(3)*p(3) )
c
                call MKZERO( 5, vatint )
                call MKZERO( 5, vesint )
                call MKZERO( 5, vxcint )
                call MKZERO( 5, excint )
c
                nonzer = 0
c
c   "vlrad" computes the contributions of the xc-potentials
c   associated with atom nuc, as well as the difference between
c   the actual electrostatic potential due to the electrons of
c   atom nuc and 2.*abs(znuc)/r
c
                call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $           exparg,pmag,alsum,radmsh(1,ktyp),wkrad )
c
                if( nrmx.gt.0 )then
c
c    Evaluate radial integrals
c
                  nrlen = nrmx - nrmn + 1
                  call VLRAD( vatint, vesint,vxcint,excint,
     $             lmx,alsum,pmag,exparg,
     $             nrlen, radmsh(nrmn,ktyp),radwt(nrmn,ktyp),
     $             vatrad(nrmn,ktyp),vesrad(nrmn,ktyp),
     $             vxcrad(nrmn,ktyp),excrad(nrmn,ktyp),
     $             wkrad, wkrad(1, 8),wkrad(1, 9),wkrad(1,10),
     $                    wkrad(1,11),wkrad(1,12),wkrad(1,13) )
                else
                  if( nonzer.eq.0 ) goto 310
                endif
c
c    Angular loops:
c
                shlfac = fourpi*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
                norbshj = norbshj0
                i = i0
                do 420 li=li0+1,li0+norbshi
                  i = i + 1
                  lis = li - li0
                  norbshj = norbshj + joff
                  j = j0
                  do 410 lj=lj0+1,lj0+norbshj
                    j = j + 1
                    ljs = lj - lj0
c
c    Evaluate contribution to vij(i,j)
c
                    call VLANG( vatij,vatint,
     $               vesij,vesint, vxcij,vxcint, excij,excint,
     $               p,li,lj,ci(1,li),cj(1,lj) )
c
                    fac = shlfac / ( vnorm(li)*vnorm(lj) )
c
c    Multiply by bf factors, and increment shell-block Hamiltonian
                    vesij = fac*vesij
                    vatij = fac*vatij
                    vxcij = fac*vxcij
                    excij = fac*excij
                    vessh(lis,ljs) = vessh(lis,ljs) + vesij
                    vijsh(lis,ljs) = vijsh(lis,ljs) + vatij
                    vxcsh(lis,ljs) = vxcsh(lis,ljs) + vxcij
                    excsh(lis,ljs) = excsh(lis,ljs) + excij
c
c    Close angular loops:
c
  410             continue
  420           continue
c
c    Close alpha loops:
c
  310         continue
  320       continue
c
            norbshj = norbshj0
            i = i0
            do 425 li=li0+1,li0+norbshi
              i = i + 1
              lis = li - li0
              norbshj = norbshj + joff
              j = j0
              do 415 lj=lj0+1,lj0+norbshj
                j = j + 1
                ljs = lj - lj0
                vesij = vessh(lis,ljs)
                vatij = vijsh(lis,ljs)
                vxcij = vxcsh(lis,ljs)
                excij = excsh(lis,ljs)
                do  k=1,nk
                  vloc1(k,j,i) = vloc1(k,j,i) + vatij*coskr(k,jlat)
                enddo
                do  k=1,nk
                  vloc2(k,j,i) = vloc2(k,j,i) + vxcij*coskr(k,jlat)
                enddo
                do  k=1,nk
                  vloc3(k,j,i) = vloc3(k,j,i) + excij*coskr(k,jlat)
                enddo
c on dec&sun, sep loops faster (stride?) - 13nov96-pas
                if( i.ne.j )then
                  do  k=1,nk
                    vloc1(k,i,j) = vloc1(k,i,j) - vatij*sinkr(k,jlat)
                  enddo
                  do  k=1,nk
                    vloc2(k,i,j) = vloc2(k,i,j) - vxcij*sinkr(k,jlat)
                  enddo
                  do  k=1,nk
                    vloc3(k,i,j) = vloc3(k,i,j) - excij*sinkr(k,jlat)
                  enddo
                endif
c
  415         continue
  425       continue
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
c    Close nucleus lattice loop:
c
 3000 continue
c
      if( nprocs .gt. 1 )then
        call TIMER( 'vlocmat compute' )
c
c       Each local processor has full Hamiltonian over its share
c       of atomic potentials.
c       Need to sum these to master to get total local Hams.
c       Use a broadcast from master to create barrier between reduces.
c       For now, broadcast total Ham to all processors when done;
c       later, since I/O restricted to master, could dispense with that.
c
c       Use the vij() as buffer for all reduces.
c       Use a broadcast between reduces to prevent pile-up at master
c
        call MPREDUC8( master, nmat, vloc1, vij, icomm )
C        tckt = master
C        lenmsg = 1
C        call MPBCAST8( master, lenmsg, tckt, icomm )
C        call MPI_Reduce( vloc1, vij, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc1' )
c
        call DCOPY( nmat, vij,1, vloc1,1 )
c
        call MPREDUC8( master, nmat, vloc2, vij, icomm )
C        call MPBCAST8( master, 1, tckt, icomm )
C        call MPI_Reduce( vloc2, vij, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc2' )
c
        call DCOPY( nmat, vij,1, vloc2,1 )
c
        call MPREDUC8( master, nmat, vloc3, vij, icomm )
C        call MPBCAST8( master, 1, tckt, icomm )
C        call MPI_Reduce( vloc3, vij, nmat, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc3' )
c
        call DCOPY( nmat, vij,1, vloc3,1 )
c
        call TIMER( 'vlocmat Reduce ' )
Cc
Cc       And broadcast matrices to all processors
Cc
C        call MPBCAST8( master, nmat, vloc1, icomm )
C        call MPI_Bcast( vloc1, nmat, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Bcast vloc1' )
Cc
C        call MPBCAST8( master, nmat, vloc2, icomm )
C        call MPI_Bcast( vloc2, nmat, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Bcast vloc2' )
Cc
C        call MPBCAST8( master, nmat, vloc3, icomm )
C        call MPI_Bcast( vloc3, nmat, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Bcast vloc3' )
Cc
C        call TIMER( 'vlocmat Bcast  ' )
      endif
c
c Pack up contributions to Hamiltonians from local nuc potentials
c We are sorting into (ijk) from (kji) here
c
      call VLJISRT( norb,nk, vij  , vloc1 )
      call VLJISRT( norb,nk, vloc1, vloc2 )
      call VLJISRT( norb,nk, vloc2, vloc3 )
c
c    That's all Folks!
c
      RETURN
      END
