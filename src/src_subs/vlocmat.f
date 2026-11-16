c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCMAT
c
c
      subroutine VLOCMAT( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all but the 2-center I:N terms (done in "vlocmii")
c
c Written: P.A. Schultz, for v2.43, based on the routine "enij"
c          originally written by P.J. Feibelman
c
c Revision history:
c  17Oct08-PAS/2.62: stripe-distributed-parallel
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
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  vloc0(nk,*),vxc0(nk,*),exc0(nk,*)
      DIMENSION  vtmp(nk,*)
c Input arrays:
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
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
c
c Local declarations:
c
      LOGICAL    Lmatreduce
      LOGICAL    keepsh, diagshl
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      DIMENSION  rnuc(3),ri(3),rj(3), p(3), ci(10,9),cj(10,9)
      DIMENSION  vatint(5),vesint(5),vxcint(5),excint(5)
c
      DIMENSION  vnorm(9)
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/VLOCMAT: stripe-parallel'
      Lmatreduce = .false.
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
      nmat = nk*norb*norb
      nmatpr = mymatpr*ncplx*nk
c
      if( lstout.ne.0 )  write(IWR,9000)
 9000 format(/3x,'nuc  klat iatm',9x,'E(Z-Z)',12x,'E(E-Z)',12x,'TOT')
 9001 format(1x,'ENRI',3i3,2x,3f20.10)
 9002 format(5x,a,3i3,2x,f20.10)
 9003 format(5x,a,4i3,2x,f20.10)
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
c
      alaminj = alamin(jtyp)
      emax = emxslo
      rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c     Test overlap of Atom j and short-range pot due to "nuc"
      if( rjsquar.gt.rjnmax ) goto 2000
c
      call CINIT( cj )
      call COFR( cj,rj )
c
c     Atom j shell loop:
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
      emax = emaxf
      if( almaxj.lt.alfast ) emax = emaxsl
      rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c     Test whether this radial function will contribute:
      if( rjsquar.gt.rjnmax ) goto 1799
c
c     Atom i loop:
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
c       Following helps ensure only i.ge.j done:
        iorbatn0 = 0
        ilat1 = 1
        ishl1 = 1
        if( iatm.eq.jatm )then
          ilat1 = 2
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
c       Criterion for including contribution of Atom j:
        ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
        rijnmax = (emaxsl-ex)*(alaminj+alamini+alaminn) /
     $                      ( (alaminj+alaminn)*alamini )
c
c       Lattice vector loop for Atom i:
c
        do 900 ilat=ilat1,nlat2c
c
          do  ir=1,3
            ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
          enddo
          risquar = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
          f1 = alaminj / (alaminj+alaminn)
          rijn0 = zero
          do  ir=1,3
            rijn0 = rijn0 + (f1*rj(ir)-ri(ir))**2
          enddo
c
c         Test for significant overlap of Atom j with i and nuc
          if( rijn0.gt.rijnmax ) goto 900
c
          call CINIT( ci )
          call COFR( ci,ri )
c
c         Atom i shell loop:
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
            if( almaxi.lt.alfast .and. almaxj.lt.alfast ) emax = emaxsl
c
            ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
            rijnmx = (emax-ex)*(alaminj+almini+alaminn) /
     $                       ( (alaminj+alaminn)*almini )
            if( rijn0.gt.rijnmx ) goto 799
c
            lmx = lshli + lshlj + 1
            call MKZERO( 25, vessh )
            call MKZERO( 25, vijsh )
            call MKZERO( 25, vxcsh )
            call MKZERO( 25, excsh )
c
c           Alpha loops:
c
            enijsum = zero
            coulsum = zero
            do 320 jal=1,nalj
              alj = ala(jal,jshl,jtyp)
c
c             Test whether this alpha contributes:
              rjnmx = emax*(alj+alaminn) / (alj*alaminn)
              if( rjsquar.gt.rjnmx ) goto 320
c
              fj1 = alj / (alj+alaminn)
              rijn = zero
              do  ir=1,3
                rijn = rijn + (fj1*rj(ir)-ri(ir))**2
              enddo
c
              do 310 ial=1,nali
                ali = ala(ial,ishl,ityp)
c
c               Test whether this pair of alphas contributes:
                ex = rjsquar*alj*alaminn / (alj+alaminn)
                rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                           ( (alj+alaminn)*ali )
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
c               "vlrad" computes the contributions of the xc-potentials
c               associated with atom nuc, as well as the difference between
c               the actual electrostatic potential due to the electrons of
c               atom nuc and 2.*abs(znuc)/r
c
                call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $           exparg,pmag,alsum,radmsh(1,ktyp),wkrad )
c
                if( nrmx.gt.0 )then
c
c                 Evaluate radial integrals
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
c               Angular loops:
c
                shlfac = fourpi*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
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
c                   Evaluate contribution to matrix element (i,j)
c
                    call VLANG( vatij,vatint,
     $               vesij,vesint, vxcij,vxcint, excij,excint,
     $               p,lj,li,cj(1,lj),ci(1,li) )
c
                    fac = shlfac / ( vnorm(li)*vnorm(lj) )
c
c                   Multiply by bf factors, increment shell-block Hamiltonian
                    vesij = fac*vesij
                    vatij = fac*vatij
                    vxcij = fac*vxcij
                    excij = fac*excij
                    vessh(mli,mlj) = vessh(mli,mlj) + vesij
                    vijsh(mli,mlj) = vijsh(mli,mlj) + vatij
                    vxcsh(mli,mlj) = vxcsh(mli,mlj) + vxcij
                    excsh(mli,mlj) = excsh(mli,mlj) + excij
c
c                   Close angular loops:
  410             continue
  420           continue
c
c               Close alpha loops:
  310         continue
  320       continue
c
c           Angular loops:
c
            ij0 = ijstart + jorbat0*norbati + iorbat0
            ij0 = ij0 - norbati
            mli1 = 1
            do 425 mlj=1,norbshj
              jorbat = jorbat0 + mlj
              j = j0 + mlj
c
              ij0 = ij0 + norbati
              if( diagshl ) mli1 = mlj
              do 415 mli=mli1,norbshi
                iorbat = iorbat0 + mli
                i = i0 + mli
                ij = ij0 + mli
c
                vesij = vessh(mli,mlj)
                vatij = vijsh(mli,mlj)
                vxcij = vxcsh(mli,mlj)
                excij = excsh(mli,mlj)
c
                do  k=1,nk
C                  vloc0(k,i,j) = vloc0(k,i,j) + vatij*coskr(k,ilat)
                  vloc0(k,ij) = vloc0(k,ij) + vatij*coskr(k,ilat)
                enddo
                do  k=1,nk
C                  vxc0(k,i,j) = vxc0(k,i,j) + vxcij*coskr(k,ilat)
                  vxc0(k,ij) = vxc0(k,ij) + vxcij*coskr(k,ilat)
                enddo
                do  k=1,nk
C                  exc0(k,i,j) = exc0(k,i,j) + excij*coskr(k,ilat)
                  exc0(k,ij) = exc0(k,ij) + excij*coskr(k,ilat)
                enddo
c               On dec&sun, sep loops faster (stride?) - 13nov96-pas
                if( ncplx.eq.2 .and. i.ne.j )then
                  ijc = ij + matij
                  do  k=1,nk
C                    vloc0(k,j,i) = vloc0(k,j,i) + vatij*sinkr(k,ilat)
                    vloc0(k,ijc) = vloc0(k,ijc) + vatij*sinkr(k,ilat)
                  enddo
                  do  k=1,nk
C                    vxc0(k,j,i) = vxc0(k,j,i) + vxcij*sinkr(k,ilat)
                    vxc0(k,ijc) = vxc0(k,ijc) + vxcij*sinkr(k,ilat)
                  enddo
                  do  k=1,nk
C                    exc0(k,j,i) = exc0(k,j,i) + excij*sinkr(k,ilat)
                    exc0(k,ijc) = exc0(k,ijc) + excij*sinkr(k,ilat)
                  enddo
                endif
c
c               Close angular loops:
  415         continue
  425       continue
c
  799       continue
c           Close Atom i shell loop:
            i = i0 + norbshi
  800     continue
c
c         Close Atom i lattice loop:
  900   continue
c
c       Close Atom i loop:
 1000 continue
c
 1799   continue
c       Close Atom j shell loop:
        j = j0 + norbshj
 1800 continue
c
c       Close Atom j loop:
c
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
      if( nprocs.gt.1 )  call TIMER( '  vlocmat compute' )
c
c >>>> Reassamble final matrix
c
      if( Lmatreduce .and. nprocs .gt. 1 )then
c       Each local processor has full Hamiltonian over its share
c       of atomic potentials.
c       Reduce (sum) these to master to get total Hams.
c
c       Use the vtmp() as buffer for all reduces.
c
c       Reduce full local potential matrix:
        call MPREDUC8( master, nmat, vloc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, vloc0,1 )
c
c       Reduce Vxc0 matrix:
        call MPREDUC8( master, nmat, vxc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, vxc0,1 )
c
c       Reduce Exc0 matrix:
        call MPREDUC8( master, nmat, exc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, exc0,1 )
c
        call TIMER( '  vlocmat Reduce ' )
c
c       The final output matrices are only on the master processor
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCMAT1
c
c
      subroutine VLOCMAT1( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all but the 2-center I:N terms (done in "vlocmii")
c
c Written: P.A. Schultz, for v2.43, based on the routine "enij"
c          originally written by P.J. Feibelman
c
c Revision history:
c   7Sep08-PAS/2.62: switch ij-ji to natural order
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
c  are effectively initialized here and are put out in vloc0 and vxc0.
c  They are output in vlocN(i,j,k) order, but are used as scratch
c  space in kij order, requiring some fast stepping.  If nk=1,
c  this order is identical, and can do some shortcuts.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  vloc0(nk,norb,norb)
      DIMENSION  vxc0(nk,norb,norb),exc0(nk,norb,norb)
      DIMENSION  vtmp(norb,norb,nk)
c Input arrays:
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
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
c
c Local declarations:
c
      LOGICAL    keepsh, diagshl
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      DIMENSION  rnuc(3),ri(3),rj(3), p(3), ci(10,9),cj(10,9)
      DIMENSION  vatint(5),vesint(5),vxcint(5),excint(5)
c
      DIMENSION  vnorm(9)
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/VLOCMAT1: ijkmat'
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
      call MKZERO( nmat, vloc0 )
      call MKZERO( nmat, vxc0 )
      call MKZERO( nmat, exc0 )
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
c
      alaminj = alamin(jtyp)
      emax = emxslo
      rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c     Test overlap of Atom j and short-range pot due to "nuc"
      if( rjsquar.gt.rjnmax ) goto 2000
c
      call CINIT( cj )
      call COFR( cj,rj )
c
c     Atom j shell loop:
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
      emax = emaxf
      if( almaxj.lt.alfast ) emax = emaxsl
      rjnmax = emax*(alaminj+alaminn) / (alaminj*alaminn)
c
c     Test whether this radial function will contribute:
      if( rjsquar.gt.rjnmax ) goto 1799
c
c     Atom i loop:
c
      inext = jstart
      do 1000 iatm=jatm,natm
        ityp = itypa(iatm)
        nshli = numshl(ityp)
        norbati = norba(ityp)
        istart = inext
        inext = istart + norbati
c
c       Following helps ensure only i.ge.j done:
        ilat1 = 1
        ishl1 = 1
        if( iatm.eq.jatm )then
          ilat1 = 2
          ishl1 = jshl
          do  ishl=1,ishl1-1
            lshli = lshel(ishl,ityp)
            norbshi = 2*lshli + 1
            istart = istart + norbshi
          enddo
        endif
c
        alamini = alamin(ityp)
c       Criterion for including contribution of Atom j:
        ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
        rijnmax = (emaxsl-ex)*(alaminj+alamini+alaminn) /
     $                      ( (alaminj+alaminn)*alamini )
c
c       Lattice vector loop for Atom i:
c
        do 900 ilat=ilat1,nlat2c
c
          do  ir=1,3
            ri(ir) = ratm(ir,iatm) + rlat(ir,ilat) - rnuc(ir)
          enddo
          risquar = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
          f1 = alaminj / (alaminj+alaminn)
          rijn0 = zero
          do  ir=1,3
            rijn0 = rijn0 + (f1*rj(ir)-ri(ir))**2
          enddo
c
c         Test for significant overlap of Atom j with i and nuc
          if( rijn0.gt.rijnmax ) goto 900
c
          call CINIT( ci )
          call COFR( ci,ri )
c
c         Atom i shell loop:
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
            if( almaxi.lt.alfast .and. almaxj.lt.alfast ) emax = emaxsl
c
            ex = rjsquar*alaminj*alaminn / (alaminj+alaminn)
            rijnmx = (emax-ex)*(alaminj+almini+alaminn) /
     $                       ( (alaminj+alaminn)*almini )
            if( rijn0.gt.rijnmx ) goto 799
c
            lmx = lshli + lshlj + 1
            call MKZERO( 25, vessh )
            call MKZERO( 25, vijsh )
            call MKZERO( 25, vxcsh )
            call MKZERO( 25, excsh )
c
c           Alpha loops:
c
            enijsum = zero
            coulsum = zero
            do 320 jal=1,nalj
              alj = ala(jal,jshl,jtyp)
c
c             Test whether this alpha contributes:
              rjnmx = emax*(alj+alaminn) / (alj*alaminn)
              if( rjsquar.gt.rjnmx ) goto 320
c
              fj1 = alj / (alj+alaminn)
              rijn = zero
              do  ir=1,3
                rijn = rijn + (fj1*rj(ir)-ri(ir))**2
              enddo
c
              do 310 ial=1,nali
                ali = ala(ial,ishl,ityp)
c
c               Test whether this pair of alphas contributes:
                ex = rjsquar*alj*alaminn / (alj+alaminn)
                rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                           ( (alj+alaminn)*ali )
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
c               "vlrad" computes the contributions of the xc-potentials
c               associated with atom nuc, as well as the difference between
c               the actual electrostatic potential due to the electrons of
c               atom nuc and 2.*abs(znuc)/r
c
                call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $           exparg,pmag,alsum,radmsh(1,ktyp),wkrad )
c
                if( nrmx.gt.0 )then
c
c                 Evaluate radial integrals
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
c               Angular loops:
c
                shlfac = fourpi*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
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
c                   Evaluate contribution to matrix element (i,j)
c
                    call VLANG( vatij,vatint,
     $               vesij,vesint, vxcij,vxcint, excij,excint,
     $               p,lj,li,cj(1,lj),ci(1,li) )
c
                    fac = shlfac / ( vnorm(li)*vnorm(lj) )
c
c                   Multiply by bf factors, increment shell-block Hamiltonian
                    vesij = fac*vesij
                    vatij = fac*vatij
                    vxcij = fac*vxcij
                    excij = fac*excij
                    vessh(mli,mlj) = vessh(mli,mlj) + vesij
                    vijsh(mli,mlj) = vijsh(mli,mlj) + vatij
                    vxcsh(mli,mlj) = vxcsh(mli,mlj) + vxcij
                    excsh(mli,mlj) = excsh(mli,mlj) + excij
c
c                   Close angular loops:
  410             continue
  420           continue
c
c               Close alpha loops:
  310         continue
  320       continue
c
c           Angular loops:
c
            mli1 = 1
            do 425 mlj=1,norbshj
              j = j0 + mlj
c
              if( diagshl ) mli1 = mlj
              do 415 mli=mli1,norbshi
                i = i0 + mli
c
                vesij = vessh(mli,mlj)
                vatij = vijsh(mli,mlj)
                vxcij = vxcsh(mli,mlj)
                excij = excsh(mli,mlj)
c
                do  k=1,nk
                  vloc0(k,i,j) = vloc0(k,i,j) + vatij*coskr(k,ilat)
                enddo
                do  k=1,nk
                  vxc0(k,i,j) = vxc0(k,i,j) + vxcij*coskr(k,ilat)
                enddo
                do  k=1,nk
                  exc0(k,i,j) = exc0(k,i,j) + excij*coskr(k,ilat)
                enddo
c               On dec&sun, sep loops faster (stride?) - 13nov96-pas
                if( ncplx.eq.2 .and. i.ne.j )then
                  do  k=1,nk
                    vloc0(k,j,i) = vloc0(k,j,i) + vatij*sinkr(k,ilat)
                  enddo
                  do  k=1,nk
                    vxc0(k,j,i) = vxc0(k,j,i) + vxcij*sinkr(k,ilat)
                  enddo
                  do  k=1,nk
                    exc0(k,j,i) = exc0(k,j,i) + excij*sinkr(k,ilat)
                  enddo
                endif
c
c               Close angular loops:
  415         continue
  425       continue
c
  799       continue
c           Close Atom i shell loop:
            i = i0 + norbshi
  800     continue
c
c         Close Atom i lattice loop:
  900   continue
c
c       Close Atom i loop:
 1000 continue
c
 1799   continue
c       Close Atom j shell loop:
        j = j0 + norbshj
 1800 continue
c
c       Close Atom j loop:
c
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
        call TIMER( '  vlocmat compute' )
c
c       Each local processor has full Hamiltonian over its share
c       of atomic potentials.
c       Reduce (sum) these to master to get total Hams.
c
c       Use the vtmp() as buffer for all reduces.
c
c       Reduce full local potential matrix:
        call MPREDUC8( master, nmat, vloc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, vloc0,1 )
c
c       Reduce Vxc0 matrix:
        call MPREDUC8( master, nmat, vxc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, vxc0,1 )
c
c       Reduce Exc0 matrix:
        call MPREDUC8( master, nmat, exc0, vtmp, icomm )
        call DCOPY( nmat, vtmp,1, exc0,1 )
c
        call TIMER( '  vlocmat Reduce ' )
c
c       The final output matrices are only on the master processor
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
