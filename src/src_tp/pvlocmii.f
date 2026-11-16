c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PVLOCMII
c
c
      subroutine PVLOCMII( IWR, atm0engy,
     $ norb,nk, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vij,vxc,exc, vtmp )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all 2-center II:N terms.
c
c Written: Peter A. Schultz, 20-July-1988, for v2.23e
c
c Revision history:
c  11May07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  17Oct04-PAS/2.59: task-parallel
c   7Mar02-PAS/2.52: ang constants extracted
c  22May01-PAS/2.47: condense radial scratch space
c   4Mar00-PAS/2.43: memory bugfix; rename "vl" routines;
c                    "enii" to "vlocmii"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Sep98-PAS/2.27: compact temp orbital matrix store
c  27Aug98-PAS/2.25: purge (N**7) triplet equivalences
c  11Aug98-PAS/2.24: match 3-center routines again
c  20Jul98-PAS/2.23: I:N split from ENIJ; split basis
c---------------------------------------------------------------
c
c  Angular integrals are evaluated analytically
c  Radial integrals are evaluated numerically
c
c  The temp matrix stores only really need the upper triangle,
c  but the bookkeeping gets painful given the ji loop structure
c  in this routine and the ij matrix structure.  While perhaps
c  wasteful, the current square atom block scheme is cheap and easy.
c
ctp: vlocmii task-parallelism is different - here we will parallelize
c  over atoms rather than over nuclear potentials, to reduce the
c  amount of communication needed at the end to assemble final matrix
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c output arrays:
      DIMENSION  vij(norb,norb,nk), vxc(norb,norb,nk),exc(norb,norb,nk)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $  nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  znuc(ntyp), occsh(nshld,ntyp),occij(nshld,nshld,ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  vatrad(nrd,ntyp)
      DIMENSION  vesrad(nrd,ntyp),vxcrad(nrd,ntyp),excrad(nrd,ntyp)
c scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vtmp(3,*)
c
c local declarations:
c
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      LOGICAL    keepsh
      DIMENSION  rnuc(3),ri(3), p(3), ci(10,9)
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
      emaxii = cutii
      emax = emaxii
c
      fourpi = pi + pi + pi + pi
c
ctp:start set up task parallel grid0xc
c  Get local MP info:
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c  nprocl: number of processors in this local image
      nprocl = nprocs
c  iprocl0: master/root node of set of processors in this image
      iprocl0 = master
      iprocl = iproc - iprocl0
c
c Set up task-parallel decomposition over atoms.
c  Rather than doing MOD( iatm, nproc ), will do block-wise
c  distribution to simplify fusion of final matrix.
      natmpr = natm / nprocl
      nprocbig = natm - nprocl*natmpr
      if( iprocl .lt. nprocbig ) natmpr = natmpr + 1
      iatm1 = 1 + iprocl*natmpr
      if( iprocl .ge. nprocbig ) iatm1 = iatm1 + nprocbig
      iatmn = iatm1 + natmpr - 1
      iv1locl = 1
      nvlocl = 0
c
ctp:end
      atm0enii = zero
c  Count up temp storage need (atom block diagonal):
      nvtmp = 0
      do  iatm=1,natm
        ityp = itypa(iatm)
        norbat = norba(ityp)
        nvtmp = nvtmp + norbat**2
ctp:    Count up how much will sit on this local proc.
        if( iatm.ge.iatm1 .and. iatm.le.iatmn )
     $   nvlocl = nvlocl + norbat**2
      enddo
ctp:
      nvlocl = 3*nvlocl
c
      nvtmp = 3*nvtmp
      call MKZERO( nvtmp, vtmp )
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
      ijnext = 0
      inext = 0
ctp      do 2000 iatm=1,natm
      do 2000 iatm=iatm1,iatmn
        ityp = itypa(iatm)
        nshli = numshl(ityp)
        istart = inext
        inext = istart + norba(ityp)
        norbat = norba(ityp)
        ij00 = ijnext
        ijnext = ij00 + norbat**2
c
        do  ir=1,3
          ri(ir) = ratm(ir,iatm) - rnuc(ir)
        enddo
        risquar = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
c
        alamini = alamin(ityp)
        rinmax = emax*(alamini+alamini+alaminn) /
     $              ( (alamini+alamini)*alaminn )
        almina = alamini
c
c  Test overlap of atom i and short-range pot due to "nuc"
        if( risquar.gt.rinmax ) goto 2000
c
        edelij = zero
        engyzz = zero
c
        if( klat.ne.1 .or. iatm.ne.nuc )then
c  Add nuc-iatm core repulsion to delta-energy.  Note: skipped IF
c  no nuc:iatm overlap (electronic and nuclear contributions cancel).
          zi = znuc(ityp)
          do 28 ishl=1,nshli
c
            norbshi = 2*lshel(ishl,ityp) + 1
            almini = ala(1,ishl,ityp)
            do  ial=2,nala(ishl,ityp)
              ali = ala(ial,ishl,ityp)
              if( ali.lt.almini ) almini = ali
            enddo
c
            do 26 jshl=1,ishl
c
c  Note that ishl and jshl must have same angular momentum
c  and that z contribution is in d1.le.d2 part
              occshl = norbshi*occij(jshl,ishl,ityp)
              if( occshl.eq.zero ) goto 26
c
              alminj = ala(1,jshl,ityp)
              do  jal=2,nala(jshl,ityp)
                alj = ala(jal,jshl,ityp)
                if( alj.lt.alminj ) alminj = alj
              enddo
              alsum = almini + alminj
              rismax = emax*(alsum+alaminn) / (alsum*alaminn)
              if( risquar.gt.rismax )  zi = zi - occshl
   26       continue
c
   28     continue
          zztot = znuc(ktyp)*zi
c
          engyadd = zztot / SQRT( risquar )
          atm0enii = atm0enii + engyadd
          engyzz = engyadd
          if( lstout.ge.2 ) write(IWR,9002) 'nuc,klat,iatm,E(Z-Z)=',
     $                                       nuc,klat,iatm,engyzz
        endif
c
        call CINIT( ci )
        call COFR( ci,ri )
c
        iorbat00 = 0
        i = istart
        do 1800 ishl=1,nshli
          nali = nala(ishl,ityp)
          lshli = lshel(ishl,ityp)
          li0 = lshli*lshli
          norbshi = 2*lshli + 1
          i0 = i
          iorbat0 = iorbat00
          iorbat00 = iorbat0 + norbshi
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
          alsum = almini + almina
          rinmax = emax*(alsum+alaminn) / (alsum*alaminn )
c
c  Test whether this radial function will contribute:
          if( risquar.gt.rinmax ) goto 1800
c
          jatm = iatm
          jtyp = ityp
          jstart = istart
          nshlj = ishl
c
          ij0 = ij00 + iorbat0
          j = jstart
          do 800 jshl=1,nshlj
            nalj = nala(jshl,jtyp)
            lshlj = lshel(jshl,jtyp)
            lj0 = lshlj*lshlj
            norbshj = 2*lshlj + 1
            j0 = j
            ijstart = ij0
            ij0 = ijstart + norbat*norbshj
c
            norbshj0 = norbshj
            joff = 0
            if( j0.eq.i0 )then
              norbshj0 = 0
              joff = 1
            endif
c
            alminj = ala(1,jshl,jtyp)
            do  jal=2,nalj
              alj = ala(jal,jshl,jtyp)
              if( alj.lt.alminj ) alminj = alj
            enddo
c
            rijnmax = emax*(alamini+alminj+alaminn) /
     $                   ( (alamini+alminj)*alaminn )
            if( risquar.gt.rijnmax ) goto 800
c
            keepsh = .false.
            occ = occij(ishl,jshl,ityp)
            if( occij(ishl,jshl,ityp).ne.zero )then
              keepsh = .true.
c             Double off-diags - only doing one side of triangle
              if( ishl.ne.jshl ) occ = occ + occ
            endif
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
              do 310 jal=1,nalj
                alj = ala(jal,jshl,jtyp)
c
c  Test whether this pair of alphas contributes:
                rijnmx = emax*(ali+alj+alaminn) /
     $                      ( (ali+alj)*alaminn )
                if( risquar.gt.rijnmx .and. .not. keepsh ) goto 310
c
                alsum = ali + alj
                exparg = ali*risquar + alj*risquar
                do  ir=1,3
                  p(ir) = two*(ali*ri(ir)+alj*ri(ir))
                enddo
                pmag = SQRT( p(1)*p(1)+p(2)*p(2)+p(3)*p(3) )
c
                call MKZERO( 5, vatint )
                call MKZERO( 5, vesint )
                call MKZERO( 5, vxcint )
                call MKZERO( 5, excint )
c
                nonzer = 0
                argtst = exparg - 0.25*pmag*pmag/alsum
c
                  if( argtst.gt.argmx ) goto 310
                  nonzer = 1
                  expfac = zero
                  if( exparg.lt.argmx ) expfac =  EXP( -exparg )
                  call VLRADZ( vesint,lmx,alsum,pmag,znuc(ktyp),
     $                         exparg,expfac )
c
c
c  "vlradz" puts the matrix element of 2.*znuc/r into vesint() for
c   self-energy correction term to vatint(), "vlrad" adds the
c   contributions of the xc-potential associated with atom nuc, as
c   well as the difference between the actual electrostatic
c   potential due to the electrons of atom nuc and 2.*abs(znuc)/r
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
     $               p,li,lj,ci(1,li),ci(1,lj) )
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
            ijs = ijstart
            norbshj = norbshj0
            i = i0
            do 425 li=li0+1,li0+norbshi
              i = i + 1
              lis = li - li0
              ijs = ijs + 1
              ij = ijs
c
              norbshj = norbshj + joff
              j = j0
              do 415 lj=lj0+1,lj0+norbshj
                j = j + 1
                ljs = lj - lj0
c
                vesij = vessh(lis,ljs)
                vtmp(1,ij) = vtmp(1,ij) + vijsh(lis,ljs)
                vtmp(2,ij) = vtmp(2,ij) + vxcsh(lis,ljs)
                vtmp(3,ij) = vtmp(3,ij) + excsh(lis,ljs)
c
                ij = ij + norbat
c
                if( keepsh .and. li.eq.lj )then
c subtract from delta-E contribution of spherical charge times itself
                  engyadd = -half*occ*vesij
                  atm0enii = atm0enii + engyadd
                  edelij = edelij + engyadd
                  if( lstout.ge.3 ) write(IWR,9003)
     $             '    li,lj,ial,jal,E(Z-e)=',li,lj,ial,jal,engyadd
                endif
c
  415         continue
  425       continue
c
c    Close atom j shell loop:
c
  800     j = j0 + norbshj
c
c    Close atom i shell loop:
c
 1800   i = i0 + norbshi
c
      if( lstout.gt.0 .and. (engyzz.ne.zero.or.edelij.ne.zero) )
     $  write(IWR,9001) nuc,klat,iatm,engyzz,edelij,engyzz+edelij
c
c    Close atom i loop:
c
 2000 continue
c
c    Close nucleus lattice loop:
c
 2900 continue
c
c    Close nuclear potential loop:
c
 3000 continue
ctp:start
c
      if( nprocl .gt. 1 )then
c       Merge the results onto the master/root node
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $     nvtmp,iv1locl,nvlocl, vtmp, icomm )

c        ... and broadcast it to everyone
        call MPBCAST8( iprocl0, nvtmp, vtmp, icomm )
c
c       Also hack together a sum for the atm0 term:
        nreduc = 1
        call MPREDUC8( iprocl0, nreduc, atm0enii, eiitmp, icomm )
C        call MPI_Reduce( atm0enii, eiitmp, nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   iprocl0, icomm, ierr )
        atm0enii = eiitmp
c
        call MPBCAST8( iprocl0, nreduc, atm0enii, icomm )
      endif
ctp:end
c
c Pack up contributions to Hams from matrix elements of potential nuc
c
      atm0engy = atm0engy + atm0enii
      do  k=1,nk
        i0 = 0
        ij0 = 0
        do  iatm=1,natm
          ityp = itypa(iatm)
          norbat = norba(ityp)
c
          j = i0
          ij = ij0
          do  jorbat=1,norbat
            i = j
            j = j + 1
            do  iorbat=jorbat,norbat
              i = i + 1
              ij = ij + 1
              vij(i,j,k) = vij(i,j,k) + vtmp(1,ij)
              vxc(i,j,k) = vxc(i,j,k) + vtmp(2,ij)
              exc(i,j,k) = exc(i,j,k) + vtmp(3,ij)
            enddo
            ij = ij + jorbat
          enddo
c
          i0 = i0 + norbat
          ij0 = ij0 + norbat**2
c
        enddo
c
      enddo
c
c    That's all Folks!
c
      RETURN
      END
