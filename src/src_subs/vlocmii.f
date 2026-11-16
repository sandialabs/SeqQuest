c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCMII
c
c
      subroutine VLOCMII( IWR, atm0engy,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all 2-center II:N terms.
c
c Written: Peter A. Schultz, 20-July-1988, for v2.23e
c
c Revision history:
c  17Oct08-PAS/2.62: stripe-distributed-parallel
c  20Jun07-PAS/2.60: merged tp and serial
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
c  vlocmiii/tp task-parallel is different - here we parallelize
c  over atoms rather than over nuclear potentials, to reduce the
c  amount of communication needed at the end to assemble final matrix
c
cpas: Coding caution: occij() is not symmetric!!!
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  vloc0(*)
      DIMENSION  vxc0(*), exc0(*)
c Input arrays:
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
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vtmp(3,*)
c
c Local declarations:
c
      LOGICAL    Lijkmats,Lstrmats
      LOGICAL    keepsh, diagshl
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      DIMENSION  rnuc(3),rj(3), p(3), cj(10,9)
      DIMENSION  vatint(5),vesint(5),vxcint(5),excint(5)
c
      DIMENSION  vnorm(9)
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/VLOCMII: stripe-parallel'
      atm0enii = zero
      Lijkmats = .false.
      Lstrmats = .true.
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxii = cutii
      emax = emaxii
c
      fourpi = pi + pi + pi + pi
c
c  Get local MP info:
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
c  Count up temp storage need (atom block diagonal):
      iv1locl = 1
      nvlocl = 0
      nvtmp = 0
      do  jatm=1,natm
        jtyp = itypa(jatm)
        norbat = norba(jtyp)
        nvtmp = nvtmp + norbat**2
c       Count up how much will sit on this local proc.
        if( jatm.ge.myatm1 .and. jatm.le.myatm2 )
     $   nvlocl = nvlocl + norbat**2
      enddo
c
c  This is how much sits on this local processor
      nvlocl = 3*nvlocl
      call MKZERO( nvlocl, vtmp )
c
      if( lstout.ne.0 )  write(IWR,9000)
 9000 format(/3x,'nuc  klat jatm',9x,'E(Z-Z)',12x,'E(E-Z)',12x,'TOT')
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
c     Atom j loop:
c
      ijnext = 0
      jnext = 0
ctp      do 2000 jatm=1,natm
      do 2000 jatm=myatm1,myatm2
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
        ij00 = ijnext
        ijnext = ij00 + norbatj**2
c
        do  ir=1,3
          rj(ir) = ratm(ir,jatm) - rnuc(ir)
        enddo
        rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
        alaminj = alamin(jtyp)
        rjnmax = emax*(alaminj+alaminj+alaminn) /
     $              ( (alaminj+alaminj)*alaminn )
        almina = alaminj
c
c       Test overlap of Atom j and short-range pot due to "nuc"
        if( rjsquar.gt.rjnmax ) goto 2000
c
        edelij = zero
        engyzz = zero
c
        if( klat.ne.1 .or. jatm.ne.nuc )then
c         Add nuc-jatm core repulsion to delta-energy.
c         Note: skipped IF no nuc:jatm overlap
c         (electronic and nuclear contributions cancel).
          zj = znuc(jtyp)
          do 28 jshl=1,nshlj
c
            norbshj = 2*lshel(jshl,jtyp) + 1
            alminj = ala(1,jshl,jtyp)
            do  jal=2,nala(jshl,jtyp)
              alj = ala(jal,jshl,jtyp)
              if( alj.lt.alminj ) alminj = alj
            enddo
c
            do 26 ishl=1,jshl
c
c             NOTE: (1) ishl and jshl must have same angular momentum
c              and (2) z contribution is in occij(is,js.le.is,ityp)
              occshl = norbshj*occij(ishl,jshl,jtyp)
              if( occshl.eq.zero ) goto 26
c
              almini = ala(1,ishl,jtyp)
              do  ial=2,nala(ishl,jtyp)
                ali = ala(ial,ishl,jtyp)
                if( ali.lt.almini ) almini = ali
              enddo
              alsum = almini + alminj
              rjsmax = emax*(alsum+alaminn) / (alsum*alaminn)
              if( rjsquar.gt.rjsmax )  zj = zj - occshl
   26       continue
c
   28     continue
          zztot = znuc(ktyp)*zj
c
          engyadd = zztot / SQRT( rjsquar )
          atm0enii = atm0enii + engyadd
          engyzz = engyadd
          if( lstout.ge.2 ) write(IWR,9002) 'nuc,klat,jatm,E(Z-Z)=',
     $                                       nuc,klat,jatm,engyzz
        endif
c
        call CINIT( cj )
        call COFR( cj,rj )
c
c       Atom shell loop:
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
          alsum = alminj + almina
          rjnmax = emax*(alsum+alaminn) / (alsum*alaminn )
c
c         Test whether this radial function will contribute:
          if( rjsquar.gt.rjnmax ) goto 1799
c
          iatm = jatm
          ityp = jtyp
          istart = jstart
          ishl1 = jshl
          do  ishl=1,ishl1-1
            lshli = lshel(ishl,ityp)
            norbshi = 2*lshli + 1
            istart = istart + norbshi
          enddo
c
c         Atom shell loop-2:
c
          iorbat0 = jorbat0
          ij0 = ij00 + jorbat0*norbatj + iorbat0
          i = istart
          do 800 ishl=ishl1,nshlj
            nali = nala(ishl,ityp)
            lshli = lshel(ishl,ityp)
            norbshi = 2*lshli + 1
            li0 = lshli*lshli
c
            i0 = i
            ijstart = ij0
            ij0 = ijstart + norbshi
c
            diagshl = .false.
            if( i0.eq.j0 ) diagshl = .true.
c
            almini = ala(1,ishl,ityp)
            do  ial=2,nali
              ali = ala(ial,ishl,ityp)
              if( ali.lt.almini ) almini = ali
            enddo
c
            rijnmax = emax*(alaminj+almini+alaminn) /
     $                   ( (alaminj+almini)*alaminn )
            if( rjsquar.gt.rijnmax ) goto 799
c
            keepsh = .false.
c           NOTE: z contribution is in occij(is,js.le.is,ityp)
            occ = occij(ishl,jshl,jtyp)
            if( occij(ishl,jshl,jtyp).ne.zero )then
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
c           Alpha loops:
c
            enijsum = zero
            coulsum = zero
            do 320 jal=1,nalj
              alj = ala(jal,jshl,jtyp)
c
              do 310 ial=1,nali
                ali = ala(ial,ishl,ityp)
c
c               Test whether this pair of alphas contributes:
                rijnmx = emax*(alj+ali+alaminn) /
     $                      ( (alj+ali)*alaminn )
                if( rjsquar.gt.rijnmx .and. .not. keepsh ) goto 310
c
                alsum = ali + alj
                exparg = ali*rjsquar + alj*rjsquar
                do  ir=1,3
                  p(ir) = two*(ali*rj(ir)+alj*rj(ir))
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
c               "vlradz" puts the matrix element of 2.*znuc/r into vesint() for
c                self-energy correction term to vatint(), "vlrad" adds the
c                contributions of the xc-potential associated with atom nuc,
c                as well as the difference between the actual electrostatic
c                potential due to the electrons of atom nuc and 2.*abs(znuc)/r
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
     $               p,lj,li,cj(1,lj),cj(1,li) )
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
            ijs = ijstart
            mli1 = 1
            do 425 mlj=1,norbshj
              lj = lj0 + mlj
              j = j0 + mlj
c
              ijs0 = ijs
              ijs = ijs + norbatj
              if( diagshl ) mli1 = mlj
              do 415 mli=mli1,norbshi
                li = li0 + mli
                i = i0 + mli
                ij = ijs0 + mli
c
                vesij = vessh(mli,mlj)
                vtmp(1,ij) = vtmp(1,ij) + vijsh(mli,mlj)
                vtmp(2,ij) = vtmp(2,ij) + vxcsh(mli,mlj)
                vtmp(3,ij) = vtmp(3,ij) + excsh(mli,mlj)
c
                if( keepsh .and. li.eq.lj )then
c                 Subtract from delta-E: spherical charge times itself
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
  799       continue
c           Close Atom shell loop-2:
            i = i0 + norbshi
  800     continue
c
 1799     continue
c         Close Atom shell loop-1:
          j = j0 + norbshj
 1800   continue
c
        if( lstout.gt.0 .and. (engyzz.ne.zero.or.edelij.ne.zero) )
     $    write(IWR,9001) nuc,klat,jatm,engyzz,edelij,engyzz+edelij
c
c       Close Atom j loop:
 2000 continue
c
c       Close Nucleus lattice loop:
c
 2900 continue
c
c       Close Nucleus potential loop:
c
 3000 continue
c
      if( nprocs .gt. 1 ) call TIMER( '  vlocmii compute' )
c
      if( nprocs .gt. 1 )then
c       Assemble the sum of atm0 terms:
        nreduc = 1
        call MPREDUC8( master, nreduc, atm0enii, eiitmp, icomm )
        atm0enii = eiitmp
        call MPBCAST8( master, nreduc, atm0enii, icomm )
      endif
      atm0engy = atm0engy + atm0enii
c
      if( Lijkmats )then
c       Move results into standard H(i,j,k) storage ...
        if( iproc.lt.nprocstr )then
c         Merge the results onto the master/root node
          call MPMERGER8( nprocstr, master, iproc,
     $       nvtmp,iv1locl,nvlocl, vtmp, icomm )
c          ... and broadcast it to everyone
          call MPBCAST8( master, nvtmp, vtmp, icomm )
        endif
c       Put contributions to Ham into ijk matrix:
        call IIPUTIJK( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   vloc0,vxc0,exc0, vtmp )
c
      elseif( Lstrmats )then
c       Move results into distributed stripe storage ...
        if( iproc.lt.nprocstr )
     $  call IIPUTSTR( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   myj0,myatm1,myatm2,
     $   vloc0,vxc0,exc0, vtmp )
c
      else
c       No other alternatives right now ...
        call STOPXERR( 'invalid vlocmii merge - code error' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IIPUTIJK
c
c
      subroutine IIPUTIJK( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: put diagonal atom block vlocmii into vloc(i,j,k) matrices
c Written: Peter A. Schultz, 8-October-2008, for 2.62
c Revision history:
c   8Oct08-PAS/2.62: pushed down from vlocmii
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output arrays:
      DIMENSION  vloc0(norb,norb,nk)
      DIMENSION  vxc0(norb,norb,nk),exc0(norb,norb,nk)
c Input arrays:
      DIMENSION  vtmp(3,*)
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
c
        j0 = 0
        ijnext = 0
        do 1000 jatm=1,natm
          jtyp = itypa(jatm)
          norbatj = norba(jtyp)
c
          j = j0
          ij = ijnext
          ijnext = ijnext + norbatj**2
c
          do 420 jorbat=1,norbatj
            i = j
            j = j + 1
            do 410 iorbat=jorbat,norbatj
              i = i + 1
              ij = ij + 1
c
              vloc0(i,j,k) = vloc0(i,j,k) + vtmp(1,ij)
              vxc0(i,j,k) = vxc0(i,j,k) + vtmp(2,ij)
              exc0(i,j,k) = exc0(i,j,k) + vtmp(3,ij)
c
  410       continue
            ij = ij + jorbat
  420     continue
c
          j0 = j0 + norbatj
c
 1000   continue
c
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IIPUTKIJ
c
c
      subroutine IIPUTKIJ( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: put diagonal atom block vlocmii into vloc(k,i,j) matrices
c Written: Peter A. Schultz, 8-October-2008, for 2.62
c Revision history:
c   8Oct08-PAS/2.62: pushed down from vlocmii
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output arrays:
      DIMENSION  vloc0(nk,norb,norb)
      DIMENSION  vxc0(nk,norb,norb),exc0(nk,norb,norb)
c Input arrays:
      DIMENSION  vtmp(3,*)
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c >>>> EXECUTABLE CODE:
c
        j0 = 0
        ijnext = 0
        do 1000 jatm=1,natm
          jtyp = itypa(jatm)
          norbatj = norba(jtyp)
c
          j = j0
          ij = ijnext
          ijnext = ijnext + norbatj**2
c
          do 420 jorbat=1,norbatj
            i = j
            j = j + 1
            do 410 iorbat=jorbat,norbatj
              i = i + 1
              ij = ij + 1
c
              do  k=1,nk
                vloc0(k,i,j) = vloc0(k,i,j) + vtmp(1,ij)
                vxc0(k,i,j) = vxc0(k,i,j) + vtmp(2,ij)
                exc0(k,i,j) = exc0(k,i,j) + vtmp(3,ij)
              enddo
c
  410       continue
            ij = ij + jorbat
  420     continue
c
          j0 = j0 + norbatj
c
 1000   continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IIPUTSTR
c
c
      subroutine IIPUTSTR( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ myj0,myatm1,myatm2,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: put diagonal atom block vlocmii into vloc(stripe) arrays
c Written: Peter A. Schultz, 8-October-2008, for 2.62
c Revision history:
c   none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output arrays:
      DIMENSION  vloc0(nk,*)
      DIMENSION  vxc0(nk,*), exc0(nk,*)
c Input arrays:
      DIMENSION  vtmp(3,*)
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c >>>> EXECUTABLE CODE:
c
        iivlnext = 0
        jjnext = 0
        jnext = myj0
        do 1000 jatm=myatm1,myatm2
          jtyp = itypa(jatm)
          norbatj = norba(jtyp)
          jstart = jnext
          jnext = jstart + norbatj
c
          njstripe = norb - jstart
          matjstr = norbatj*njstripe
          jjstart = jjnext
          jjnext = jjstart + matjstr*ncplx
c
          j0 = jstart
          i0 = jstart
c
          ij0 = jjstart
          ij0 = ij0 - norbatj
c
          iivl = iivlnext
          iivlnext = iivlnext + norbatj**2
c
          do 420 jorbat=1,norbatj
            j = j0 + jorbat
c
            ij0 = ij0 + norbatj
            do 410 iorbat=jorbat,norbatj
              i = i0 + iorbat
              ij = ij0 + iorbat
              iivl = iivl + 1
c
              do  k=1,nk
                vloc0(k,ij) = vloc0(k,ij) + vtmp(1,iivl)
                vxc0(k,ij) = vxc0(k,ij) + vtmp(2,iivl)
                exc0(k,ij) = exc0(k,ij) + vtmp(3,iivl)
              enddo
c
  410       continue
            iivl = iivl + jorbat
  420     continue
c
 1000   continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCMII1
c
c
      subroutine VLOCMII1( IWR, atm0engy,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat  ,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ vloc0,vxc0,exc0, vtmp )
c---------------------------------------------------------------
c Purpose: calculate electron-nuclear matrix elements for the
c          local part of the spherical reference potentials
c          for all 2-center II:N terms.
c
c Written: Peter A. Schultz, 20-July-1988, for v2.23e
c
c Revision history:
c   7Sep08-PAS/2.62: switch ij-ji to natural order
c  20Jun07-PAS/2.60: merged tp and serial
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
c  vlocmiii/tp task-parallel is different - here we parallelize
c  over atoms rather than over nuclear potentials, to reduce the
c  amount of communication needed at the end to assemble final matrix
c
cpas: Coding caution: occij() is not symmetric!!!
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
c Input arrays:
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
c Scratch arrays:
      DIMENSION  wkrad(nrd,*)
      DIMENSION  vtmp(3,*)
c
c Local declarations:
c
      LOGICAL    keepsh, diagshl
      DIMENSION  vessh(5,5),vijsh(5,5),vxcsh(5,5),excsh(5,5)
      DIMENSION  rnuc(3),rj(3), p(3), cj(10,9)
      DIMENSION  vatint(5),vesint(5),vxcint(5),excint(5)
c
      DIMENSION  vnorm(9)
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/VLOCMII1: ijkmat'
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxii = cutii
      emax = emaxii
c
      fourpi = pi + pi + pi + pi
c
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
c  Rather than doing MOD( jatm, nproc ), will do block-wise
c  distribution to simplify fusion of final matrix.
      natmpr = natm / nprocl
      nprocbig = natm - nprocl*natmpr
      nprocsml = nprocl - nprocbig
c     Here, I assume that iprocl is sequential beginning with 0.
      if( iprocl .ge. nprocsml ) natmpr = natmpr + 1
      jatm1 = 1 + iprocl*natmpr
      if( iprocl .ge. nprocsml ) jatm1 = jatm1 - nprocsml
      jatmn = jatm1 + natmpr - 1
      iv1locl = 1
      nvlocl = 0
c
      atm0enii = zero
c  Count up temp storage need (atom block diagonal):
      nvtmp = 0
      do  jatm=1,natm
        jtyp = itypa(jatm)
        norbat = norba(jtyp)
        nvtmp = nvtmp + norbat**2
c       Count up how much will sit on this local proc.
        if( jatm.ge.jatm1 .and. jatm.le.jatmn )
     $   nvlocl = nvlocl + norbat**2
      enddo
c
c  This is how much sits on this local processor
      nvlocl = 3*nvlocl
c
      nvtmp = 3*nvtmp
      call MKZERO( nvtmp, vtmp )
c
      if( lstout.ne.0 )  write(IWR,9000)
 9000 format(/3x,'nuc  klat jatm',9x,'E(Z-Z)',12x,'E(E-Z)',12x,'TOT')
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
c     Atom j loop:
c
      ijnext = 0
      jnext = 0
ctp      do 2000 jatm=1,natm
      do 2000 jatm=jatm1,jatmn
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
        ij00 = ijnext
        ijnext = ij00 + norbatj**2
c
        do  ir=1,3
          rj(ir) = ratm(ir,jatm) - rnuc(ir)
        enddo
        rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
        alaminj = alamin(jtyp)
        rjnmax = emax*(alaminj+alaminj+alaminn) /
     $              ( (alaminj+alaminj)*alaminn )
        almina = alaminj
c
c       Test overlap of Atom j and short-range pot due to "nuc"
        if( rjsquar.gt.rjnmax ) goto 2000
c
        edelij = zero
        engyzz = zero
c
        if( klat.ne.1 .or. jatm.ne.nuc )then
c         Add nuc-jatm core repulsion to delta-energy.
c         Note: skipped IF no nuc:jatm overlap
c         (electronic and nuclear contributions cancel).
          zj = znuc(jtyp)
          do 28 jshl=1,nshlj
c
            norbshj = 2*lshel(jshl,jtyp) + 1
            alminj = ala(1,jshl,jtyp)
            do  jal=2,nala(jshl,jtyp)
              alj = ala(jal,jshl,jtyp)
              if( alj.lt.alminj ) alminj = alj
            enddo
c
            do 26 ishl=1,jshl
c
c             NOTE: (1) ishl and jshl must have same angular momentum
c              and (2) z contribution is in occij(is,js.le.is,ityp)
              occshl = norbshj*occij(ishl,jshl,jtyp)
              if( occshl.eq.zero ) goto 26
c
              almini = ala(1,ishl,jtyp)
              do  ial=2,nala(ishl,jtyp)
                ali = ala(ial,ishl,jtyp)
                if( ali.lt.almini ) almini = ali
              enddo
              alsum = almini + alminj
              rjsmax = emax*(alsum+alaminn) / (alsum*alaminn)
              if( rjsquar.gt.rjsmax )  zj = zj - occshl
   26       continue
c
   28     continue
          zztot = znuc(ktyp)*zj
c
          engyadd = zztot / SQRT( rjsquar )
          atm0enii = atm0enii + engyadd
          engyzz = engyadd
          if( lstout.ge.2 ) write(IWR,9002) 'nuc,klat,jatm,E(Z-Z)=',
     $                                       nuc,klat,jatm,engyzz
        endif
c
        call CINIT( cj )
        call COFR( cj,rj )
c
c       Atom shell loop:
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
          alsum = alminj + almina
          rjnmax = emax*(alsum+alaminn) / (alsum*alaminn )
c
c         Test whether this radial function will contribute:
          if( rjsquar.gt.rjnmax ) goto 1799
c
          iatm = jatm
          ityp = jtyp
          istart = jstart
          ishl1 = jshl
          do  ishl=1,ishl1-1
            lshli = lshel(ishl,ityp)
            norbshi = 2*lshli + 1
            istart = istart + norbshi
          enddo
c
c         Atom shell loop-2:
c
          iorbat0 = jorbat0
          ij0 = ij00 + jorbat0*norbatj + iorbat0
          i = istart
          do 800 ishl=ishl1,nshlj
            nali = nala(ishl,ityp)
            lshli = lshel(ishl,ityp)
            norbshi = 2*lshli + 1
            li0 = lshli*lshli
c
            i0 = i
            ijstart = ij0
            ij0 = ijstart + norbshi
c
            diagshl = .false.
            if( i0.eq.j0 ) diagshl = .true.
c
            almini = ala(1,ishl,ityp)
            do  ial=2,nali
              ali = ala(ial,ishl,ityp)
              if( ali.lt.almini ) almini = ali
            enddo
c
            rijnmax = emax*(alaminj+almini+alaminn) /
     $                   ( (alaminj+almini)*alaminn )
            if( rjsquar.gt.rijnmax ) goto 799
c
            keepsh = .false.
c           NOTE: z contribution is in occij(is,js.le.is,ityp)
            occ = occij(ishl,jshl,jtyp)
            if( occij(ishl,jshl,jtyp).ne.zero )then
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
c           Alpha loops:
c
            enijsum = zero
            coulsum = zero
            do 320 jal=1,nalj
              alj = ala(jal,jshl,jtyp)
c
              do 310 ial=1,nali
                ali = ala(ial,ishl,ityp)
c
c               Test whether this pair of alphas contributes:
                rijnmx = emax*(alj+ali+alaminn) /
     $                      ( (alj+ali)*alaminn )
                if( rjsquar.gt.rijnmx .and. .not. keepsh ) goto 310
c
                alsum = ali + alj
                exparg = ali*rjsquar + alj*rjsquar
                do  ir=1,3
                  p(ir) = two*(ali*rj(ir)+alj*rj(ir))
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
c               "vlradz" puts the matrix element of 2.*znuc/r into vesint() for
c                self-energy correction term to vatint(), "vlrad" adds the
c                contributions of the xc-potential associated with atom nuc,
c                as well as the difference between the actual electrostatic
c                potential due to the electrons of atom nuc and 2.*abs(znuc)/r
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
     $               p,lj,li,cj(1,lj),cj(1,li) )
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
            ijs = ijstart
            mli1 = 1
            do 425 mlj=1,norbshj
              lj = lj0 + mlj
              j = j0 + mlj
c
              ijs0 = ijs
              ijs = ijs + norbatj
              if( diagshl ) mli1 = mlj
              do 415 mli=mli1,norbshi
                li = li0 + mli
                i = i0 + mli
                ij = ijs0 + mli
c
                vesij = vessh(mli,mlj)
                vtmp(1,ij) = vtmp(1,ij) + vijsh(mli,mlj)
                vtmp(2,ij) = vtmp(2,ij) + vxcsh(mli,mlj)
                vtmp(3,ij) = vtmp(3,ij) + excsh(mli,mlj)
c
                if( keepsh .and. li.eq.lj )then
c                 Subtract from delta-E: spherical charge times itself
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
  799       continue
c           Close Atom shell loop-2:
            i = i0 + norbshi
  800     continue
c
 1799     continue
c         Close Atom shell loop-1:
          j = j0 + norbshj
 1800   continue
c
        if( lstout.gt.0 .and. (engyzz.ne.zero.or.edelij.ne.zero) )
     $    write(IWR,9001) nuc,klat,jatm,engyzz,edelij,engyzz+edelij
c
c       Close Atom j loop:
 2000 continue
c
c       Close Nucleus lattice loop:
c
 2900 continue
c
c       Close Nucleus potential loop:
c
 3000 continue
c
      if( nprocl .gt. 1 )then
c       Merge the results onto the master/root node
        call MPMERGER8( nprocl, iprocl0, iprocl,
     $     nvtmp,iv1locl,nvlocl, vtmp, icomm )

c        ... and broadcast it to everyone
        call MPBCAST8( iprocl0, nvtmp, vtmp, icomm )
c
c       Assemble the sum of atm0 terms:
        nreduc = 1
        call MPREDUC8( iprocl0, nreduc, atm0enii, eiitmp, icomm )
        atm0enii = eiitmp
c
        call MPBCAST8( iprocl0, nreduc, atm0enii, icomm )
      endif
c
c Pack up contributions to Hams from matrix elements of potential nuc
c
      atm0engy = atm0engy + atm0enii
c
        i0 = 0
        ijnext = 0
        do  iatm=1,natm
          ityp = itypa(iatm)
          norbat = norba(ityp)
c
          j = i0
          ij = ijnext
          do  jorbat=1,norbat
            i = j
            j = j + 1
            do  iorbat=jorbat,norbat
              i = i + 1
              ij = ij + 1
              do  k=1,nk
                vloc0(k,i,j) = vloc0(k,i,j) + vtmp(1,ij)
                vxc0(k,i,j) = vxc0(k,i,j) + vtmp(2,ij)
                exc0(k,i,j) = exc0(k,i,j) + vtmp(3,ij)
              enddo
            enddo
            ij = ij + jorbat
          enddo
c
          i0 = i0 + norbat
          ijnext = ijnext + norbat**2
c
        enddo
c
c    That's all Folks!
c
      RETURN
      END
