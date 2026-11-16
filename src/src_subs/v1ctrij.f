c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> V1CTRIJ
c
c
      subroutine V1CTRIJ( idmatfl,ilocfl, conv,
     $ norb,nk, natm,ntyp,nshld,nald, noad, nlat,
c output quantities
     $ v,  xp,   dmatrd,   dcountr,rh1coef,
     $ exc1c, cv,ce,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,ioffa, occij,
     $ lmx1ctr,almnv1c, nrad,nrd, radmsh,radwt,
     $ ratm,rlat, vk, coskr,sinkr,
     $ neq3d,ieqn,ieqi,ieqj,ipermut,
     $ vradsr,v1csum,velem, pradsr,p1csum,pelem,xradsr,x1csum,xelem,
     $ besstor,radfac,w1,w2,w3,w4,w5,w6,
     $ ves1c, pot1c,xen1c,
     $ nsyms,dsym,naofsym,
     $ alxcmin,nalfxc,nalfxcd,alfxc,alsq,ipvt, nrxcfit,vxcrad,excrad,
     $ nang, angylm,angpts,angwts, nfitd,nfityp,radwf,
     $ rholoc,vxcloc,excloc, vlm,elm,rfac, vg,eg,
     $ orb,grorb,grgrorb, wf1cwk1,wf1cwk2,wf1cwk3,wf1cwk4,wf1cwk5,
     $ orbj,grorbj,grgrorj, gofac,ggofac,
     $ eswork, dmatwk )
c---------------------------------------------------------------
c Purpose: compute matrix elements off of local 1-center pots.
c
c Written: Peter J. Feibelman
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c   7Mar02-PAS/2.52: ang constants extracted
c  24Jul01-PAS/2.49: packed dmat
c   2Jul01-PAS/2.48: replace writes
c   1Mar00-PAS/2.43: bugfix for non-xcfit
c  21Dec99-PAS/2.42: cosmetic clean-up
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  29Jun98-PAS/2.23: split basis
c  28Apr97-PAS/2.17: naming conventions
c  22Feb93-PAS
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c      DATA ARRAYS:
c
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  occij(*), alamin(ntyp), norba(ntyp),ioffa(natm)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  ratm(3,natm),rlat(3,nlat)
c ... only used in blkvofn
      DIMENSION  vk(3,nk)
c ... only used in v1ctrij
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
c
c atom symmetry data
c ... only in v1ctrij:
      DIMENSION  ieqn(neq3d),ieqi(neq3d),ieqj(neq3d),ipermut(*)
c ... only in blkvofn:
      DIMENSION  dsym(5,5,3,*),naofsym(natm,*)
c
c ... used in blkvofn for xc-fit:
      DIMENSION  alxcmin(ntyp),nalfxc(ntyp),alfxc(nalfxcd,ntyp)
      DIMENSION  alsq(nalfxcd,nalfxcd,3,ntyp),ipvt(nalfxcd,3,ntyp)
      DIMENSION  nrxcfit(ntyp),vxcrad(nrd,ntyp),excrad(nrd,ntyp)
      DIMENSION  angylm(nang,9),angpts(3,nang),angwts(nang)
      DIMENSION  nfityp(ntyp),radwf(nang,nrd,noad,nfitd)
c ... from blkvofn into v1ctrij:
      DIMENSION  ves1c(nrd,10,natm)
c ... from blkvofn into v1ctrij, xc-fit:
      DIMENSION  pot1c(nrd,10,natm),xen1c(nrd,10,natm)
c
c---------------------------------------------------------------
c
c      OUTPUT ARRAYS:
c
c ... out of v1ctrij
      DIMENSION  v(norb,norb,nk),xp(norb,norb,nk)
      DIMENSION  dmatrd(2,norb,norb,nk)
complex dmatrd
c ... out of blkvofn
      DIMENSION  dcountr(10,natm),rh1coef(10,noad,natm)
c ... out of blkvofn, xc-fit
      DIMENSION  cv(10,nalfxcd,natm),ce(10,nalfxcd,natm)
c
c---------------------------------------------------------------
c
c      SCRATCH ARRAYS:
c
c ... only in v1ctrij
      DIMENSION  vradsr(7,10,neq3d),v1csum(neq3d,49),velem(neq3d)
      DIMENSION  besstor(nrd,7),radfac(nrd,7)
      DIMENSION  w1(nrd),w2(nrd),w3(nrd),w4(nrd),w5(nrd),w6(nrd)
c ... only in v1ctrij, for xc-fit:
      DIMENSION  pradsr(7,10,neq3d),p1csum(neq3d,49),pelem(neq3d)
      DIMENSION  xradsr(7,10,neq3d),x1csum(neq3d,49),xelem(neq3d)
c ... used only in blkvofn:
      DIMENSION  orb(nlat),grorb(3,nlat),grgrorb(3,3,nlat)
      DIMENSION  wf1cwk1(3,nlat),wf1cwk2(3,3,nlat),
     $           wf1cwk3(3,nlat),wf1cwk4(nlat),wf1cwk5(nlat)
      DIMENSION  eswork(nrd), dmatwk(noad,noad,nlat)
c ... only in blkvofn, for xc-fit:
      DIMENSION  rholoc(nrd,nang),vxcloc(nrd,nang),excloc(nrd,nang)
      DIMENSION  orbj(noad),grorbj(3,noad),grgrorj(3,3,noad)
      DIMENSION  gofac(nang),ggofac(nang)
      DIMENSION  vlm(nrd),elm(nrd),rfac(nrd,nalfxcd)
      DIMENSION  vg(nalfxcd,9),eg(nalfxcd,9)
c---------------------------------------------------------------
c
c      LOCAL DECLARATIONS:
c
      LOGICAL    pjfsrc
      DATA  pjfsrc / .true. /
      LOGICAL    ntest
      DIMENSION  rnuc(3),ri(3),rj(3)
      DIMENSION  p(3), ci(10,9),cj(10,9),cigrgr(10,3,3)
c
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c---------------------------------------------------------------
c
c      this program calculates the matrix elements of 1-center
c      local potentials of the form v(r)* an angular factor
c      which may be 1, x, y, z, xy, yz, zx, x**2, y**2, z**2,
c      depending on the value of lnuc, which varies from 1 through 10.
c
c      these local potentials are evaluated, given the delta
c      density matrix "dmat" as input. in order to guarantee
c      convergence of the lattice sums, the local potentials are
c      constructed so as to have no long-ranged part.
c
c      the radial wave functions are assumed to have been
c      normalized (this is done in routine orbtl, where
c      the cala coefficients have been rescaled). the normali-
c      zation of the angular wave functions is done in the present
c      routine via the division by the quantities "vnorm"
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c >>>> EXECUTABLE CODE:
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Get angular normalizations:
      call VNORMS9( vnorm )
      emax = - LOG( conv )
      argmx = half*cutexp
      fourpi = pi + pi + pi + pi
c
      call MKZERO( 10*noad*natm, rh1coef )
c
c    Read momentum space complex density matrix
c
      REWIND( unit=idmatfl )
cpas: temp work-around/expand new packed dmat into old local usage
      do  k=1,nk
        call READBIG( idmatfl, norb*norb, v )
        do  j=1,norb
          do  i=1,j-1
            dmatrd(1,i,j,k) = v(j,i,1)
            dmatrd(2,i,j,k) = v(i,j,1)
          enddo
            dmatrd(1,j,j,k) = v(j,j,1)
            dmatrd(2,j,j,k) = zero
          do  i=j+1,norb
            dmatrd(1,i,j,k) = v(i,j,1)
            dmatrd(2,i,j,k) = -v(j,i,1)
          enddo
        enddo
      enddo
      call MKZERO( norb*norb, v )
c
c     Check if have local density file
      if( ilocfl.ne.0 ) REWIND( unit=ilocfl )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c   Compute short-range local potentials:
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    Nucleus loop:
c
      do 5 nuc=1,natm
        ktyp = itypa(nuc)
        if( lmx1ctr(ktyp).lt.0 ) goto 5
c
c      Construct 1-center potentials from nucleus "nuc"
c
        nrxc = nrxcfit(ktyp)
        nrho = nrxc
        if( ilocfl.gt.0 )then
          read(ilocfl)  nrho,mang
          if( nrho.lt.0 .or. nrho.gt.nrd .or. mang.ne.nang )then
c           We have a problem - kill this one.
            call FLGETIWR( IWR )
            write(IWR,*)'ilocfl problem, nuc=',nuc,nrho,mang,nrd,nang
            read(ilocfl) dum
            ilocfl = -ilocfl
            nrho = nrxc
          endif
        endif
c
c      update of 1/90:
c
c      BLKVOFN now also determines the density in the neighborhood
c      of atom nuc on a spherical mesh, in order to compute the
c      exchange-correlation potential and energy density on that mesh.
c      The l,m components of these functions of position are
c      then projected out, for l=0,1 and 2. Then, the radial
c      functions that result from the projection are fit to
c      linear combinations of Gaussians, and the fitted quantities
c      are used to account for rapid variation in x-c
c      quantities more accurately than previously.
c
        itypxc = nfityp(ktyp)
        ktypxc = MAX( 1, itypxc )
        call BLKVOFN( ilocfl, nuc,ktyp, lmx1ctr(ktyp),almnv1c,
     $   norb,nk, natm,ntyp,nshld,nald, noad, nlat,
     $   itypa, numshl,lshel,nala,ala,cala,  norba,ioffa, occij,
     $   nrad(ktyp),nrd,radmsh(1,ktyp),radwt(1,ktyp),
     $   ratm,rlat,vk, nsyms,dsym,naofsym,
     $   dmatrd,
     $   dcountr(1,nuc),rh1coef(1,1,nuc), cv(1,1,nuc),ce(1,1,nuc),
     $    ves1c(1,1,nuc), pot1c(1,1,nuc),xen1c(1,1,nuc),
     $   alxcmin(ktyp),nalfxc(ktyp),nalfxcd,alfxc(1,ktyp),
     $    alsq(1,1,1,ktyp),ipvt(1,1,ktyp),
     $   nrxc, vxcrad(1,ktyp),excrad(1,ktyp),
     $   itypxc, nang, angylm,angpts,angwts, radwf(1,1,1,ktypxc),
     $   nrho, rholoc,vxcloc,excloc, vlm,elm,rfac, vg,eg,
     $   orb,grorb,grgrorb, wf1cwk1,wf1cwk2,wf1cwk3,wf1cwk4,wf1cwk5,
     $   orbj,grorbj,grgrorj, gofac,ggofac,
     $   eswork, dmatwk, v )
c
        if( ilocfl.lt.0 ) ilocfl = -ilocfl
c
    5 continue
c
      call MKZERO( 2*nk , v )
c
cpas      call TIMER('After blkvofn    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c   Compute matrix elements off of short-range local potentials
c   indexed by atom "nuc"
c       ves1c() - delta-electrostatic potential
c       pot1c() - fast delta xc-potential
c       xen1c() - fast delta xc energy density
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    Nucleus loop:
c
      do 3000 nuc=1,natm
c
      ktyp = itypa(nuc)
      if( lmx1ctr(ktyp).lt.0 ) goto 3000
c
c    Evaluate the 1-center potentials from nucleus "nuc"
c    It is assumed that vofn is the same for all nuclei
c      equivalent to nuc within a lattice translation
c
      alaminn = almnv1c(ktyp)
c
c    Nuclear lattice loop:
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
      alamini = alamin(ityp)
      rinmax = emax*(alamini+alaminn) / (alamini*alaminn)
c
c    Test overlap of atom i and short-range pot due to "nuc"
c
      if( risquar.gt.rinmax ) goto 2000
c
      call CINIT( ci )
      call COFR( ci,ri )
      if( lmx1ctr(ityp).eq.2 ) call GRCOFR( cigrgr,ri )
c
      iorbat0 = 0
      do 1800 ishl=1,nshli
      nali = nala(ishl,ityp)
      lorbi = lshel(ishl,ityp)
      norbshi = 2*lorbi + 1
      li0 = lorbi*lorbi
      i0 = istart + iorbat0
c find min alpha for this radial function:
      almini = ala(1,ishl,ityp)
      do  ial=1,nali
        ali = ala(ial,ishl,ityp)
        if( ali.lt.almini ) almini = ali
      enddo
      rinmax = emax*(almini+alaminn) / (almini*alaminn)
      if( .not. pjfsrc )then
        alamini = almini
        if( risquar .gt.rinmax ) goto 1800
      endif
c
c    Atom j loop:
c
      jnext = 0
      do 1000 jatm=1,iatm
        jtyp = itypa(jatm)
        jstart = jnext
        jnext = jstart + norba(jtyp)
        nshlj = numshl(jtyp)
c following helps ensure only j.le.i done:
        if( jatm.eq.iatm ) nshlj = ishl
c
c  If any previous "equivalent" triplets, v1cij has already
c  been computed. accordingly, skip this triple of atoms
c
        call EQTRIP( ntest, iateq,jateq,iatm,jatm,nuc,
     $   natm, itypa,ratm )
        if( .not. ntest )then
c         See if jatm,iatm,nuc is equivalent to an earlier triple
c          - again iateq must be greater than or equal to jateq
          call EQTRIP( ntest, iateq,jateq,jatm,iatm,nuc,
     $     natm, itypa,ratm )
        endif
        if( ntest ) goto 1000
c
c    Compile list of equivalent atom triplets
c
        ipoint = 0
        if( nuc.lt.natm )then
          do 50 nuceq=nuc+1,natm
            call GTTRIP( ntest, iateq,jateq,nuceq,iatm,jatm,nuc,
     $       natm, itypa,ratm )
            if( ntest )then
              ipoint = ipoint + 1
              ipermut(ipoint) = 0
            else
            call GTTRIP( ntest, iateq,jateq,nuceq,jatm,iatm,nuc,
     $       natm, itypa,ratm )
              if( ntest )then
                ipoint = ipoint + 1
                ipermut(ipoint) = 1
              endif
            endif
            if( ntest )then
c   (iateq,jateq,nuceq) equivalent to (iatm,jatm,nuc)
              ieqi(ipoint) = iateq
              ieqj(ipoint) = jateq
              ieqn(ipoint) = nuceq
            endif
   50     continue
        endif
c
c  Criterion for including contribution of atom j:
        alaminj = alamin(jtyp)
        ex = risquar*alamini*alaminn / (alamini+alaminn)
        rijnmax = (emax-ex)*(alamini+alaminj+alaminn) /
     $                ((alamini+alaminn)*alaminj)
c
c    Lattice vector loop:
c
        do 900 jlat=1,nlat
c
          do  jr=1,3
            rj(jr) = ratm(jr,jatm) + rlat(jr,jlat) - rnuc(jr)
          enddo
          rjsquar = rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3)
c
c    Test for significant overlap of this atom j with i and nuc
c
          f1 = alamini/(alaminn+alamini)
          rijn = zero
          do  jr=1,3
            rijn = rijn + (f1*ri(jr)-rj(jr))**2
          enddo
          if( rijn.gt.rijnmax ) goto 900
c
          call CINIT( cj )
          call COFR( cj,rj )
c
          jorbat0 = 0
          do 800 jshl=1,nshlj
            nalj = nala(jshl,jtyp)
            lorbj = lshel(jshl,jtyp)
            lj0 = lorbj*lorbj
            norbshj = 2*lorbj+1
            norbshj0 = norbshj
            joff = 0
            j0 = jstart + jorbat0
            if( j0.eq.i0 )then
              norbshj0 = 0
              joff = 1
            endif
c
            lmx = lorbi + lorbj + 1
c
            call MKZERO( neq3d*norbshi*norbshj, v1csum )
            call MKZERO( neq3d*norbshi*norbshj, p1csum )
            call MKZERO( neq3d*norbshi*norbshj, x1csum )
c
c    Alpha loops:
c
            do 320 ial=1,nali
              ali = ala(ial,ishl,ityp)
c      test whether this alpha contributes
              rinmx = emax*(ali+alaminn) / (ali*alaminn)
              if( risquar.gt.rinmx ) goto 320
c
              fi1 = ali/(ali+alaminn)
              rijn = zero
              do  jr=1,3
                rijn = rijn + (fi1*ri(jr)-rj(jr))**2
              enddo
c
              do 310 jal=1,nalj
                alj = ala(jal,jshl,jtyp)
c          test whether this alpha contributes
                ex = risquar*ali*alaminn / (ali+alaminn)
                rijnmx = (emax-ex)*(ali+alj+alaminn) /
     $                ((ali+alaminn)*alj)
                if( rijn.gt.rijnmx ) goto 310
c
                alsum = ali + alj
                exparg = ali*risquar + alj*rjsquar
                do  ir=1,3
                  p(ir) = two*(ali*ri(ir)+alj*rj(ir))
                enddo
                pmag = SQRT( p(1)*p(1)+p(2)*p(2)+p(3)*p(3) )
c
                call NRSCUT( nrmn,nrmx,nrad(ktyp),argmx,
     $                       exparg,pmag,alsum, radmsh(1,ktyp),w1 )
                if( nrmx.eq.0 ) goto 310
c
                nr0 = nrmn - 1
                nrlen = nrmx - nrmn + 1
c
c     Evaluate radial integrals
c
                call MKZERO( 70*neq3d, vradsr )
                call MKZERO( 70*neq3d, pradsr )
                call MKZERO( 70*neq3d, xradsr )
c
                call SRINT( vradsr,ves1c,pradsr,pot1c,xradsr,xen1c,
     $           alsum,pmag,exparg, lmx1ctr(ktyp),lmx,
     $           natm,nuc,ipoint,ieqn,
     $           nrlen,nr0,nrd, radmsh(nrmn,ktyp),radwt(nrmn,ktyp),
     $           besstor,radfac,w1,w2,w3,w4,w5,w6 )
c
c    Evaluate contribution to delta hamiltonian
c
                efac0 = fourpi*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
                lij = 0
                norbshj = norbshj0
                do 420 li=li0+1,li0+norbshi
                  norbshj = norbshj + joff
                  do 410 lj=lj0+1,lj0+norbshj
                    lij = lij + 1
c
                    call ONECTR( velem,vradsr,pelem,pradsr,xelem,xradsr,
     $               ipoint,neq3d,lmx1ctr(ktyp),
     $               alsum,p,li,lj,ci(1,li),cj(1,lj) )
c
                    fac = efac0/(vnorm(li)*vnorm(lj))
                    do  ieq=1,ipoint+1
                      v1csum(ieq,lij) = v1csum(ieq,lij) + fac*velem(ieq)
                      p1csum(ieq,lij) = p1csum(ieq,lij) + fac*pelem(ieq)
                      x1csum(ieq,lij) = x1csum(ieq,lij) + fac*xelem(ieq)
                    enddo
c
  410             continue
  420           continue
c
  310         continue
  320       continue
c
c    Multiply v1csum by bloch phase factor and stuff hamiltonian
c
            norbshj = norbshj0
            i = i0
            lij = 0
            do 520 iorbat=iorbat0+1,iorbat0+norbshi
              i = i + 1
              norbshj = norbshj + joff
              j = j0
              do 510 jorbat=jorbat0+1,jorbat0+norbshj
                lij = lij + 1
                j = j + 1
c
                istuff = i
                jstuff = j
                sign = one
c
                do 500 ieq=1,ipoint+1
                  if( istuff.eq.jstuff )then
                    dmfac = one
                  else
                    dmfac = two
                  endif
                  do  k=1,nk
                    dmatr = dmfac*dmatrd(1,istuff,jstuff,k)
                    exc1c = exc1c + dmatr*x1csum(ieq,lij)*coskr(k,jlat)
                    xp(istuff,jstuff,k) = xp(istuff,jstuff,k) +
     $               p1csum(ieq,lij)*coskr(k,jlat)
                    v(istuff,jstuff,k) = v(istuff,jstuff,k) +
     $               v1csum(ieq,lij)*coskr(k,jlat)
                  enddo
                  if( i.ne.j )then
                    do  k=1,nk
                      dmfac = sign*two*dmatrd(2,jstuff,istuff,k)
                      exc1c = exc1c-dmfac*x1csum(ieq,lij)*sinkr(k,jlat)
                      xp(jstuff,istuff,k) = xp(jstuff,istuff,k) -
     $                 sign*p1csum(ieq,lij)*sinkr(k,jlat)
                      v(jstuff,istuff,k) = v(jstuff,istuff,k) -
     $                 sign*v1csum(ieq,lij)*sinkr(k,jlat)
                    enddo
                  endif
                  if( ieq.le.ipoint )then
                    if( ipermut(ieq).eq.0 )then
                      istuff = ioffa(ieqi(ieq)) + iorbat
                      jstuff = ioffa(ieqj(ieq)) + jorbat
                      sign = one
                    else
                      istuff = ioffa(ieqi(ieq)) + jorbat
                      jstuff = ioffa(ieqj(ieq)) + iorbat
                      sign = -one
                    endif
                  endif
  500           continue
c
  510         continue
  520       continue
c
  800     jorbat0 = jorbat0 + norbshj
c
c    Close center j lattice loop:
c
  900   continue
c
c    Close center j loop:
c
 1000 continue
c
c    Close shell, center i loops:
c
 1800 iorbat0 = iorbat0 + norbshi
 2000 continue
c
c    Close nucleus loops:
c
 2900 continue
 3000 continue
c
c    That's all Folks!
c
      RETURN
      END
