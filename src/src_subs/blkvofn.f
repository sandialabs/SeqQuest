c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BLKVOFN
c
c
      subroutine BLKVOFN( ilocrho, nuc,ktyp, lmx1ctr,almnv1c,
     $ norb,nk, natm,ntyp,nshld,nald, noad, nlat,
     $ itypa, numshl,lshel,nala,ala,cala,  norba,ioffa, occij,
     $ nrs,nrd,radmsh,radwt,
     $ ratm,rlat,vk, nsyms,dsym,naofsym,
     $ dmatrd,
     $ dcountr,rh1coef, cv,ce,  ves1c, pot1c,xen1c,
     $ alxcmin,nalfxc,nalfxcd,alfxc,  alsq,ipvt,
     $ nrxc, vxcrad,excrad,
     $ itypxc, nang, angylm,angpts,angwts, radwf,
     $ nrho, rholoc,vxcloc,excloc, vlm,elm,rfac, vg,eg,
     $ orb,grorb,grgrorb, wf1cwk1,wf1cwk2,wf1cwk3,wf1cwk4,wf1cwk5,
     $ orbj,grorbj,grgrorj,gofac,ggofac,
     $ erfwork, dmatwk, cphi )
c---------------------------------------------------------------
c Purpose: construct local 1-center potentials for nuc given dmat
c
c Revision history:
c  21Jun07-PAS/2.60: removed dsisl calls to clean dependencies
c  11Jun05-PAS/2.59: File units
c   7Mar02-PAS/2.52: ang constants extracted
c   2Jul01-PAS/2.48: replace undirected writes, xcfcn->dftlda
c   1Mar00-PAS/2.43: non-xcfit bugfix
c  21Dec99-PAS/2.42: cosmetic clean-up
c   9Dec98-PAS/2.29: promote source to explicit double precision
c  29Jul98-PAS/2.23: split basis occs
c  28Apr97-PAS/2.17: naming conventions
c  16Oct94-PAS:  reorged for speed,clarity/all-but-real-space-opts
c---------------------------------------------------------------
c
c  There are two contributions to local 1-ctr potentials at nuc from
c  change in density matrix dmat(i,j):
c  First, i and j are both centered on nuc. this term is:
c          sum(   dmat(i,j)* phi(i,r) *phi(j,r)  )
c  and the other term is a two-center contribution (j not on "nuc"):
c          sum( 2.*dmat(i,j) * phi(i,r)* phi(j,r(nuc)) )
c  To avoid divergences in performing the required lattice sums, we
c  make sure that the zeroth, first and second moments of the "1-ctr"
c  charge density vanish exactly.  Thus potential we calculate here
c  falls off rapidly as a function of distance from nucleus "nuc."
c
c  On output, array "ves1c" = short-range es potential for nuc
c  Also computes the quantities "rh1coef" used in "rho1ctr" to
c  subtract the density corresponding to potentials computed here
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c data arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  norba(ntyp),ioffa(natm), occij(nshld,nshld,ntyp)
      DIMENSION  almnv1c(ntyp)
      DIMENSION  radmsh(nrs),radwt(nrs)
      DIMENSION  ratm(3,natm),rlat(3,nlat),vk(3,nk)
      DIMENSION  dsym(5,5,3,*),naofsym(natm,*)
      DIMENSION  dmatrd(2,norb,norb,nk)
complex dmatrd
c xc-fit data
      DIMENSION  alfxc(nalfxc),alsq(nalfxcd,nalfxcd,3),ipvt(nalfxcd,3)
      DIMENSION  vxcrad(nrxc),excrad(nrxc)
      DIMENSION  angylm(nang,9),angpts(3,nang),angwts(nang)
      DIMENSION  radwf(nang,nrd,noad)
c output arrays:
      DIMENSION  dcountr(10),rh1coef(10,noad)
      DIMENSION  cv(10,nalfxc),ce(10,nalfxc)
      DIMENSION  ves1c(nrd,10), pot1c(nrd,10),xen1c(nrd,10)
c scratch arrays:
      DIMENSION  rholoc(nrho,nang),vxcloc(nrxc,nang),excloc(nrxc,nang)
      DIMENSION  vlm(nrxc),elm(nrxc),rfac(nrs,nalfxc)
      DIMENSION  vg(nalfxcd,9),eg(nalfxcd,9)
      DIMENSION  orb(nlat),grorb(3,nlat),grgrorb(3,3,nlat)
      DIMENSION  wf1cwk1(3,nlat),wf1cwk2(3,3,nlat),
     $           wf1cwk3(3,nlat),wf1cwk4(nlat),wf1cwk5(nlat)
      DIMENSION  orbj(noad),grorbj(3,noad),grgrorj(3,3,noad)
      DIMENSION  gofac(nang),ggofac(nang)
      DIMENSION  erfwork(nrs), dmatwk(noad,noad,nlat)
c scratch array used complex only in ROTGRP, declare real here:
      DIMENSION  cphi(*)
c---------------------------------------------------------------
c
c  LOCAL DECLARATIONS:
c
c scratch arrays:
      DIMENSION  eslongr(10),eslfac(10)
      DIMENSION  rij(3)
      DIMENSION  erfiwk(64),eswk(64,2)
c data constant arrays:
      DIMENSION  number(5)
      DATA  number / 1,4,10,20,35 /
c npots has no. of pots in 0<L<3 from i*j orb products
      DIMENSION  npots(9,9)
      DATA  npots
     s/ 1, 1,1,1, 1,1,1, 2,3,
     x  1, 1,1,1, 0,0,0, 0,0,
     y  1, 1,1,1, 0,0,0, 0,0,
     z  1, 1,1,1, 0,0,0, 0,0,
     d  1, 0,0,0, 0,0,0, 0,0,
     d  1, 0,0,0, 0,0,0, 0,0,
     d  1, 0,0,0, 0,0,0, 0,0,
     D  2, 0,0,0, 0,0,0, 0,0,
     D  3, 0,0,0, 0,0,0, 0,0 /
c
c lpots has L<3 pot. no. of off-diag orb prods of ip<8
      DIMENSION  lpots(7,7),lpt(4)
      DATA  lpots
     s/ 1, 2,3,4,  5,6,7,
     x  2, 8,5,7,  0,0,0,
     y  3, 5,9,6,  0,0,0,
     z  4, 7,6,10, 0,0,0,
     d  5, 0,0,0,  0,0,0,
     d  6, 0,0,0,  0,0,0,
     d  7, 0,0,0,  0,0,0/
      DATA  lpt / 8,9,10,8 /
      DIMENSION  lpotgr(3,4),lpgrgr(3,3)
      DATA  lpotgr/2,3,4,
     $             8,5,7,
     $             5,9,6,
     $             7,6,10/
      DATA  lpgrgr/8,5,7,
     $             5,9,6,
     $             7,6,10/
c
      DIMENSION  wght(3,3)
      DATA  wght / 1.,0.,0.,  1.,-1.,0.,  -1.,-1.,2. /
c
      DIMENSION  vnorm(9)
c
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c >>>> EXECUTABLE CODE:
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Get angular normalizations:
      call VNORMS9( vnorm )
      lp1 = lmx1ctr + 1
      lp1sq = lp1*lp1
      lpotmax = number(lp1)
c
c    Zero potentials of nucleus "nuc"
c
      do  lpot=1,lpotmax
        eslongr(lpot) = zero
        call MKZERO( nrs, pot1c(1,lpot) )
        call MKZERO( nrs, xen1c(1,lpot) )
        call MKZERO( nrs, ves1c(1,lpot) )
      enddo
c
c    Initialize accumulators for Taylor series coefficients.
c
      do  iorbat=1,norba(ktyp)
        orbj(iorbat) = zero
        do  ir=1,3
          grorbj(ir,iorbat) = zero
          do  jr=1,3
            grgrorj(ir,jr,iorbat) = zero
          enddo
        enddo
      enddo
c
c    Initialize density of atom "nuc" on spherical mesh:
c
      if( ilocrho.gt.0 )then
c       Initialize density with tails from other atoms
        read(ilocrho)  rholoc
      else
c       Start with a blank slate
        do  nr=1,nrxc
          do  iang=1,nang
            rholoc(nr,iang) = zero
          enddo
        enddo
      endif
c
c    Set error function for alpha minimum
c
      alvmin = almnv1c(ktyp)
      alxmin =  MAX( alxcmin, almnv1c(ktyp) )
      rtalvmn = SQRT( alvmin )
      do  nr=1,nrs
        erfwork(nr) = DERF( rtalvmn*radmsh(nr) ) / radmsh(nr)
      enddo
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Begin Atom j loop
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do 1000 jatm=1,natm
        jtyp = itypa(jatm)
c
c    Sum density matrix over surface brillouin zone
c
        call MKZERO( nlat*noad**2, dmatwk )
        do 20 klat=1,nlat
c
          do  jr=1,3
            rij(jr) = ratm(jr,nuc) - ratm(jr,jatm) - rlat(jr,klat)
          enddo
c
          call ROTGRP( dmatwk(1,1,klat), nuc,jatm,rij,
     $     norb,nk, natm,ntyp,nshld, noad,
     $     itypa,numshl,lshel, ioffa, vk,rlat(1,klat),
     $     nsyms,dsym,naofsym,
     $     dmatrd,cphi )
c
   20   continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              begin rh1coef calculation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        jorbat = 0
        do 75 jshl=1,numshl(jtyp)
          nalj = nala(jshl,jtyp)
          lorbj = lshel(jshl,jtyp)
          lj = lorbj**2
          norbshj = 2*lorbj + 1
          do 65 jorbsh=1,norbshj
            jorbat = jorbat + 1
            lj = lj + 1
c
c    Obtain orb, grorb, and grgrorb for all nlat:
c
            call WF1COFR( orb,grorb,grgrorb,
     $       nuc,jatm, lj,jtyp,jshl,nalj,
     $       natm,ntyp,nshld,nald, nlat, ala,cala, ratm,rlat,
     $       wf1cwk1,wf1cwk2,wf1cwk3,wf1cwk4,wf1cwk5 )
c
            iorbat = 0
            do 70 ishl=1,numshl(ktyp)
              nali = nala(ishl,ktyp)
              lorbi = lshel(ishl,ktyp)
              norbshi = 2*lorbi + 1
              do 60 iorbsh=1,norbshi
                iorbat = iorbat + 1
c
                do 50 klat=1,nlat
c
                  if( jatm.eq.nuc .and. klat.eq.1 )then
c
c    Accumulate 1-atom part of density to "rholoc"
c
                    if( itypxc.ne.0 )then
                      dm = dmatwk(iorbat,jorbat,1)
                      do  iang=1,nang
                        do  nr=1,nrxc
                          rholoc(nr,iang) = rholoc(nr,iang) + dm*
     $                     radwf(iang,nr,iorbat)*radwf(iang,nr,jorbat)
                        enddo
                      enddo
                    endif
c                   Skip rest as jatm and nuc are same atom
                    goto 50
c
                  endif
c
c    Add approximate contributions to "rholoc" for 2-atom case
c
c  Accumulate coefficients for Taylor series expansions about site iatm
c  of orbitals j, multplied by 2*dmat(i,j).  Note 2nd derivative in
c  in Taylor series includes factor of 1/2, cancelling this factor of 2
c  This information is used to complete "rholoc" calculation after
c  completion of the "do 1000" loop, below
c
                  dm = dmatwk(iorbat,jorbat,klat)
                  dmm = dm + dm
                  orbj(iorbat) = orbj(iorbat) + dmm*orb(klat)
                  do  ir=1,3
                    grorbj(ir,iorbat) = grorbj(ir,iorbat)
     $               + dmm*grorb(ir,klat)
                    do  jr=1,3
                      grgrorj(ir,jr,iorbat) = grgrorj(ir,jr,iorbat)
     $                  + dm*grgrorb(ir,jr,klat)
                    enddo
                  enddo
c
c Accumulate 2-atom contributions to "rh1coef" ...
c
                  if( lorbi .gt. lmx1ctr ) goto 60
c
c  ... from orbital j at site
c
                  rh1coef(1,iorbat) = rh1coef(1,iorbat) +
     $             dmatwk(iorbat,jorbat,klat)*orb(klat)
c
                  if( lorbi+1 .gt. lmx1ctr ) goto 50
c
c  ... from gradient of orbital j
c
                  do  lgrad=1,3
                    rh1coef(1+lgrad,iorbat) = rh1coef(1+lgrad,iorbat)
     $               + dmatwk(iorbat,jorbat,klat)*grorb(lgrad,klat)
                  enddo
c
                  if( lorbi+2 .gt. lmx1ctr ) goto 50
c
c  ... from grad(grad(orbital j).
c      (since lmx1ctr restricted to .le.2, lorbi must=0 to use)
c
                  lcoef = 4
                  do  lgr=1,3
                    do  llgr=1,lgr
                      if( llgr.eq.lgr )then
                        offd = half
                      else
                        offd = one
                      endif
                      lcoef = lcoef + 1
                      rh1coef(lcoef,iorbat) =rh1coef(lcoef,iorbat)+offd*
     $                 dmatwk(iorbat,jorbat,klat)*grgrorb(llgr,lgr,klat)
                    enddo
                  enddo
c
c    Close lattice vector loop:
c
   50           continue
c
c    Close i-orbital loops:
c
   60         continue
   70       continue
c
c    Close j-orbital loops:
c
   65     continue
   75   continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               end rh1coef calculation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                begin vofn calculation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( jatm.eq.nuc )then
c
c **********************  1 atom case  **********************
c
      iorbat0 = 0
      do 182 ishl=1,numshl(ktyp)
        nali = nala(ishl,ktyp)
        lorbi = lshel(ishl,ktyp)
        norbshi = 2*lorbi + 1
        li0 = lorbi*lorbi
c
        jorbat0 = 0
        do 181 jshl=1,numshl(ktyp)
          nalj = nala(jshl,ktyp)
          lorbj = lshel(jshl,ktyp)
          norbshj = 2*lorbj + 1
          lj0 = lorbj*lorbj
c
          lorbsum = lorbi + lorbj
c
          if( lorbsum .gt. lmx1ctr ) goto 181
c
          esll = zero
          esl0 = zero
          do  nr=1,nrs
            eswk(nr,1) = zero
            eswk(nr,2) = zero
          enddo
c
          if( ishl.gt.jshl )then
            occ = occij(ishl,jshl,ktyp)
          else
            occ = occij(jshl,ishl,ktyp)
          endif
c
c    Alpha loops:
c
          do 162 ial=1,nali
            ali = ala(ial,ishl,ktyp)
            cali = cala(ial,ishl,ktyp)
c
            do 161 jal=1,nalj
              alj = ala(jal,jshl,ktyp)
              calj = cala(jal,jshl,ktyp)
c
              als = ali + alj
              if( als.lt.alvmin ) goto 161
c
              twoals = als + als
              rtals = SQRT( als )
              rtpial = rtpi/rtals
              expfac = two/rtpial
              alfac = cali*calj*(rtpial**3)/(twoals**lorbsum)
c
              esll = esll + alfac
              esl0 = esl0 + twoals*alfac
c
c  Last statements account for fact that x**2 potential (y**2 and
c  z**2 too) equivalent to a (d/dx)**2 plus an s-like potential
c
              alfac = two*alfac
              alfm2 = twoals*alfac
c generate and save unvectorizeable functions of r:
              do  nr=1,nrs
                r = radmsh(nr)
                erfiwk(nr) = DERF(rtals*r)/r
              enddo
c
c  Calculate the radial parts of the potentials:
c
              if( lorbsum.eq.0 )then
                do  nr=1,nrs
                  esfac = erfiwk(nr) - erfwork(nr)
                  eswk(nr,1) = eswk(nr,1) + alfac*esfac
                enddo
              elseif( lorbsum.eq.1 )then
                do  nr=1,nrs
                  esfac = erfiwk(nr) - erfwork(nr)
                  eswk(nr,1) = eswk(nr,1) + alfac*esfac
                enddo
              elseif( lorbsum.eq.2 )then
                do  nr=1,nrs
                  r = radmsh(nr)
                  esfac = erfiwk(nr) - erfwork(nr)
                  eswk(nr,1) = eswk(nr,1) + alfac*esfac
                  eswk(nr,2) = eswk(nr,2) + alfm2*esfac
                enddo
              endif
c
c    Close alpha loops:
c
  161       continue
  162     continue
c
c    Orbital loops:
c
          iorbat = iorbat0
          li = li0
          do 172 iorbsh=1,norbshi
            iorbat = iorbat + 1
            li = li + 1
c
            jorbat = jorbat0
            lj = lj0
            do 171 jorbsh=1,norbshj
              jorbat = jorbat + 1
              lj = lj + 1
c
c  Determine how many potentials this pair of orbitals contributes
c  to, what l-values, and with what weight
c
              numpot = npots(li,lj)
              if( numpot.eq.0 ) goto 171
c
              dmatr = dmatwk(iorbat,jorbat,1)
              if( li.eq.lj ) dmatr = dmatr - occ
              dmfac = dmatr/(vnorm(li)*vnorm(lj))
c
              efacl = dmfac*esll
              efac0 = dmfac*esl0
c
              if( numpot.eq.1 )then
                lpot = lpots(li,lj)
                if( lpot.gt.lpotmax ) goto 171
                vfac = dmfac
                eslongr(lpot) = eslongr(lpot) + efacl
                if( lpot.lt.8 )then
                  do  nr=1,nrs
                    ves1c(nr,lpot) = ves1c(nr,lpot) + vfac*eswk(nr,1)
                  enddo
                else
                  eslongr(1) = eslongr(1) + efac0
                  do  nr=1,nrs
                    ves1c(nr,lpot) = ves1c(nr,lpot) + vfac*eswk(nr,1)
                    ves1c(nr,1) = ves1c(nr,1) + vfac*eswk(nr,2)
                  enddo
                endif
              else
                do 155 npt=1,numpot
                  wt = wght(npt,numpot)
                  lpot = lpt(npt)
                  if( lpot.gt.lpotmax ) goto 171
                  if( wt.eq.zero ) goto 155
                  vfac = wt*dmfac
                  eslongr(lpot) = eslongr(lpot) + wt*efacl
                  if( lpot.lt.8 )then
                    do  nr=1,nrs
                      ves1c(nr,lpot) = ves1c(nr,lpot) + vfac*eswk(nr,1)
                    enddo
                  else
                    eslongr(1) = eslongr(1) + wt*efac0
                    do  nr=1,nrs
                      ves1c(nr,lpot) = ves1c(nr,lpot) + vfac*eswk(nr,1)
                      ves1c(nr,1) = ves1c(nr,1) + vfac*eswk(nr,2)
                    enddo
                  endif
  155           continue
              endif
c
c    Close orbital loops:
c
  171       continue
  172     continue
c
c    Close radial function loops:
c
  181   jorbat0 = jorbat0 + norbshj
  182 iorbat0 = iorbat0 + norbshi
c
      endif
c
c **********************  2 atom case  **********************
c
      iorbat0 = 0
      do 282 ishl=1,numshl(ktyp)
        nali = nala(ishl,ktyp)
        lorbi = lshel(ishl,ktyp)
        norbshi = 2*lorbi + 1
        li0 = lorbi*lorbi
c
c    Alpha i loop:
c
      do 262 ial=1,nali
        ali = ala(ial,ishl,ktyp)
        cali = cala(ial,ishl,ktyp)
        if( ali.lt.alvmin ) goto 262
c
        do  lpot=1,10
          eslfac(lpot) = zero
        enddo
c
c    Get value of jorbat-th orbital at nucleus "nuc"
c
        twoali = ali + ali
        jorbat = 0
        do 281 jshl=1,numshl(jtyp)
          nalj = nala(jshl,jtyp)
          lorbj = lshel(jshl,jtyp)
          lj = lorbj*lorbj
          norbshj = 2*lorbj + 1
          do 272 jorbsh=1,norbshj
            jorbat = jorbat + 1
            lj = lj + 1
c
c    Obtain orb, grorb, and grgrorb for all nlat
c
            call WF1COFR( orb,grorb,grgrorb,
     $       nuc,jatm, lj,jtyp,jshl,nalj,
     $       natm,ntyp,nshld,nald, nlat, ala,cala, ratm,rlat,
     $       wf1cwk1,wf1cwk2,wf1cwk3,wf1cwk4,wf1cwk5 )
c
c    Orbital i loop:
c
            iorbat = iorbat0
            li = li0
            do 271 iorbsh=1,norbshi
              iorbat = iorbat + 1
              li = li + 1
c
c    Loop for lattice vector j
c
              do 290 klat=1,nlat
c
c    Skip if jatm and nuc are same atom
c
                if( klat.eq.1 .and. jatm.eq.nuc ) goto 290
c
                dmfac = two*dmatwk(iorbat,jorbat,klat)/vnorm(li)
                rtali = SQRT( ali )
                rtpial = rtpi/rtali
                cofac = dmfac*two*cali*(rtpial**3)/(twoali**lorbi)
                coeff = cofac*orb(klat)
c
c    Accumulate counter-charge:
c
                numpot = npots(li,1)
                if( numpot.eq.1 )then
                  lpot = lpots(li,1)
                  if( lpot.gt.lpotmax ) goto 272
                  eslfac(lpot) = eslfac(lpot) + coeff
                else
                  do  npt=1,numpot
                    wt = wght(npt,numpot)
                    lpot = lpt(npt)
                    if( lpot.gt.lpotmax ) goto 272
                    eslfac(lpot) = eslfac(lpot) + coeff*wt
                  enddo
                endif
c
                if( lorbi+1 .le. lmx1ctr )then
c
c      ... from gradient of orbital j at nucleus i
c
                  do  lgrad=1,3
                    coeff = cofac*grorb(lgrad,klat)/twoali
                    lpot = lpotgr(lgrad,li)
                    eslfac(lpot) = eslfac(lpot) + coeff
                  enddo
                endif
c
                if( lorbi+2 .le. lmx1ctr )then
c
c      ... from grad(grad(orbital j).
c
                  do  lgr=1,3
                    do  llgr=1,lgr
                      if( llgr.eq.lgr )then
                        offd = half
                      else
                        offd = one
                      endif
                      coeff =cofac*offd*grgrorb(llgr,lgr,klat)/twoali**2
                      lpot = lpgrgr(llgr,lgr)
                      eslfac(lpot) = eslfac(lpot) + coeff
                    enddo
                  enddo
                endif
c
c    Close lattice vector loop:
c
  290         continue
c
c    Close orbital i loop:
c
  271       continue
c
c    Close orbital j loop:
c
  272     continue
c
c    Close jshl loop:
c
  281   continue
c
        do  lpot=1,10
          if( eslfac(lpot).ne.zero ) goto 420
        enddo
        goto 262
c
c    Add potential to ves1c:
c
  420   sesl = zero
        do  lpot=8,10
          sesl = sesl + twoali*eslfac(lpot)
        enddo
c
        do  nr=1,nrs
          erfiwk(nr) = DERF(rtali*radmsh(nr))/radmsh(nr) - erfwork(nr)
        enddo
c
        eslfac(1) = eslfac(1) + sesl
        do 440 lpot=1,10
          eslp = eslfac(lpot)
          if( eslp.ne.zero )then
            eslongr(lpot) = eslongr(lpot) + half*eslp
            do  nr=1,nrs
              ves1c(nr,lpot) = ves1c(nr,lpot) + eslp*erfiwk(nr)
            enddo
          endif
  440   continue
c
c    Close alpha i loop:
c
  262 continue
c
c    Close ishl loop:
c
  282 iorbat0 = iorbat0 + norbshi
c
c  Calculate coefficients of counter charges to cancel
c  long-ranged portions of the 1-center potential.
c
        do 950 lpot=1,lpotmax
          dcountr(lpot) = eslongr(lpot)
          if( lpot.lt.8 ) goto 950
          dcountr(1) = dcountr(1) - two*alvmin*eslongr(lpot)
  950   continue
c
c    Close atom j loop:
c
 1000 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Finished with electrostatic potentials
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Check if we have an xc-fit potential to build:
      if( itypxc.eq.0 ) RETURN
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Complete accumulation of rholoc, the density
c        associated with atom "nuc" on a spherical mesh.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize xc-fit Gaussian coefficients:
      do  ialfxc=1,nalfxc
        do  ilxc=1,10
          cv(ilxc,ialfxc) = zero
          ce(ilxc,ialfxc) = zero
        enddo
      enddo
c
c  To compute contribution of terms where one atom is not "nuc,"
c  multiply Taylor series coefficients by angular factors.
c
      do 1070 iorbat=1,norba(ktyp)
        do  iang=1,nang
          x = angpts(1,iang)
          y = angpts(2,iang)
          z = angpts(3,iang)
          twoy = y + y
          twoz = z + z
          gofac(iang) = grorbj(1,iorbat)*x + grorbj(2,iorbat)*y +
     $     grorbj(3,iorbat)*z
          ggofac(iang) = x*( x*grgrorj(1,1,iorbat) +
     $     twoy*grgrorj(2,1,iorbat) + twoz*grgrorj(3,1,iorbat) ) +
     $     y*( y*grgrorj(2,2,iorbat) + twoz*grgrorj(3,2,iorbat) ) +
     $     z*z*grgrorj(3,3,iorbat)
        enddo
c
c  Multiply by orbital i and radial factors then accumulate density
c
        do  nr=1,nrxc
          r = radmsh(nr)
          rsq = r*r
          do  iang=1,nang
            orbi = radwf(iang,nr,iorbat)
            orbfac = orbj(iorbat) + r*gofac(iang) + rsq*ggofac(iang)
            rholoc(nr,iang) = rholoc(nr,iang) + orbi*orbfac
          enddo
        enddo
c
 1070 continue
c
c  (1/4/90)-PJF: In this attempt, I deliberately do not add to rholoc
c  the contributions from the tails of the neighboring atoms.
c  This contribution would require computing additional dmatwk matrix
c  elements, and the resulting contribution to the density is small.
c  20Feb93-PAS: tails of nearby atoms approximated by adding
c  density from initial occupations (computed by ATMRHO in setup)
c
c  At this point, the density on the spherical mesh around atom "iatm"
c  has been accumulated in "rholoc".  Now evaluate the exchange
c  correlation potential and energy density.
c
      rhat0 = rholoc(1,1)
      rhat = zero
      rnorm = zero
      do  nr=1,nrxc
        do  iang=1,nang
          rhat = rhat + radwt(nr)*rholoc(nr,iang)*angwts(iang)
          rnorm = rnorm + radwt(nr)*angwts(iang)*angylm(iang,1)
          call DFTLDA( rholoc(nr,iang), vxcloc(nr,iang),excloc(nr,iang))
        enddo
      enddo
c
c  Now project out l,m components of vxcloc and excloc using values
c  of real ylm's on unit sphere in array angylm.
c
      lval = 1
      do 1200 lm=1,lp1sq
c
        if( lm.eq.1 )then
          do  ialfxc=1,nalfxc
            do  nr=1,nrxc
              r = radmsh(nr)
              rfac(nr,ialfxc) = radwt(nr)*EXP( -alfxc(ialfxc)*r*r )
            enddo
          enddo
          vn1 = vnorm(1)
c  Subtract isolated atom spherical parts
          do  nr=1,nrxc
            vlm(nr) = -vn1*vxcrad(nr)
            elm(nr) = -vn1*excrad(nr)
          enddo
        else
          do  nr=1,nrxc
            vlm(nr) = zero
            elm(nr) = zero
          enddo
        endif
c
        do  iang=1,nang
          angfac = angwts(iang)*angylm(iang,lm)
          do  nr=1,nrxc
            vlm(nr) = vlm(nr) + angfac*vxcloc(nr,iang)
            elm(nr) = elm(nr) + angfac*excloc(nr,iang)
          enddo
        enddo
cpas:start
        if( lm.eq.1 )then
          vreq0 = vxcloc(1,1)
          ereq0 = excloc(1,1)
          vleq0 = vlm(1) / vnorm(1)
          eleq0 = elm(1) / vnorm(1)
C         write(IWR,1313)  nuc,rhat,rhat0,vreq0,ereq0,vleq0,eleq0
 1313 format(1x,'NUC',i2,' elns=',2f12.8/1x,'vxc(0),exc(0)=',4f12.8)
        endif
cpas:end
c
c  Now least squares fit vlm and elm to linear combinations of
c  gaussians times r**l, by solving AX=B.  Compute B vectors:
c
        do  ialfxc=1,nalfxc
          vg(ialfxc,lm) = zero
          eg(ialfxc,lm) = zero
          do  nr=1,nrxc
            vg(ialfxc,lm) = vg(ialfxc,lm) + vlm(nr)*rfac(nr,ialfxc)
            eg(ialfxc,lm) = eg(ialfxc,lm) + elm(nr)*rfac(nr,ialfxc)
          enddo
        enddo
c
c   Using saved matrices for this atom type, solve AX=B problem to
c   determine best fit coefficients for vxc and exc.
c
clib: disable calls to linpack/dsisl
        call STOPXERR( 'DSIABLED - removed linpack/dsisl calls' )
clib        call DSISL(alsq(1,1,lval),nalfxcd,nalfxc,ipvt(1,lval),vg(1,lm))
clib        call DSISL(alsq(1,1,lval),nalfxcd,nalfxc,ipvt(1,lval),eg(1,lm))
c
        if( lm.eq.1 .or. lm.eq.4 )then
          lval = lval + 1
          do  ialfxc=1,nalfxc
            do  nr=1,nrxc
              rfac(nr,ialfxc) = rfac(nr,ialfxc)*radmsh(nr)
            enddo
          enddo
        endif
c
 1200 continue
c
c  This ends lm loop.  Output is set of coefficients "vg" and "eg" for
c  the xc potential and energy density about atom "nuc."
c
c  Assemble rapidly varying components of Vxc and Exc on radial mesh,
c  using fit coefficients corresponding to the Gaussians that vary
c  rapidly.  That is, only use the coefficients in the just completed
c  fit corresponding to alphas greater than alv1cmn.
c
      do 1500 ialfxc=1,nalfxc
        alfi = alfxc(ialfxc)
        twoalfi = alfi + alfi
        if( alfi.lt.alxmin ) goto 1500
c
c      First reorganize coeffs so that potentials are proportional to:
c      1,x,y,z,xy,yz,zx,x**2,y**2,z**2
c
        cv( 1,ialfxc) = vg(ialfxc,1)/vnorm(1)
        ce( 1,ialfxc) = eg(ialfxc,1)/vnorm(1)
c
        if( lpotmax.gt.1 )then
          cv( 2,ialfxc) = -vg(ialfxc,2)/vnorm(2)
          cv( 3,ialfxc) = -vg(ialfxc,3)/vnorm(2)
          cv( 4,ialfxc) =  vg(ialfxc,4)/vnorm(2)
          ce( 2,ialfxc) = -eg(ialfxc,2)/vnorm(2)
          ce( 3,ialfxc) = -eg(ialfxc,3)/vnorm(2)
          ce( 4,ialfxc) =  eg(ialfxc,4)/vnorm(2)
        endif
c
        if( lpotmax.gt.4 )then
          cv( 5,ialfxc) =  vg(ialfxc,5)/vnorm(5)
          cv( 6,ialfxc) = -vg(ialfxc,6)/vnorm(5)
          cv( 7,ialfxc) = -vg(ialfxc,7)/vnorm(5)
          vg1334 = vg(ialfxc,8)/vnorm(8)
          vg1135 = vg(ialfxc,9)/vnorm(9)
          cv( 8,ialfxc) =   vg1334 - vg1135
          cv( 9,ialfxc) = - vg1334 - vg1135
          cv(10,ialfxc) =   vg1135 + vg1135
          ce( 5,ialfxc) =  eg(ialfxc,5)/vnorm(5)
          ce( 6,ialfxc) = -eg(ialfxc,6)/vnorm(5)
          ce( 7,ialfxc) = -eg(ialfxc,7)/vnorm(5)
          eg1334 = eg(ialfxc,8)/vnorm(8)
          eg1135 = eg(ialfxc,9)/vnorm(9)
          ce( 8,ialfxc) =   eg1334 - eg1135
          ce( 9,ialfxc) = - eg1334 - eg1135
          ce(10,ialfxc) =   eg1135 + eg1135
        endif
c
        do  nr=1,nrs
          r = radmsh(nr)
          rfac(nr,ialfxc) = EXP( -alfxc(ialfxc)*r*r )
        enddo
c
        do 1450 lpot=1,lpotmax
          do  nr=1,nrs
            ves1c(nr,lpot) = ves1c(nr,lpot)
     $       + cv(lpot,ialfxc)*rfac(nr,ialfxc)
            pot1c(nr,lpot) = pot1c(nr,lpot)
     $       + cv(lpot,ialfxc)*rfac(nr,ialfxc)
            xen1c(nr,lpot) = xen1c(nr,lpot)
     $       + ce(lpot,ialfxc)*rfac(nr,ialfxc)
          enddo
c
          if( lpot.eq.1 .or. lpot.eq.4 )then
            do  nr=1,nrs
c might need a minus here:
              rfac(nr,ialfxc) = rfac(nr,ialfxc)/twoalfi
            enddo
          endif
c
 1450   continue
c
c    Close loop over fit alphas:
c
 1500 continue
c
c    That's all Folks!
c
      RETURN
      END
