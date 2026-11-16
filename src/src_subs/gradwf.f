c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRADWF
c
c
      subroutine GRADWF( ibox, ncplx, ndim, iatm,ishl,iorbsh, almin,
     $ nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $ itypa, numshl,lshel,nala,ala,cala, ratm, rlat, coskr,sinkr,
     $ nrbox, boxrad,boxctr,
     $ rv, rsq, sumg,sum2ag, dpdx,dpdy,dpdz,
     $ dphix,dphiy,dphiz,  dphiuv )
c---------------------------------------------------------------
c Purpose: grid wave function gradients for force calculation
c
c Written: Peter A. Schultz, 28-May-1997, for v2.20 (boxes)
c
c Revision history:
c  23Jan09-PAS/2.63: k-parallel dmat matrices
c   7Mar02-PAS/2.52: ang constants extracted
c  21Dec99-PAS/2.42: cosmetic clean-up
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  14Jun99-PAS/2.36: rework grid force memory
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  dphix,dphiy,dphiz,  dphiuv
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  ratm(3,iatm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  boxctr(3), rv(3,nrbox)
c output arrays:
      DIMENSION  dphix(*),dphiy(*),dphiz(*)
      DIMENSION  dphiuv(ndim,ndim,*)
c scratch arrays:
      DIMENSION  rsq(nrbox), sumg(nrbox),sum2ag(nrbox)
      DIMENSION  dpdx(nrbox),dpdy(nrbox),dpdz(nrbox)
c local declarations:
      DIMENSION  ratom(3), rlast(3),rshift(3)
c
      DIMENSION  vnorm(5)
      DATA  zero,one,two, four / 0.d0,1.d0,2.d0, 4.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm)
c
      ityp = itypa(iatm)
c
      lshl = lshel(ishl,ityp)
      nali = nala(ishl,ityp)
c
      vn = vnorm(lshl+1)
      angfac = one
c
      do  jr=1,3
        rlast(jr) = zero
      enddo
c
c    Clear out wave function array:
c
          lenv = ncplx*nkloc*nrbox
          call MKSZERO( lenv, dphix )
          call MKSZERO( lenv, dphiy )
          call MKSZERO( lenv, dphiz )
          if( ndim.ne.0 )then
            call MKSZERO( ndim*ndim*lenv, dphiuv )
          endif
c
          irorb = 0
c
          almins = ala(1,ishl,ityp)
          do  ial=2,nali
            ali = ala(ial,ishl,ityp)
            if( ali.lt.almin ) almins = ali
          enddo
          rngmax = SQRT( cutgrid / (almin +almins) )
c
c  Cycle through lattice vectors:
c
          do 900 ilat=1,nlat
c
            rngbox = zero
            do  jr=1,3
              ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
              rngbox = rngbox + (ratom(jr)-boxctr(jr))**2
            enddo
            rngbox = SQRT( rngbox ) - boxrad
            if( rngbox.gt.rngmax ) goto 900
c
            rsqbox = rngbox**2
            if( rngbox.lt.zero ) rsqbox = zero
c
            do  jr=1,3
              rshift(jr) = ratom(jr) - rlast(jr)
            enddo
            do  jr=1,3
              rlast(jr) = ratom(jr)
            enddo
c
c  Shift coordinates to this lattice image:
c
            do  nr=1,nrbox
              rv(1,nr) = rv(1,nr) - rshift(1)
              rv(2,nr) = rv(2,nr) - rshift(2)
              rv(3,nr) = rv(3,nr) - rshift(3)
              rsq(nr) = rv(1,nr)**2 + rv(2,nr)**2 + rv(3,nr)**2
              sumg(nr) = zero
              sum2ag(nr) = zero
            enddo
c
            ndo = 0
            do 70 ial=1,nali
              alf = ala(ial,ishl,ityp)
              twoalc = alf + almin
              rsqcut = cutgrid/twoalc
              if( rsqbox.gt.rsqcut ) goto 70
c
              twoalf = alf + alf
              coef = cala(ial,ishl,ityp)/vn
              ndo = ndo + 1
              do  nr=1,nrbox
                expar2 = coef*EXP( -alf*rsq(nr) )
                sumg(nr) = sumg(nr) + expar2
                sum2ag(nr) = sum2ag(nr) + twoalf*expar2
              enddo
   70       continue
            if( ndo.eq.0 ) goto 900
c
            goto (100,110,120), lshl+1
c
c        -----  s-state  -----
c
  100       continue
c       S-function / iorbsh = 1
            do  nr=1,nrbox
              dpdx(nr) = - rv(1,nr)*sum2ag(nr)
              dpdy(nr) = - rv(2,nr)*sum2ag(nr)
              dpdz(nr) = - rv(3,nr)*sum2ag(nr)
            enddo
            goto 400
c
c        -----  p-state  -----
c
  110       continue
            goto (111,112,113) iorbsh
c
  111       continue
c       X-function / iorbsh = 1
            do  nr=1,nrbox
              p = rv(1,nr)
              dpdx(nr) = sumg(nr) - rv(1,nr)*p*sum2ag(nr)
              dpdy(nr) =          - rv(2,nr)*p*sum2ag(nr)
              dpdz(nr) =          - rv(3,nr)*p*sum2ag(nr)
            enddo
            goto 400
c
  112       continue
c       Y-function / iorbsh = 2
            do  nr=1,nrbox
              p = rv(2,nr)
              dpdx(nr) =          - rv(1,nr)*p*sum2ag(nr)
              dpdy(nr) = sumg(nr) - rv(2,nr)*p*sum2ag(nr)
              dpdz(nr) =          - rv(3,nr)*p*sum2ag(nr)
            enddo
            goto 400
c
  113       continue
c       Z-function / iorbsh = 3
            do  nr=1,nrbox
              p = rv(3,nr)
              dpdx(nr) =          - rv(1,nr)*p*sum2ag(nr)
              dpdy(nr) =          - rv(2,nr)*p*sum2ag(nr)
              dpdz(nr) = sumg(nr) - rv(3,nr)*p*sum2ag(nr)
            enddo
            goto 400
c
c        -----  d-state  -----
c
  120       continue
            goto (121,122,123,124,125) iorbsh
c
  121       continue
c       XY-function / iorbsh = 1
            do  nr=1,nrbox
              x = rv(1,nr)
              y = rv(2,nr)
              z = rv(3,nr)
              dpdx(nr) = y*(sumg(nr) - x*x*sum2ag(nr))
              dpdy(nr) = x*(sumg(nr) - y*y*sum2ag(nr))
              dpdz(nr) =           - x*y*z*sum2ag(nr)
            enddo
            goto 400
c
  122       continue
c       YZ-function / iorbsh = 2
            do  nr=1,nrbox
              x = rv(1,nr)
              y = rv(2,nr)
              z = rv(3,nr)
              dpdx(nr) =           - y*z*x*sum2ag(nr)
              dpdy(nr) = z*(sumg(nr) - y*y*sum2ag(nr))
              dpdz(nr) = y*(sumg(nr) - z*z*sum2ag(nr))
            enddo
            goto 400
c
  123       continue
c       ZX-function / iorbsh = 3
            do  nr=1,nrbox
              x = rv(1,nr)
              y = rv(2,nr)
              z = rv(3,nr)
              dpdx(nr) = z*(sumg(nr) - x*x*sum2ag(nr))
              dpdy(nr) =           - z*x*y*sum2ag(nr)
              dpdz(nr) = x*(sumg(nr) - z*z*sum2ag(nr))
            enddo
            goto 400
c
  124       continue
c       X2-Y2 function / iorbsh = 4
            angfac = vnorm(3)/vnorm(4)
            do  nr=1,nrbox
              x = rv(1,nr)
              y = rv(2,nr)
              z = rv(3,nr)
              d = x*x - y*y
              dpdx(nr) = x*( two*sumg(nr) - d*sum2ag(nr))
              dpdy(nr) = y*(-two*sumg(nr) - d*sum2ag(nr))
              dpdz(nr) =                - z*d*sum2ag(nr)
            enddo
            goto 400
c
  125       continue
c       Z2-(X2+Y2)/2 = 3R2-Z2 function / iorbsh = 5
            angfac = vnorm(3)/vnorm(5)
            do  nr=1,nrbox
              x = rv(1,nr)
              y = rv(2,nr)
              z = rv(3,nr)
              d = two*z*z - x*x - y*y
              dpdx(nr) = x*(-two*sumg(nr) - d*sum2ag(nr))
              dpdy(nr) = y*(-two*sumg(nr) - d*sum2ag(nr))
              dpdz(nr) = z*(four*sumg(nr) - d*sum2ag(nr))
            enddo
            goto 400
c
c  Compute wave-function derivatives
c
  400       continue
c
            if( ncplx.eq.1 )then
c
              ir0 = irorb
              do  nr=1,nrbox
                dphix(ir0+nr) = dphix(ir0+nr) + angfac*dpdx(nr)
                dphiy(ir0+nr) = dphiy(ir0+nr) + angfac*dpdy(nr)
                dphiz(ir0+nr) = dphiz(ir0+nr) + angfac*dpdz(nr)
              enddo
c
              ir0 = irorb
              if( ndim.eq.1 )then
                do  nr=1,nrbox
                  dphiuv(1,1,ir0+nr) = dphiuv(1,1,ir0+nr)
     $             + angfac*dpdx(nr)*rv(1,nr)
                enddo
              elseif( ndim.eq.2 )then
                do  nr=1,nrbox
                  dphiuv(1,1,ir0+nr) = dphiuv(1,1,ir0+nr)
     $             + angfac*dpdx(nr)*rv(1,nr)
                  dphiuv(2,1,ir0+nr) = dphiuv(2,1,ir0+nr)
     $             + angfac*dpdy(nr)*rv(1,nr)
                  dphiuv(1,2,ir0+nr) = dphiuv(1,2,ir0+nr)
     $             + angfac*dpdx(nr)*rv(2,nr)
                  dphiuv(2,2,ir0+nr) = dphiuv(2,2,ir0+nr)
     $             + angfac*dpdy(nr)*rv(2,nr)
                enddo
              elseif( ndim.eq.3 )then
                do  nr=1,nrbox
                  dphiuv(1,1,ir0+nr) = dphiuv(1,1,ir0+nr)
     $             + angfac*dpdx(nr)*rv(1,nr)
                  dphiuv(2,1,ir0+nr) = dphiuv(2,1,ir0+nr)
     $             + angfac*dpdy(nr)*rv(1,nr)
                  dphiuv(3,1,ir0+nr) = dphiuv(3,1,ir0+nr)
     $             + angfac*dpdz(nr)*rv(1,nr)
                  dphiuv(1,2,ir0+nr) = dphiuv(1,2,ir0+nr)
     $             + angfac*dpdx(nr)*rv(2,nr)
                  dphiuv(2,2,ir0+nr) = dphiuv(2,2,ir0+nr)
     $             + angfac*dpdy(nr)*rv(2,nr)
                  dphiuv(3,2,ir0+nr) = dphiuv(3,2,ir0+nr)
     $             + angfac*dpdz(nr)*rv(2,nr)
                  dphiuv(1,3,ir0+nr) = dphiuv(1,3,ir0+nr)
     $             + angfac*dpdx(nr)*rv(3,nr)
                  dphiuv(2,3,ir0+nr) = dphiuv(2,3,ir0+nr)
     $             + angfac*dpdy(nr)*rv(3,nr)
                  dphiuv(3,3,ir0+nr) = dphiuv(3,3,ir0+nr)
     $             + angfac*dpdz(nr)*rv(3,nr)
                enddo
c
              endif
c
            elseif( ncplx.eq.2 )then
c
c First, compute d/dv(phi) terms need for force:
              ir = irorb
              do 450 k=1,nkloc
                ik = nk0 + k
                facr = angfac*coskr(ik,ilat)
                faci = angfac*sinkr(ik,ilat)
                do  nr=1,nrbox
                  dphix(ir+1) = dphix(ir+1) + facr*dpdx(nr)
                  dphix(ir+2) = dphix(ir+2) + faci*dpdx(nr)
                  dphiy(ir+1) = dphiy(ir+1) + facr*dpdy(nr)
                  dphiy(ir+2) = dphiy(ir+2) + faci*dpdy(nr)
                  dphiz(ir+1) = dphiz(ir+1) + facr*dpdz(nr)
                  dphiz(ir+2) = dphiz(ir+2) + faci*dpdz(nr)
                  ir = ir + 2
                enddo
  450         continue
c
c and compute d/dv(phi)*u terms needed for stress calculation
              ir = irorb
              if( ndim.eq.1 )then
                do 451 k=1,nkloc
                  ik = nk0 + k
                  facr = angfac*coskr(ik,ilat)
                  faci = angfac*sinkr(ik,ilat)
                  do  nr=1,nrbox
                    dphiuv(1,1,ir+1) = dphiuv(1,1,ir+1)
     $               + facr*dpdx(nr)*rv(1,nr)
                    dphiuv(1,1,ir+2) = dphiuv(1,1,ir+2)
     $                 + faci*dpdx(nr)*rv(1,nr)
                    ir = ir + 2
                  enddo
  451           continue
              elseif( ndim.eq.2 )then
                do 452 k=1,nkloc
                  ik = nk0 + k
                  facr = angfac*coskr(ik,ilat)
                  faci = angfac*sinkr(ik,ilat)
                  do  nr=1,nrbox
                    dphiuv(1,1,ir+1) = dphiuv(1,1,ir+1)
     $               + facr*dpdx(nr)*rv(1,nr)
                    dphiuv(1,1,ir+2) = dphiuv(1,1,ir+2)
     $                 + faci*dpdx(nr)*rv(1,nr)
                    dphiuv(2,1,ir+1) = dphiuv(2,1,ir+1)
     $               + facr*dpdy(nr)*rv(1,nr)
                    dphiuv(2,1,ir+2) = dphiuv(2,1,ir+2)
     $                 + faci*dpdy(nr)*rv(1,nr)
                    dphiuv(1,2,ir+1) = dphiuv(1,2,ir+1)
     $               + facr*dpdx(nr)*rv(2,nr)
                    dphiuv(1,2,ir+2) = dphiuv(1,2,ir+2)
     $                 + faci*dpdx(nr)*rv(2,nr)
                    dphiuv(2,2,ir+1) = dphiuv(2,2,ir+1)
     $               + facr*dpdy(nr)*rv(2,nr)
                    dphiuv(2,2,ir+2) = dphiuv(2,2,ir+2)
     $                 + faci*dpdy(nr)*rv(2,nr)
                    ir = ir + 2
                  enddo
  452           continue
              elseif( ndim.eq.3 )then
                do 453 k=1,nkloc
                  ik = nk0 + k
                  facr = angfac*coskr(ik,ilat)
                  faci = angfac*sinkr(ik,ilat)
                  do  nr=1,nrbox
                    dphiuv(1,1,ir+1) = dphiuv(1,1,ir+1)
     $               + facr*dpdx(nr)*rv(1,nr)
                    dphiuv(1,1,ir+2) = dphiuv(1,1,ir+2)
     $                 + faci*dpdx(nr)*rv(1,nr)
                    dphiuv(2,1,ir+1) = dphiuv(2,1,ir+1)
     $               + facr*dpdy(nr)*rv(1,nr)
                    dphiuv(2,1,ir+2) = dphiuv(2,1,ir+2)
     $                 + faci*dpdy(nr)*rv(1,nr)
                    dphiuv(3,1,ir+1) = dphiuv(3,1,ir+1)
     $               + facr*dpdz(nr)*rv(1,nr)
                    dphiuv(3,1,ir+2) = dphiuv(3,1,ir+2)
     $                 + faci*dpdz(nr)*rv(1,nr)
                    dphiuv(1,2,ir+1) = dphiuv(1,2,ir+1)
     $               + facr*dpdx(nr)*rv(2,nr)
                    dphiuv(1,2,ir+2) = dphiuv(1,2,ir+2)
     $                 + faci*dpdx(nr)*rv(2,nr)
                    dphiuv(2,2,ir+1) = dphiuv(2,2,ir+1)
     $               + facr*dpdy(nr)*rv(2,nr)
                    dphiuv(2,2,ir+2) = dphiuv(2,2,ir+2)
     $                 + faci*dpdy(nr)*rv(2,nr)
                    dphiuv(3,2,ir+1) = dphiuv(3,2,ir+1)
     $               + facr*dpdz(nr)*rv(2,nr)
                    dphiuv(3,2,ir+2) = dphiuv(3,2,ir+2)
     $                 + faci*dpdz(nr)*rv(2,nr)
                    dphiuv(1,3,ir+1) = dphiuv(1,3,ir+1)
     $               + facr*dpdx(nr)*rv(3,nr)
                    dphiuv(1,3,ir+2) = dphiuv(1,3,ir+2)
     $                 + faci*dpdx(nr)*rv(3,nr)
                    dphiuv(2,3,ir+1) = dphiuv(2,3,ir+1)
     $               + facr*dpdy(nr)*rv(3,nr)
                    dphiuv(2,3,ir+2) = dphiuv(2,3,ir+2)
     $                 + faci*dpdy(nr)*rv(3,nr)
                    dphiuv(3,3,ir+1) = dphiuv(3,3,ir+1)
     $               + facr*dpdz(nr)*rv(3,nr)
                    dphiuv(3,3,ir+2) = dphiuv(3,3,ir+2)
     $                 + faci*dpdz(nr)*rv(3,nr)
                    ir = ir + 2
                  enddo
  453           continue
              endif
c
            endif
c
c    Close lattice loop:
c
  900     continue
c
c    Reconstruct original coordinates:
c
      do  nr=1,nrbox
        rv(1,nr) = rv(1,nr) + rlast(1)
        rv(2,nr) = rv(2,nr) + rlast(2)
        rv(3,nr) = rv(3,nr) + rlast(3)
      enddo
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
