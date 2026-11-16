c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> XCORB
c
c
      subroutine XCORB( nang,nrd,noad,nfitd,
     $ natm,ntyp,nshld,nald, itypa, numshl,lshel,nala,ala,cala,
     $ nrad,radmsh, angpts, nfityp, radwf,
     $ gauss,rsq )
c---------------------------------------------------------------
c Purpose: calculates orbitals on radial mesh for those atom
c          types with 1-ctr potentials, for the purpose of
c          fitting of the rapid part of xc-potentials
c
c Written: Version 1.0, 12/19/89, by PJF
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted, arg list fixed
c   1Mar00-PAS/2.43: bugfix for non-xcfit
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c  20Apr93-PAS: fix underflow, ang dims, clean/opt code
c---------------------------------------------------------------
c
c      theta-phi points are chosen according to
c      A. D. McLaren's 72 point formula, which guarantees
c      accuracy through l=14. See, Math. Comp. 17, 361 (1963).
c      Radial points are given via array radmsh
c
c 20Apr93-PAS: should invert nr/nang dims in array radwf and
c              invert corresponding loops.
c              this would pay off in the d-loops below
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c data arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp)
      DIMENSION  angpts(*)
      DIMENSION  nfityp(ntyp)
c output arrays:
      DIMENSION  radwf(nang,nrd,noad,nfitd)
c scratch arrays:
      DIMENSION  gauss(*),rsq(*)
c
c local declarations:
      DIMENSION  vnorm(5)
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
c     Loop on atom types:
      do 1000 ityp=1,ntyp
c
c      skip this type if its alphas are small
c
        itypxc = nfityp(ityp)
        if( itypxc.eq.0 ) goto 1000
c
c       Loop on shells of atom type ityp
        ioa = 0
        do 800 ishl=1,numshl(ityp)
          lshl = lshel(ishl,ityp)
c
          nrmax = nrad(ityp)
          do  nr=1,nrmax
            r = radmsh(nr,ityp)
            rsq(nr) = r*r
            gauss(nr) = zero
          enddo
c
c  Compute radial factor for this radial orbital
c
          do 45 ial=1,nala(ishl,ityp)
            coeff = cala(ial,ishl,ityp)
            ali = ala(ial,ishl,ityp)
c  Underflow protection cut-off stuff:
            rsqmax = 0.4*cutexp/ali
            nrmx = nrmax
            do  nr=nrmax,1,-1
              if( rsq(nr).lt.rsqmax ) goto 39
              nrmx = nrmx - 1
            enddo
   39       continue
            do  nr=1,nrmx
              gauss(nr) = gauss(nr) + coeff* EXP( -ali*rsq(nr) )
            enddo
   45     continue
c
c  Compute angular factors for this radial orbital
c
          if( lshl.eq.0 )then
c           S-ORBITAL
            do  nr=1,nrmax
              gauss(nr) = gauss(nr)/vnorm(1)
            enddo
            do  nr=1,nrmax
              do  iang=1,nang
                radwf(iang,nr,ioa+1,itypxc) = gauss(nr)
              enddo
            enddo
c
          elseif( lshl.eq.1 )then
c           P-ORBITALS
            do  nr=1,nrmax
              gauss(nr) = gauss(nr)*radmsh(nr,ityp)/vnorm(2)
            enddo
            do  nr=1,nrmax
              iang0 = 0
              do  iang=1,nang
                radwf(iang,nr,ioa+1,itypxc) = angpts(iang0+1)*gauss(nr)
                radwf(iang,nr,ioa+2,itypxc) = angpts(iang0+2)*gauss(nr)
                radwf(iang,nr,ioa+3,itypxc) = angpts(iang0+3)*gauss(nr)
                iang0 = iang0 + 3
              enddo
            enddo
c
          elseif( lshl.eq.2 )then
c           D-ORBITALS
            do  nr=1,nrmax
              gauss(nr) = gauss(nr)*rsq(nr)/vnorm(3)
            enddo
            vn34 = vnorm(3)/vnorm(4)
            vn35 = vnorm(3)/vnorm(5)
            do  nr=1,nrmax
              iang0 = 0
              do  iang=1,nang
                x = angpts(iang0+1)
                y = angpts(iang0+2)
                z = angpts(iang0+3)
                xsq = x*x
                ysq = y*y
                zsq = z*z
                xsqysq = (xsq-ysq)*vn34
                dzsq = (zsq+zsq-xsq-ysq)*vn35
                radwf(iang,nr,ioa+1,itypxc) = x*y*gauss(nr)
                radwf(iang,nr,ioa+2,itypxc) = y*z*gauss(nr)
                radwf(iang,nr,ioa+3,itypxc) = z*x*gauss(nr)
                radwf(iang,nr,ioa+4,itypxc) = xsqysq*gauss(nr)
                radwf(iang,nr,ioa+5,itypxc) = dzsq*gauss(nr)
                iang0 = iang0 + 3
              enddo
            enddo
c
          endif
c
          ioa = ioa + 2*lshl + 1
c
c    Close shell loop:
c
  800   continue
c
c    Close atom loop:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
