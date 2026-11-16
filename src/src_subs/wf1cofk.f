c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WF1COFK
c
c
      subroutine WF1COFK( orbofk,  norb,nk, natm,ntyp,nshld,nald, nlat,
     $ itypa, numshl,lshel,nala,ala,cala, norba, ratm,rlat,coskr,sinkr )
c---------------------------------------------------------------
c Purpose: calculate quantities orbofk = k-space value of the
c     i-th orbital at the n-th nucleus.
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c  13Apr94-PAS/2.10
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c output array:
      DIMENSION  orbofk(2,natm,norb,nk)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  norba(ntyp)
      DIMENSION  ratm(3,natm),rlat(3,nlat),coskr(nk,*),sinkr(nk,*)
c
c local declarations:
      DIMENSION  ri(3),rji(3), oint(5),ph(2)
c
      DIMENSION  vnorm(5)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  ph / 1.d0,0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
      argmx = half*cutexp
c
      call MKZERO( 2*natm*norb*nk, orbofk )
c
c    Loop on remote atom site:
c
      inext = 0
      do 200 iatm=1,natm
        ityp = itypa(iatm)
        istart = inext
        inext = istart + norba(ityp)
c
c    Lattice vector loop:
c
        do 190 ilat=1,nlat
          do  ir=1,3
            ri(ir) = ratm(ir,iatm) + rlat(ir,ilat)
          enddo
c
          i0 = istart
          do 180 ishl=1,numshl(ityp)
            nali = nala(ishl,ityp)
            lshl = lshel(ishl,ityp)
            norbsh = 2*lshl + 1
            almini = ala(1,ishl,ityp)
            do  ial=1,nali
              almini = MIN(almini,ala(ial,ishl,ityp))
            enddo
c
c    Loop on site about which we compute 1-ctr density
c
            do 100 jatm=1,natm
              if( jatm.eq.iatm .and. ilat.eq.1 ) goto 100
c
              do  jr=1,3
                rji(jr) = ratm(jr,jatm) - ri(jr)
              enddo
              rjisq = rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3)
              if( almini*rjisq.gt.argmx ) goto 100
c
c      evaluate contribution of orbital i at site jatm
c
              sumg = zero
              do  ial=1,nali
                alrsq = ala(ial,ishl,ityp)*rjisq
                if( alrsq.le.argmx )
     $           sumg = sumg + cala(ial,ishl,ityp)*EXP( -alrsq )
              enddo
c
c install angular factor
c
              goto (40,41,42), lshl+1
c     s-fcn
   40         oint(1) = sumg/vnorm(1)
              goto 50
c     p-fcn:  x;  y;  z
   41         fac = sumg/vnorm(2)
              oint(1) = rji(1) * fac
              oint(2) = rji(2) * fac
              oint(3) = rji(3) * fac
              goto 50
c     d-fcn:  xy;  yz;  zx;  x2-y2;  2z2-x2-y2
   42         fac = sumg/vnorm(3)
              oint(1) = rji(1)*rji(2) * fac
              oint(2) = rji(2)*rji(3) * fac
              oint(3) = rji(3)*rji(1) * fac
              fac = sumg/vnorm(4)
              oint(4) = (rji(1)**2-rji(2)**2) * fac
              fac = sumg/vnorm(5)
              oint(5) = (two*rji(3)**2-rji(1)**2-rji(2)**2) * fac
c
   50         continue
c
c    Alpha loop:
c
              i = i0
              do  is=1,norbsh
                i = i + 1
                do  k=1,nk
                  ph(1) = coskr(k,ilat)
                  ph(2) = sinkr(k,ilat)
                  orbofk(1,jatm,i,k)=orbofk(1,jatm,i,k) + oint(is)*ph(1)
                  orbofk(2,jatm,i,k)=orbofk(2,jatm,i,k) + oint(is)*ph(2)
                enddo
              enddo
c
  100       continue
c
            i0 = i0 + norbsh
  180     continue
  190   continue
  200 continue
c
c    That's all Folks!
c
      RETURN
      END
