c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHO1CTR
c
c
      subroutine RHO1CTR( IWR, irh1orb,icpctfl,icountr,
     $ norb,nk, natm,ntyp,nshld,nald, noad,
     $ itypa, numshl,lshel,nala,ala,cala,  norba,ioffa, occij,
     $ lmx1ctr,almnv1c, nsyms,dsym,naofsym,
     $ hh,vk, n1r,n2r,n3r,nptr, weight, mngr3,mxgr3,ngr,
     $ dcountr,rh1coef, dmatwk,
     $ rho, dmatofk, r,cnrho )
c---------------------------------------------------------------
c Purpose: Computes and subtracts from grid density rho()
c   rapidly-varying "1-ctr density" due to delta density matrix
c   whose electrostatic potential has already been done in
c   VLOCxxx routines, and must not be repeated in FFT.
c
c Written: P.J. Feibelman
c
c Revision history:
c  19Dec99-PAS/2.42: clean-up; read/write big grid records
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  29Jul98-PAS/2.23: split basis occs
c   1May97-PAS/2.17: naming conventions, new grid coord rtn
c  15Dec95-PAS/2.14: shell-by-shell grid ranges, cleanup
c  29Aug94-PAS/2.11: cleanup, rationalize args.
c   6Apr90-PAS: make more efficient; use less disk
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  norba(ntyp),ioffa(natm), occij(nshld,nshld,ntyp)
      DIMENSION  hh(3,3),vk(3,nk), dsym(5,5,3,*),naofsym(natm,*)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  dcountr(10,natm),rh1coef(10,noad,natm)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
      DIMENSION  dmatofk(2,norb,norb,nk)
c input/output:
      DIMENSION  rho(nptr)
c scratch arrays:
      DIMENSION  cnrho(nptr),r(3,nptr)
      DIMENSION  dmatwk(noad,noad)
c local declarations:
      DIMENSION  r0(3),rij(3)
      DIMENSION  number(3)
      DATA  number / 1,4,10 /
      DATA  zero,two / 0.d0,2.d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
      r0(1) = zero
      r0(2) = zero
      r0(3) = zero
      call MKRBOX( hh,r0, n1r,n2r,n3r, r )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              One atom contributions
c * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      i1next = 1
      do 100 iatm=1,natm
        ityp = itypa(iatm)
        nshl = numshl(ityp)
c I need this for grid indexing arrays:
        iorb1 = i1next
        do  ishl=1,nshl
          i1next = i1next + 2*lshel(ishl,ityp) + 1
        enddo
        if( lmx1ctr(ityp).lt.0 ) goto 100
c
        al1cmin = almnv1c(ityp)
c
        do  jr=1,3
          rij(jr) = zero
        enddo
c
c        Sum density matrix over surface Brillouin zone
c
        call MKZERO( noad**2, dmatwk )
        call ROTGRP( dmatwk, iatm,iatm,rij,
     $   norb,nk, natm,ntyp,nshld, noad,
     $   itypa, numshl,lshel, ioffa, vk,rij,
     $   nsyms,dsym,naofsym,
     $   dmatofk,cnrho )
c -->    dmat-i  cphi-s
c
        iorbat0 = 0
        do 85 ishl=1,nshl
          lshli = lshel(ishl,ityp)
          li0 = lshli*lshli
          nali = nala(ishl,ityp)
          norbshi = 2*lshli + 1
c
          nrshi = ngr(iorb1+iorbat0)
          imin3sh = mngr3(iorb1+iorbat0)
c
          alimax = ala(nali,ishl,ityp)
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.gt.alimax )  alimax = ali
          enddo
c
          jorbat0 = 0
          do 80 jshl=1,ishl
            lshlj = lshel(jshl,ityp)
            lj0 = lshlj*lshlj
            nalj = nala(jshl,ityp)
            norbshj = 2*lshlj + 1
            lsum = lshli + lshlj
c
            nrshj = ngr(iorb1+jorbat0)
            jmin3sh = mngr3(iorb1+jorbat0)
c
            nrs = MIN( nrshi,nrshj )
            ijmin3 = MAX( imin3sh,jmin3sh )
            nr0 = n1r*n2r*(ijmin3-1)
c
            aljmax = ala(nalj,jshl,ityp)
            do  jal=1,nalj
              alj = ala(jal,jshl,ityp)
              if( alj.gt.aljmax )  aljmax = alj
            enddo
            alsmx = alimax + aljmax
c
            if( lsum .gt.lmx1ctr(ityp) .or. alsmx.lt.al1cmin ) goto 80
c
            norbj0 = norbshj
            joff = 0
            if( jshl.eq.ishl )then
              norbj0 = 0
              joff = 1
            endif
c
            occ = occij(ishl,jshl,ityp)
c
            li = li0
            norbshj = norbj0
            do 75 iorbsh=1,norbshi
              iorbat = iorbat0 + iorbsh
              li = li + 1
c
              lj = lj0
              norbshj = norbshj + joff
              do 70 jorbsh=1,norbshj
                jorbat = jorbat0 + jorbsh
                lj = lj + 1
                if( iorbat.eq.jorbat )then
                  dmatr = dmatwk(iorbat,jorbat) - occ
                elseif( li.eq.lj )then
                  dmatr = two*( dmatwk(iorbat,jorbat) - occ )
                else
                  dmatr = two*dmatwk(iorbat,jorbat)
                endif
c
C               read(irh1orb)  (cnrho(nr),nr=1,nrs)
                call READBIG( irh1orb, nrs, cnrho )
                do  nr=1,nrs
                  rho(nr0+nr) = rho(nr0+nr) - dmatr*cnrho(nr)
                enddo
c
c    Close orbital loops:
c
   70         continue
   75       continue
c
c    Close radial function loops:
c
   80     jorbat0 = jorbat0 + norbshj
   85   iorbat0 = iorbat0 + norbshi
c
c    Close center loops:
c
  100 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                 two atom contributions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      i1next = 1
      do 200 iatm=1,natm
        ityp = itypa(iatm)
        nshl = numshl(ityp)
        iorb1 = i1next
        do  ishl=1,nshl
          i1next = i1next + 2*lshel(ishl,ityp) + 1
        enddo
        if( lmx1ctr(ityp).lt.0 ) goto 200
c
        iorbat0 = 0
        do 180 ishl=1,nshl
          lshl = lshel(ishl,ityp)
          nali = nala(ishl,ityp)
          norbsh = 2*lshl + 1
c
          nrs = ngr(iorb1+iorbat0)
          imin3 = mngr3(iorb1+iorbat0)
          nr0 = n1r*n2r*(imin3-1)
c
          iorbat = iorbat0
          iorbat0 = iorbat0 + norbsh
c
          alimax = ala(nali,ishl,ityp)
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.gt.alimax ) alimax = ali
          enddo
          if( lshl.gt.lmx1ctr(ityp).or.alimax.lt.almnv1c(ityp) )goto 180
c
          do 170 iorbsh=1,norbsh
            iorbat = iorbat + 1
c
c      orbofk contribution; read i-th orbital
c
C           read(icpctfl)  (cnrho(nr),nr=1,nrs)
            call READBIG( icpctfl, nrs, cnrho )
c
c      compute contributions to delta rho ...
c
c       ... from orbitals at neighboring sites
c
            dmijwk = two*rh1coef(1,iorbat,iatm)
            do  nr=1,nrs
              rho(nr0+nr) = rho(nr0+nr) - dmijwk*cnrho(nr)
            enddo
c
            if( lshl+1 .gt. lmx1ctr(ityp) ) goto 170
c
c       ... from grad orbital j at siteiatm
c
            do 150 lgr=1,3
C             read(icpctfl)  (cnrho(nr),nr=1,nrs)
              call READBIG( icpctfl, nrs, cnrho )
              dmijwk = two*rh1coef(1+lgr,iorbat,iatm)
              do  nr=1,nrs
                rho(nr0+nr) = rho(nr0+nr) - dmijwk*cnrho(nr)
              enddo
  150       continue
c
            if( lshl+2 .gt. lmx1ctr(ityp) ) goto 170
c
c       ... from grad(grad orbital j) at site i
c
            lcoef = 4
            do 162 lgr=1,3
              do 161 llgr=1,lgr
                lcoef = lcoef + 1
C               read(icpctfl)  (cnrho(nr),nr=1,nrs)
                call READBIG( icpctfl, nrs, cnrho )
                dmijwk = two*rh1coef(lcoef,iorbat,iatm)
                do  nr=1,nrs
                  rho(nr0+nr) = rho(nr0+nr) - dmijwk*cnrho(nr)
                enddo
  161         continue
  162       continue
c
c    Close  orbital loop:
c
  170     continue
c
c    Close radial function loop:
c
  180   continue
c
c    Close center loop:
c
  200 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Add back counter-charge associated with center i that
c     cancels the l=0,1,2 multipole moments of the 1-ctr charge
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( lstout.gt.0 ) write(IWR,*) 'RHO1CTR:'
      iorb = 1
      do 300 iatm=1,natm
        ityp = itypa(iatm)
        nrat = ngr(iorb)
        imin3 = mngr3(iorb)
        do  ishl=1,numshl(ityp)
          if( ngr(iorb).gt.nrat )then
            nrat = ngr(iorb)
            imin3 = mngr3(iorb)
          endif
          iorb = iorb + 2*lshel(ishl,ityp) + 1
        enddo
        if( lmx1ctr(ityp).lt.0 ) goto 300
c
        nr0 = n1r*n2r*(imin3-1)
        norbsh = number(lmx1ctr(ityp)+1)
        if( lstout.gt.0 )
     $  write(IWR,9300) iatm,norbsh,(dcountr(li,iatm),li=1,norbsh)
 9300   format(1x,'iatm=',i3,'; li=1-',i2,' dcountr values:'/(1p5e15.8))
c
        do 280 iorbsh=1,norbsh
C         read(icountr)  (cnrho(nr),nr=1,nrat)
          call READBIG( icountr, nrat, cnrho )
c
          do  nr=1,nrat
            rho(nr0+nr) = rho(nr0+nr) + dcountr(iorbsh,iatm)*cnrho(nr)
          enddo
c
  280   continue
c
  300 continue
c
c    That's all Folks!
c
      RETURN
      END
