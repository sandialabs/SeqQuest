c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CMPTORB
c
c
      subroutine CMPTORB( icpctfl,
     $ natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat,maxwf,
     $ itypa, numshl,lshel,nala,ala,cala, ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ r, rv,rsq,gsum, wf )
c---------------------------------------------------------------
c Purpose: compute compact orbitals on mesh for two-atom terms
c          in RHO1CTR.
c
c Written: Peter A. Schultz, based on PJFs original routine
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted
c  21Jun01-PAS/2.48: replace STOPs
c  19Dec99-PAS/2.42: clean-up; read/write big grid records
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   1May97-PAS/2.17: naming conventions, on-off criterion
c  15Dec95-PAS/2.14: shell-by-shell ranges, cleanup
c  17Sep93-PAS: more stringent cutoffs installed
c   5Apr90-PAS: to speed up and bring in line with ORBBULK
c---------------------------------------------------------------
c
c  In order to give a "1-center density" that has the symmetry
c  of the crystal, these bloch orbitals all need to correspond
c  to k(parallel)=0.
c
c  On input, coordinate mesh assumed to be stored in r()
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  ratm(3,natm),hh(3,3),rlat(3,nlat)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
c input/output array:
      DIMENSION  r(3,nptr)
c scratch arrays:
      DIMENSION  rv(3,nptr),rsq(nptr),gsum(nptr)
      DIMENSION  wf(*)
c local declarations:
      DIMENSION  ratom(3)
c
      DIMENSION  vnorm(5)
      DIMENSION  lwrite(3)
      DATA  lwrite / 10,12,5 /
      DATA  zero,half,three / 0.d0,0.5d0,3.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
      argmx = half*cutgaus
c
c    Center loop:
c
      i1next = 1
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
        nshl = numshl(ityp)
        iorb1 = i1next
        do  ishl=1,nshl
          i1next = i1next + 2*lshel(ishl,ityp) + 1
        enddo
        if( lmx1ctr(ityp).lt.0 ) goto 1000
c
        almin = ala(1,1,ityp)
        almax = ala(1,1,ityp)
        do  ishl=1,nshl
          do  ial=1,nala(ishl,ityp)
            if( ala(ial,ishl,ityp).lt.almin ) almin = ala(ial,ishl,ityp)
            if( ala(ial,ishl,ityp).gt.almax ) almax = ala(ial,ishl,ityp)
          enddo
        enddo
c
        al1c = almax
        do  ishl=1,nshl
          do  ial=1,nala(ishl,ityp)
            ali = ala(ial,ishl,ityp)
            if( ali.ge.almnv1c(ityp) .and. ali.lt.al1c ) al1c = ali
          enddo
        enddo
        twoalc = al1c + almin
        rsqtop = cutgrid/twoalc
c
        nlatm = nlat
        iorbat0 = 0
        do 800 ishl=1,nshl
          lshl = lshel(ishl,ityp)
          nali = nala(ishl,ityp)
          norbsh = 2*lshl + 1
c
          nrs = ngr(iorb1+iorbat0)
          imin3 = mngr3(iorb1+iorbat0)
          nr0 = n1r*n2r*(imin3-1)
c
          iorbat0 = iorbat0 + norbsh
c
          alimax = ala(nali,ishl,ityp)
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.gt.alimax )  alimax = ali
          enddo
          if( lshl.gt.lmx1ctr(ityp).or.alimax.lt.almnv1c(ityp) )goto 800
c
          al1c = almax
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.ge.almnv1c(ityp) .and. ali.lt.al1c ) al1c = ali
          enddo
          twoalc = al1c + almin
          r2cutsh = cutgrid/twoalc
c
          nrset = lwrite(lshl+1)*nrs
          if( nrset.gt.2*maxwf )  call STOPXERR( 'cmptorb1' )
c
          do  nr=1,nrset
            wf(nr) = zero
          enddo
c
          nwrite = 0
          nlatmx = 1
          do 900 ilat=1,nlatm
c
            do  jr=1,3
              ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
            enddo
c
            do  nr=1,nrs
              rv(1,nr) = r(1,nr0+nr) - ratom(1)
              rv(2,nr) = r(2,nr0+nr) - ratom(2)
              rv(3,nr) = r(3,nr0+nr) - ratom(3)
              rsq(nr) = rv(1,nr)**2 + rv(2,nr)**2 + rv(3,nr)**2
              gsum(nr) = zero
            enddo
            rsqmin = rsq(1)
            rsqmax = rsq(1)
            do  nr=1,nrs
              if( rsq(nr).lt.rsqmin ) rsqmin = rsq(nr)
              if( rsq(nr).gt.rsqmax ) rsqmax = rsq(nr)
            enddo
            if( rsqmin.lt.rsqtop .and. ilat.gt.nlatmx ) nlatmx = ilat
            if( rsqmin.gt.r2cutsh ) goto 900
c
            do 290 ial=1,nali
              alf = ala(ial,ishl,ityp)
              if( alf.lt.almnv1c(ityp) ) goto 290
c
              twoalf = alf + alf
              twoalc = alf + almin
              rsqcut = cutgrid/twoalc
              if( rsqmin.gt.rsqcut ) goto 290
c
              coef = cala(ial,ishl,ityp)/vnorm(lshl+1)
              if( twoalf*rsqmax.lt.argmx )then
                do  nr=1,nrs
                  gsum(nr) = gsum(nr) + coef*EXP( -alf*rsq(nr) )
                enddo
              else
                do  nr=1,nrs
                  if( rsq(nr).lt.rsqcut )
     $            gsum(nr) = gsum(nr) + coef*EXP( -alf*rsq(nr) )
                enddo
              endif
  290       continue
c
            if( lshl.eq.0 )then
c
              if( lmx1ctr(ityp).eq.0 )then
c
                nwrite = 1
                do  nr=1,nrs
                  wf(      nr) = wf(      nr) + gsum(nr)
                enddo
c
              elseif( lmx1ctr(ityp).eq.1 )then
c
                nwrite = 4
                do  nr=1,nrs
                  wf(      nr) = wf(      nr) + gsum(nr)
                  wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*rv(1,nr)
                  wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*rv(2,nr)
                  wf(3*nrs+nr) = wf(3*nrs+nr) + gsum(nr)*rv(3,nr)
                enddo
c
              else
c
                nwrite = 10
                do  nr=1,nrs
                  wf(      nr) = wf(      nr) + gsum(nr)
                  wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*rv(1,nr)
                  wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*rv(2,nr)
                  wf(3*nrs+nr) = wf(3*nrs+nr) + gsum(nr)*rv(3,nr)
                  wf(4*nrs+nr) = wf(4*nrs+nr) + gsum(nr)*
     $             rv(1,nr)*rv(1,nr)
                  wf(5*nrs+nr) = wf(5*nrs+nr) + gsum(nr)*
     $             rv(2,nr)*rv(1,nr)
                  wf(6*nrs+nr) = wf(6*nrs+nr) + gsum(nr)*
     $             rv(2,nr)*rv(2,nr)
                  wf(7*nrs+nr) = wf(7*nrs+nr) + gsum(nr)*
     $             rv(3,nr)*rv(1,nr)
                  wf(8*nrs+nr) = wf(8*nrs+nr) + gsum(nr)*
     $             rv(3,nr)*rv(2,nr)
                  wf(9*nrs+nr) = wf(9*nrs+nr) + gsum(nr)*
     $             rv(3,nr)*rv(3,nr)
                enddo
              endif
c
            elseif( lshl.eq.1 .and. lshl.lt.lmx1ctr(ityp) )then
c
              nwrite = 12
              do  nr=1,nrs
                wf(       nr) = wf(       nr) + gsum(nr)*rv(1,nr)
                wf(   nrs+nr) = wf(   nrs+nr) + gsum(nr)*rv(1,nr)*
     $           rv(1,nr)
                wf( 2*nrs+nr) = wf( 2*nrs+nr) + gsum(nr)*rv(1,nr)*
     $           rv(2,nr)
                wf( 3*nrs+nr) = wf( 3*nrs+nr) + gsum(nr)*rv(1,nr)*
     $           rv(3,nr)
                wf( 4*nrs+nr) = wf( 4*nrs+nr) + gsum(nr)*rv(2,nr)
                wf( 5*nrs+nr) = wf( 5*nrs+nr) + gsum(nr)*rv(2,nr)*
     $           rv(1,nr)
                wf( 6*nrs+nr) = wf( 6*nrs+nr) + gsum(nr)*rv(2,nr)*
     $           rv(2,nr)
                wf( 7*nrs+nr) = wf( 7*nrs+nr) + gsum(nr)*rv(2,nr)*
     $           rv(3,nr)
                wf( 8*nrs+nr) = wf( 8*nrs+nr) + gsum(nr)*rv(3,nr)
                wf( 9*nrs+nr) = wf( 9*nrs+nr) + gsum(nr)*rv(3,nr)*
     $           rv(1,nr)
                wf(10*nrs+nr) = wf(10*nrs+nr) + gsum(nr)*rv(3,nr)*
     $           rv(2,nr)
                wf(11*nrs+nr) = wf(11*nrs+nr) + gsum(nr)*rv(3,nr)*
     $           rv(3,nr)
              enddo
c
            elseif( lshl.eq.1 .and. lshl.ge.lmx1ctr(ityp) )then
c
              nwrite = 3
              do  nr=1,nrs
                wf(      nr) = wf(      nr) + gsum(nr)*rv(1,nr)
                wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*rv(2,nr)
                wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*rv(3,nr)
              enddo
c
            elseif( lshl.eq.2 )then
c
              if( lshl.lt.lmx1ctr(ityp) ) call STOPXERR( 'cmptorb2' )
              ang4 = vnorm(3)/vnorm(4)
              ang5 = vnorm(3)/vnorm(5)
c
              nwrite = 5
              do  nr=1,nrs
                wf(      nr) = wf(      nr) + gsum(nr)*
     $           rv(1,nr)*rv(2,nr)
                wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*
     $           rv(2,nr)*rv(3,nr)
                wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*
     $           rv(3,nr)*rv(1,nr)
                wf(3*nrs+nr) = wf(3*nrs+nr) + gsum(nr)*
     $           (rv(1,nr)**2 - rv(2,nr)**2)*ang4
                wf(4*nrs+nr) = wf(4*nrs+nr) + gsum(nr)*
     $           (three*rv(3,nr)**2 - rsq(nr))*ang5
              enddo
c
            else
              call STOPXERR( 'cmptorb3' )
            endif
c
c    Close lattice vector loop:
c
  900     continue
          nlatm = nlatmx
c
c    Save short-range orbitals to disk:
c
          if( nwrite.gt.0 )then
            nrset0 = 0
            do  nwr=1,nwrite
C             write(icpctfl)  (wf(nr),nr=nrset0+1,nrset0+nrs)
              call WRITBIG( icpctfl, nrs, wf(nrset0+1) )
              nrset0 = nrset0 + nrs
            enddo
          endif
c
c    Close radial function, then center loop:
c
  800   continue
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
