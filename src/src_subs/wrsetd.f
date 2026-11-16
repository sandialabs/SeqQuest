c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRSETD
c
c
      subroutine WRSETD( iset0fl, iprgm,jprgm, nxcd,
     $ norb,nk, natm,ntyp,nrd, noad, nlat,nlat1c,nlat2c,nlat3c,
     $ atm0engy, xc0val,  vatrad,vesrad,vxcrad,excrad, gntcomb,
     $ rhoatom,nearest, near2c, s1atom, orbint,orbofk,
     $ nang,nfitd, angylm,nfityp,radwf )
c---------------------------------------------------------------
c Purpose: write out (almost) generally dimensioned iset0fl file
c
c Written: Peter A. Schultz, 17-February-1994
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c  21Jun01-PAS/2.48: replace STOPs, writes
c   3May99-PAS/2.35: use near2c to determine 2c defect use
c  17Dec98-PAS/2.30: take all force stuff out of setup phase
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  22May98-PAS/2.22: bugfix
c  28Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  xc0val(nxcd)
      DIMENSION  vatrad(nrd,ntyp),vesrad(nrd,ntyp),
     $           vxcrad(nrd,ntyp),excrad(nrd,ntyp)
      DIMENSION  gntcomb(8100)
      DIMENSION  rhoatom(natm,3),nearest(3,natm)
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  orbint(2,norb,nk),orbofk(2,natm,norb,nk)
c for stress routines:
      DIMENSION  angylm(nang,9),nfityp(ntyp),radwf(nang,nrd,noad,nfitd)
c
c  Dump to iset0fl:
c
      REWIND( unit=iset0fl )
      write(iset0fl)  iprgm,jprgm
      write(iset0fl)  nxcd, norb,nk, natm,ntyp,nrd, noad, near2c
      write(iset0fl)  nlat, nlat1c,nlat2c,nlat3c
      write(iset0fl)  atm0engy,xc0val
      write(iset0fl)  vatrad,vesrad,vxcrad,excrad,gntcomb
      if( near2c.eq.2 )then
        write(iset0fl)  rhoatom,nearest,s1atom,orbint,orbofk
      else
        write(iset0fl)  rhoatom,nearest,s1atom
      endif
      write(iset0fl)  nang,nfitd
      write(iset0fl)  angylm,nfityp,radwf
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDSETD
c
c
      ENTRY RDSETD(iset0fl, iprgm,jprgm, nxcd,
     $ norb,nk, natm,ntyp,nrd, noad, nlat,nlat1c,nlat2c,nlat3c,
     $ atm0engy, xc0val,  vatrad,vesrad,vxcrad,excrad, gntcomb,
     $ rhoatom,nearest, near2c, s1atom, orbint,orbofk,
     $ nang,nfitd, angylm,nfityp,radwf )
c
c  Retrieve from file:
c
      REWIND( unit=iset0fl )
      read(iset0fl)  iprgm,jprgm
      if( iprgm.lt.230 ) call STOPXERR( 'isetdvrs/restart version' )
      read(iset0fl)  mxcd, morb,mk, matm,mtyp,mrd, moad, mear2c
c check the dimensioning:
      if( morb.ne.norb .or. mk.ne.nk .or. matm.ne.natm .or.
     $    mtyp.ne.ntyp .or. mrd.ne.nrd .or. moad.ne.noad )then
        call FLGETIWR( IERRFL )
        write(IERRFL,*) 'DATA:  norb, nk, natm, ntyp, nrd, noad'
        write(IERRFL,*) 'input:',norb,nk, natm,ntyp,nrd, noad
        write(IERRFL,*) 'saved:',morb,mk, matm,mtyp,mrd, moad
        call STOPXERR( 'isetdims/restart file dim mismatch' )
      endif
      read(iset0fl)  nlat, nlat1c,nlat2c,nlat3c
      read(iset0fl)  atm0engy,(xc0val(i),i=1,mxcd)
      read(iset0fl)  vatrad,vesrad,vxcrad,excrad, gntcomb
      if( mear2c.ne.near2c )then
        call STOPXERR( 'near2c  /input near2c .ne. disk near2c' )
      elseif( near2c.eq.2 )then
        read(iset0fl)  rhoatom,nearest,s1atom,orbint,orbofk
      else
        read(iset0fl)  rhoatom,nearest,s1atom
      endif
      read(iset0fl)  mang,mfitd
      if( mang.ne.nang .or. mfitd.ne.nfitd ) call STOPXERR( 'isetd1c' )
      read(iset0fl)  angylm,nfityp,radwf
c
c    That's all Folks!
c
      RETURN
      END
