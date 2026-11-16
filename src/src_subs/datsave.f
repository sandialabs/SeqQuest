c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DATSAVE
c
c
      subroutine DATSAVE( norb, natm,ntyp,nshld,nald,noad,
     $ itypa, numshl,lshel,nala,ala,cala, znuc,occ,
     $ ratm, rprim,nlat,rlat, hh, n1r,n2r,n3r,
     $ orig, rscale,
     $ lshl2, nal2,al2,cal2 )
c---------------------------------------------------------------
c Purpose: write to disk file geometry and basis set data
c          so as to (1) give imp code something it can check
c          against, and (2) give density fitting routine info
c          necessary to attempt fit of output rho with
c          overlapping spherical atoms
c
c Written: Peter A. Schultz, 9-December-1989
c
c Revision History:
c  26Jan09-PAS/2.62j: convert rhosave into subroutine
c   6Dec01-PAS/2.51: file mgmt
c   2Jul01-PAS/2.48: replace writes
c   7Jul00-PAS/2.45: unified scaling
c  19Dec99-PAS/2.42: cleanup; read/write big grid records
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  22May98-PAS/2.22: bugfix
c   5Mar98-PAS/2.21: give data for charge plotting
c  24Apr97-PAS/2.17: naming conventions
c  15Dec95-PAS/2.14: shell-by-shell grid ranges
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  znuc(ntyp), occ(noad,ntyp)
      DIMENSION  ratm(3,natm), rprim(3,3), rlat(3,nlat), hh(3,3)
      DIMENSION  orig(3),rscale(3)
c Scratch arrays:
      DIMENSION  lshl2(nshld*ntyp),nal2(nshld*ntyp)
      DIMENSION  al2(nald*nshld*ntyp),cal2(nald*nshld*ntyp)
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
C      OPEN(unit=idatsav,file='lcao.data',status='unknown',
C     $ form='unformatted',err=1300)
      call FLOPENB( idatsav, 'data' )
c
c  Compact the arrays with big holes in 'em:
c
      nrf = 0
      nal = 0
      do 1000 ityp=1,ntyp
c
        do 100 ishl=1,numshl(ityp)
          nrf = nrf + 1
          lshl2(nrf) = lshel(ishl,ityp)
          nal2(nrf) = nala(ishl,ityp)
c
          do 10 ial=1,nala(ishl,ityp)
            nal = nal + 1
            cal2(nal) = cala(ial,ishl,ityp)
            al2(nal) = ala(ial,ishl,ityp)
   10     continue
  100   continue
 1000 continue
c
c  Write out the good stuff:
c
      write(idatsav)  natm,ntyp,nrf,nal,n1r,n2r,n3r,nlat,noad
      write(idatsav)  itypa,numshl,occ,znuc
      write(idatsav)  (lshl2(irf),irf=1,nrf),(nal2(irf),irf=1,nrf)
      write(idatsav)  (cal2(ial),ial=1,nal),(al2(ial),ial=1,nal)
      write(idatsav)  hh,rprim,ratm,rlat,orig,rscale
c
      call FLCLOSE( idatsav )
c
      write(IWR,*) '>>>>> geom/basis data saved to file'
c
      RETURN
c
c >>>> ERROR HANDLING:
c
 1300 write(IWR,*)  ' '
      write(IWR,*)  '*****  ERROR  *****'
      write(IWR,*)  '>>>>> cannot open unit #',idatsav
      write(IWR,*)  '>>>>> geometry/basis data NOT saved'
      write(IWR,*)  ' '
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOSAVE
c
c
      subroutine RHOSAVE( nptr, delrho )
c---------------------------------------------------------------
c Purpose: write density to disk file
c
c Written: Peter A. Schultz, 9-December-1989
c
c Revision History:
c  26Jan09-PAS/2.62j: converted rhosave into subroutine
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  delrho(*)
c
c >>>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
C      OPEN(unit=idatsav,file='lcao.data',status='old',
C     $ form='unformatted', err=1310 )
      call FLOPENB( idatsav, 'data' )
c
c read through header data:
      read(idatsav)   dum
      read(idatsav)  idum
      read(idatsav)  idum
      read(idatsav)   dum
      read(idatsav)   dum
c output the delta-density:
c     write(idatsav)  (delrho(i),i=1,nptr)
      call WRITBIG( idatsav, nptr, delrho )
c
      call FLCLOSE( idatsav )
c
      write(IWR,*) '>>>>> grid delta-density saved to file'
c
      RETURN
c
 1310 write(IWR,*)  ' '
      write(IWR,*)  '*****  ERROR  *****'
      write(IWR,*)  '>>>>> cannot open unit #',idatsav
      write(IWR,*)  '>>>>> geometry/basis data NOT saved'
      write(IWR,*)  ' '
c
      RETURN
      END
