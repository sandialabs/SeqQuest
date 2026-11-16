c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DOSDUMP
c
c
      subroutine DOSDUMP( ipostfl,
     $ nspin, norb,nstate,nk,ncplx, natm,ntyp,nshld,nald,
     $ itypa,numshl,lshel,
     $ ala,cala,nala,
     $ temp, efermi, eigval,eigpop,wtk,veck, typnm )
c---------------------------------------------------------------
c Purpose: write out data necessary to compute DOS and properly
c          identify the data
c
c Written: 30-May-2001, for v2.47
c
c Revision history:
c   9Sep14-PAS/2.66: cosmetics (units/params)
c   3May12-ACP&AHE/2.63:  output ityp, 11Oct12-PAS: fix format
c  29May10-AHE&PAS/2.61k: fix formating problems
c   9Apr07-AHE&PAS/2.59c: new output for dosdump put into master
c  15Sep05-AHE/2.58: placed contracted basis set information for
c                    matrix element calculations
c  11Dec01-PAS/2.51: add spin, pop, temp output; self-doc labels
c  31Oct01-PAS/2.50: naming convention change
c  30May01-PAS/2.47: modified Art Edwards code for output
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp), lshel(nshld,ntyp)
      DIMENSION nala(nshld,ntyp),ala(nald,nshld,ntyp),
     $     cala(nald,nshld,ntyp)
      DIMENSION  efermi(*)
      DIMENSION  eigval(*),eigpop(*)
      DIMENSION  wtk(nk),veck(3,nk)
      CHARACTER*(*)  typnm(ntyp)
c
c local declarations:
c     basis function labels
      CHARACTER*6  orblbl(9)
      DATA  orblbl / 's     ','x     ','y     ','z     ',
     $      'xy    ','yz    ','zx    ','x2-y2 ','3z2-r2' /
c
c >>>> EXECUTABLE CODE:
c
 9000 format(a)
c Output the parameters that define the size of everything:
      write(ipostfl,9000) 'norb,ntyp,nstate,nk,ncplx,nspin'
      write(ipostfl,9011)  norb,ntyp,nstate,nk,ncplx,nspin
 9011 format( 6(1x,i6) )
c
c Output fermi level(s) and temperature for occupying states:
      if( nspin.eq.1 )then
        write(ipostfl,9000) 'fermi level and electron temperature (Ry)'
      else
        write(ipostfl,9000) 'fermi levels (up/dn) and temperature (Ry)'
      endif
      write(ipostfl,9021)  (efermi(ispin),ispin=1,nspin), temp
 9021 format( 5(1x,f14.8) )
 9023 format( 3(1x,f14.8) )
c
c Output weights for k-vectors ...
      write(ipostfl,9000) 'weights for k-points'
      write(ipostfl,9021)  (wtk(k), k=1,nk)
c
c Output k-vectors (in cartesian coordinates)
      write(ipostfl,9000) 'k-points (bohr^-1)'
      write(ipostfl,9023)  ((veck(ic,ik),ic=1,3),ik=1,nk)
c  ... and output the eigenvalues:
      nksp = nstate*nk*nspin
      write(ipostfl,9000)
     $  'eigval (((i,k,isp),i=1,nstate),k=1,nk),isp=1,nspin (Ry)'
      write(ipostfl,9021) (eigval(iksp),iksp=1,nksp)
c  ... and output the level occupations:
      write(ipostfl,9000)
     $  'eigpop (((i,k,isp),i=1,nstate),k=1,nk),isp=1,nspin'
      write(ipostfl,9021) (eigpop(iksp),iksp=1,nksp)
c
c Output of orbital id information:
      write(ipostfl,9000) 'labels for basis functions and atoms'
      write(ipostfl,9041)
 9041 format(1x,'  bf# ',1x,' atm#',1x,'bftype',1x,'atomname',1x,
     &     'atomtype')
 9042 format(1x,   i6,   1x,   i5,  1x,   a6,   1x, a8, 1x, i2)
      iorb = 0
      do  iatm=1,natm
        ityp = itypa(iatm)
        nshl = numshl(ityp)
        do  ishl=1,nshl
          lshl = lshel(ishl,ityp)
          norbat = 2*lshl + 1
          lmorb = lshl**2
          do  iorbat=1,norbat
            iorb = iorb + 1
            lmorb = lmorb + 1
cACP            write(ipostfl,9042)  iorb,iatm,orblbl(lmorb),typnm(ityp)
            write(ipostfl,9042)iorb,iatm,orblbl(lmorb),typnm(ityp),ityp
          enddo
        enddo
      enddo
c
c Output the (segmented) internal basis set contractions:
      write(ipostfl,9000) 'contracted basis fcns for each atom type'
      write(ipostfl,9043)  (numshl(ityp),ityp=1,ntyp)
      do  ityp=1,ntyp
         nshl = numshl(ityp)
         write(ipostfl,9043)  (nala(ishl,ityp),ishl=1,nshl)
         write(ipostfl,9043)  (lshel(ishl,ityp),ishl=1,nshl)
         do  ishl=1,nshl
            nal = nala(ishl,ityp)
            write(ipostfl,9044)  (ala(ial,ishl,ityp),ial=1,nal)
            write(ipostfl,9044)  (cala(ial,ishl,ityp),ial=1,nal)
         enddo
      enddo
 9043 format(5(2x,i2))
 9044 format(5(2x,f15.8))
c
c    That's all Folks!
c
      RETURN
      END
