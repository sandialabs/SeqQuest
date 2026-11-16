c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DOTRACE
c
c
      subroutine DOTRACE( trace, norb, dmat,vmat )
c---------------------------------------------------------------
c Purpose: take trace of matrix with (complex) density matrix
c
c Written: Peter A. Schultz, 7-June-1993
c
c Revision history:
c   5Dec07-PAS/2.60c: split off entry DO0TRACE from DOTRACE
c  23Jul01-PAS/2.49: compact dmat; wtk norm convention
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  23Jul98-PAS/2.23: split basis
c  29Apr97-PAS/2.17: spherical atom trace do0trace installed
c---------------------------------------------------------------
c
c Variables:
c   trace  - output trace = Tr[dmat*vmat]
c   norb   - number of orbitals/dimension of matrices
c   dmat(i,j) - input complex density matrix
c   vmat(i,j) - matrix elements of quantity to be traced
c       (i,j real part i>j; imag part j>i)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c Input arrays:
      DIMENSION  dmat(norb,norb)
      DIMENSION  vmat(norb,norb)
c
c >>>> EXECUTABLE CODE:
c
      diag = zero
      do  i=1,norb
        diag = diag + dmat(i,i)*vmat(i,i)
      enddo
c
      offd = zero
      do  i=2,norb
        do  j=1,i-1
          offd = offd + dmat(i,j)*vmat(i,j) + dmat(j,i)*vmat(j,i)
        enddo
      enddo
c
      trace = trace + diag + offd + offd
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DO0TRACE
c
c
      subroutine DO0TRACE( trace0,  norb, natm,ntyp,nshld,
     $ itypa,numshl,lshel, occij, wtk, vmat )
c---------------------------------------------------------------
c Purpose: take trace of matrix with (complex) density matrix
c
c Written: Peter A. Schultz, 7-June-1993 (in DOTRACE)
c
c Revision history:
c   5Dec07-PAS/2.60c: split off entry DO0TRACE from DOTRACE
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c Input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  occij(nshld,nshld,ntyp)
c
      DIMENSION  vmat(norb,norb)
c
c >>>> EXECUTABLE CODE:
c
c  Compute trace of spherical-atoms = Tr[ dmat-atom0 x V]:
c
      tii = zero
      tij = zero
      i0 = 0
      do 100 iatm=1,natm
        ityp = itypa(iatm)
        nshli = numshl(ityp)
        istart = i0
        i0 = istart
        do 90 ishl=1,nshli
          lshli = lshel(ishl,ityp)
          norbshi = 2*lshli + 1
          j0 = istart
          do 80 jshl=1,ishl
            lshlj = lshel(jshl,ityp)
            norbshj = 2*lshlj + 1
            occ = occij(ishl,jshl,ityp)
            if( occ.ne.zero )then
              sum = zero
              do  iorbsh=1,norbshi
                sum = sum + occ*vmat(i0+iorbsh,j0+iorbsh)
              enddo
              if( ishl.eq.jshl )then
                tii = tii + sum
              else
                tij = tij + sum
               endif
            endif
            j0 = j0 + norbshj
   80     continue
          i0 = i0 + norbshi
   90   continue
  100 continue
c
c  I need to double off-diags as only tracing triangle
      trace0 = trace0 + wtk*(tii + tij + tij)
c
      RETURN
      END
