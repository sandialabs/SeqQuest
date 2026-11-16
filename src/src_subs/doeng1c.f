c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DOENG1C
c
c
      subroutine DOENG1C( ves1ctr,ves1c0,
     $ nk,norb,natm,ntyp,nshld, itypa,numshl,lshel, occij,wtk,
     $ ves, dmat )
c---------------------------------------------------------------
c Purpose: get 1-center es energies
c
c Written: Peter A. Schultz, 20-May-2001 for v2.47
c
c Revision history:
c  24Jul01-PAS/2.49: packed dmat/cpas=temp workaround
c  20May01-PAS/2.47: extracted from DOENGY
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  occij(nshld,nshld,ntyp), wtk(nk)
cpas  DIMENSION  ves(norb,norb,nk), dmat(norb,norb,nk)
      DIMENSION  ves(norb,norb,nk), dmat(2,norb,norb,nk)
c
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      ves1ctr = zero
      ves1c0  = zero
      do  k=1,nk
c
c  Non-spherical electrostatic energy (from 1-center stuff)
cpas        call DOTRACE( ves1ctr, norb, dmat(1,1,k), ves(1,1,k) )
c------------------
      diag = zero
      do  i=1,norb
        diag = diag + dmat(1,i,i,k)*ves(i,i,k)
      enddo
c
      offd = zero
      do  i=2,norb
        do  j=1,i-1
          offd = offd
     $         + ( dmat(1,i,j,k) + dmat(1,j,i,k) ) * ves(i,j,k)
     $         + (-dmat(2,i,j,k) + dmat(2,j,i,k) ) * ves(j,i,k)
        enddo
      enddo
c
      ves1ctr = ves1ctr + diag + offd
c------------------
        call DO0TRACE( ves1c0,  norb, natm,ntyp,nshld,
     $   itypa,numshl,lshel, occij,wtk(k), ves(1,1,k) )
c
      enddo
c
c    That's all Folks!
c
      RETURN
      END
