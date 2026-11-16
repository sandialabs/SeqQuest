c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMORDR
c
c
      subroutine SYMORDR( natm,nsyma,nsymg, isymgrd,rmatsym,naofsym )
c---------------------------------------------------------------
c Purpose: reorder symmetries to put symmetries in reduced grid
c          symmetry group at start of full atomic symmetry group
c
c Written: Peter A. Schultz, 6-January-1999, for v2.30(redo sym)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Variables:
c   natm  - # of atoms
c   nsyma - # of (non-id) elements in atomic symmetry group
c   nsymg - # of (non-id) elements in grid symmetry group
c   isymgrd() - list of elements in grid symmetry group
c   rmatsym() - cartesian representations of symmetry operators
c   naofsym() - mapping of atoms to other atoms by symmetry
c
c  Note: that naofsym begins with the identity, while the
c  cartesian symmetry operations skip the identity, so that
c  the indexing for naofsym and rmatsym differ by 1.
c  Also note that the location of the grid symmetry elements
c  in the total list is strictly increasing, and the scheme
c  used here depends on this fact.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c input/output arrays:
      DIMENSION  isymgrd(*)
      DIMENSION  rmatsym(12,*), naofsym(natm,*)
c
c >>>> EXECUTABLE CODE
c
      if( nsymg.lt.1 ) RETURN
c
      do 100 isymg=1,nsymg
        isyma = isymgrd(isymg)
        if( isyma.ne.isymg )then
          do  iatm=1,natm
            itmp = naofsym(iatm,isymg+1)
            naofsym(iatm,isymg+1) = naofsym(iatm,isyma+1)
            naofsym(iatm,isyma+1) = itmp
          enddo
          do  ij=1,12
            tmp = rmatsym(ij,isymg)
            rmatsym(ij,isymg) = rmatsym(ij,isyma)
            rmatsym(ij,isyma) = tmp
          enddo
        endif
  100 continue
c
c    That's all Folks!
c
      RETURN
      END
