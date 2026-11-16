c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRBOX
c
c
      subroutine WRBOX( iwffil, ncplx,nk, nrbox,nobox, wf )
c---------------------------------------------------------------
c Purpose: put/get grid orbitals to/from disk
c
c Written: Peter A. Schultz, 28-May-1997, for v2.20 (boxes)
c
c Revision history:
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wf
c
      DIMENSION  wf(*)
c
      if( nobox.eq.0 ) RETURN
c
      nrorb = ncplx*nrbox
      nr0 = 0
      do  iorb=1,nobox
        do  k=1,nk
          write(iwffil)  (wf(ir),ir=nr0+1,nr0+nrorb)
          nr0 = nr0 + nrorb
        enddo
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDBOX
c
c
      entry RDBOX( iwffil, ncplx,nk,nrbox,nobox, wf )
c
      if( nobox.eq.0 ) RETURN
c
      nrorb = ncplx*nrbox
      nr0 = 0
      do  iorb=1,nobox
        do  k=1,nk
          read (iwffil)  (wf(ir),ir=nr0+1,nr0+nrorb)
          nr0 = nr0 + nrorb
        enddo
      enddo
c
      RETURN
      END
