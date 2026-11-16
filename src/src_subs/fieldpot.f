c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FIELDPOT
c
c
      subroutine FIELDPOT( efield, r0,nda,ndb,ndc,hh, vfield )
c---------------------------------------------------------------
c Purpose: generate potential on regular grid from electric field
c
c Written: Peter A. Schultz and Jamil Tahir-Kheli
c          17-September-2003 for v2.57
c
c Revision history:
c  17Sep03-PAS/2.57: pushed into subroutine; 3-vector E-field
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
c   efield - electric field, (Ry/bohr)
      DIMENSION  efield(3)
c   r0 - position of grid origin w.r.t. field origin.
      DIMENSION  r0(3)
c   hh(3,3) = grid interval vectors
      DIMENSION  hh(3,3)
c
c Output: potential (Ry) of E-field on grid
      DIMENSION  vfield(*)
c
c Local:
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      dpotc = efield(1)*hh(1,3) + efield(2)*hh(2,3) + efield(3)*hh(3,3)
      dpotb = efield(1)*hh(1,2) + efield(2)*hh(2,2) + efield(3)*hh(3,2)
      dpota = efield(1)*hh(1,1) + efield(2)*hh(2,1) + efield(3)*hh(3,1)
c Do not know how I will want to use this offset yet, so blank it out:
c      potc0 = -( efield(1)*r0(1) + efield(2)*r0(2) + efield(3)*r0(3) )
      potc0 = zero
c
      nr = 0
      do 300 idc=1,ndc
        potb0 = potc0
c
        do 200 idb=1,ndb
          potr = potb0
c
          do 100 ida=1,nda
            nr = nr + 1
            vfield(nr) = vfield(nr) + potr
            potr = potr + dpota
  100     continue
c
          potb0 = potb0 + dpotb
  200   continue
c
        potc0 = potc0 + dpotc
  300 continue
c
      RETURN
c
c    That's all Folks!
c
      END
