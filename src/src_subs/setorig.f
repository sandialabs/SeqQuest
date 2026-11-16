c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SETORIG
c
c
      subroutine SETORIG( dorig,orig,rsym, origws,rchrg,
     $ ndim, dolvatm,rprim,rscale )
c---------------------------------------------------------------
c Purpose: initialize the various origins we have
c
c Written: Peter A. Schultz, 23-May-2001, for v2.47
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      LOGICAL    dolvatm
      DIMENSION  rprim(3,3),rscale(3)
      DIMENSION  dorig(3),orig(3),rsym(3)
      DIMENSION  origws(3),rchrg(3)
c
c local declarations
      DIMENSION  dorig0(3),rsym0(3)
      DIMENSION  vec0(3)
      DATA       vec0 / 0.d0,0.d0,0.d0 /
      DATA       half / 0.5d0 /
c
c >>>> EXECUTABLE CODE:
c
c  Transform supercell origin offset into internal coordinates
c  NB: do NOT shift by orig(), i.e., itself: pass zero offset vector
      call RINTRNL( ndim,  1 , dolvatm, rprim,vec0,rscale,
     $ dorig,dorig0 )
      do  i=1,3
        dorig(i) = dorig0(i)
        orig(i) = dorig(i)
      enddo
c  Offset origin by half lattice vectors in non-periodic directions:
      do  id=ndim+1,3
        orig(1) = orig(1) - half*rprim(1,id)
        orig(2) = orig(2) - half*rprim(2,id)
        orig(3) = orig(3) - half*rprim(3,id)
      enddo
c  Transform center of symmetry into internal coordinates
      call RINTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $ rsym,rsym0 )
      do  i=1,3
        rsym(i) = rsym0(i)
      enddo
c
c  Set wigner-seitz origin to initial charge location
c    (or center of unit cell if not charged)
c    Note that we will convert to internal units another time ...
      origws(1) = rchrg(1)
      origws(2) = rchrg(2)
      origws(3) = rchrg(3)
c
c    That's all Folks!
c
      RETURN
      END
