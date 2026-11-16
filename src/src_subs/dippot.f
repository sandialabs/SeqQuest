c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DIPPOT
c
c
      subroutine DIPPOT( ion_opt, dipole,rdipole, ndim, esion,
     $ rprim, n1r,n2r,n3r,  hh, weight, ndbox,nbox,
     $ rhion, vion, r )
c---------------------------------------------------------------
c Purpose: Compute potential on grid from dipole correction
c          treating background counterdipole as gaussians.
c
c Written: Peter A. Schultz,  2-Movember-1998 (for v2.27)
c
c Revision history:
c   8Mar02-PAS/2.52: clean unused variables/literal constants d0
c   7Feb99-PAS/2.42: bugfix for bulk dipole range
c  29Mar99-PAS/2.32: to accomodate periodic systems
c  15Dec98-PAS/2.30: put into boxes
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Nov98-PAS/2.27: position of dipole included
c   5Nov98-PAS/2.27: 2D slab dipole correction built
c---------------------------------------------------------------
c
c  Notes:
c   What I will do here to cancel the grid charge dipole is
c   construct a charge density with the same dipole from a pair
c   of Gaussians, of opposite sign, centered about the origin,
c   and separated by the distance 2*rdip, where rdip is specified
c   in this routine.  I will take this dipole density and remove
c   it from the delta-density, leaving the fft density dipole-free.
c   This routine should only be necessary for neutral systems,
c   as the dipole is corrected in charged systems by shifting the
c   origin of the neutralizing charge.  2Nov98-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  dipole(3), rdipole(3)
c  grid generating primitive vector, and unit cell vectors
      DIMENSION  rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3)
c output arrays:
c  neutralizing density, and associated electrostatic potential:
      DIMENSION  rhion(*),vion(*)
c scratch arrays:
      DIMENSION  r(3,*)
c local declarations:
      DIMENSION  zchrg(2),rchrg(3,2)
      DIMENSION  orig(3), gg(3,3)
      DIMENSION  nrx(3)
c
      DATA  zero,half,one / 0.d0,0.5d0,1.d0 /
      DATA  four / 4.d0 /
c Exponent for Gaussian neutralizing charge:
      DATA  alfion / 0.30d0 /, rdip / 1.5d0 /
c   Should be small enough to be well integrable on any reasonable
c   grid but large enough to be contained within any viable unit cell.
c
c >>>> EXECUTABLE CODE:
c
      esion = zero
      dipmag = SQRT( dipole(1)**2 + dipole(2)**2 + dipole(3)**2 )
      if( ion_opt.le.0 .or. (dipmag .lt. 2.d-2) ) RETURN
c
c  We have a finite cluster system for which we wish to construct a
c  density that will cancel the dipole.
c
c Set up some constants:
      pi = four*ATAN( one )
      fourpi = pi + pi + pi + pi
      rtalf = SQRT( alfion )
      aldpi = alfion / pi
c
c Get me the net charge I will need to place in each Gaussian:
      zgaus = half*dipmag/rdip
c
c  Get range of Gaussians on the grid:
      rngalf = SQRT( cutgrid/alfion )
c
      call CROSS( hh(1,2),hh(1,3), gg(1,1) )
      call CROSS( hh(1,3),hh(1,1), gg(1,2) )
      call CROSS( hh(1,1),hh(1,2), gg(1,3) )
      do  j=1,3
        gmag = SQRT( gg(1,j)**2 + gg(2,j)**2 + gg(3,j)**2 )
        gg(1,j) = gg(1,j)/gmag
        gg(2,j) = gg(2,j)/gmag
        gg(3,j) = gg(3,j)/gmag
      enddo
c
      do  id=1,3
        rhh = ABS( gg(1,id)*hh(1,id) + gg(2,id)*hh(2,id)
     $           + gg(3,id)*hh(3,id) )
        rnx = rngalf/rhh + half
        nrx(id) = rnx
      enddo
c
c  Branch to correct location:
c
      if( ion_opt.eq.2 )then
c
c First, put center of neutralizing dipole at center of dipole:
        do  j=1,3
          orig(j) = rdipole(j)
        enddo
c
c Give charge and position of our two-Gaussian array:
        nglm = 0
        zchrg(nglm+1) = -zgaus
        zchrg(nglm+2) =  zgaus
        do  j=1,3
          rchrg(j,nglm+1) = orig(j) - rdip* dipole(j)/dipmag
          rchrg(j,nglm+2) = orig(j) + rdip* dipole(j)/dipmag
        enddo
        nglm = nglm + 2
c
CC The following check is wrong for the new WS scheme in bulk
Cc  Check position of coordinate origin on the grid:
C        do 200 iglm=nglm-1,nglm
C          do  j=1,3
C            r1(j) = rchrg(j,iglm)
C          enddo
C          call PROJLV( r1, hh(1,1), hh(1,2), hh(1,3), t )
C          do  id=1,3
C            if( t(id).lt.zero )
C     $       call STOPXERR( 'DIPPOT: LMCC gaussian outside unit cell' )
C          enddo
C  200   continue
c
c  Go ahead and compute Gaussian LMCC array:
        call IONGPOT( ndim, esion, nglm,zchrg,rchrg,
     $   rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $   rhion, vion, r )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
