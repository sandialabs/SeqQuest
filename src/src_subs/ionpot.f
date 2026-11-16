c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IONPOT
c
c
      subroutine IONPOT( ion_opt, elchrg,dipole, ndim, esion, rchrg,
     $ rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $ rhion, vion, r )
c---------------------------------------------------------------
c Purpose: Compute potential on grid from background charge
c          treating background countercharge as a gaussian, or
c          as a flat jellium countercharge (treating the potential
c          as either flat, or spherical)
c
c Written: Peter A. Schultz, 27-October-1998 (for v2.27)
c
c Revision history:
c  18May11-PAS/2.62: some cosmetics
c  14Jun05-PAS/2.59: File units
c   8Mar02-PAS/2.52: clean unused variables/literal constants d0
c  21Jun01-PAS/2.48: replace STOPs
c   8Feb00-PAS/2.42: bulk charge position bugfix
c  26Jul99-PAS/2.37: fix diagnostic/send out charge position right
c  29Mar99-PAS/2.32: adapt for dipole in non-cluster systems
c  11Jan99-PAS/2.32: charge states for non-cluster systems
c  15Dec98-PAS/2.30: put into boxes
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
c Input arrays:
      DIMENSION  dipole(3), rchrg(3)
c   grid generating primitive vector, and unit cell vectors
      DIMENSION  rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3)
c Output arrays:
c   neutralizing density, and associated electrostatic potential:
      DIMENSION  rhion(*),vion(*)
c
c Scratch arrays:
      DIMENSION  r(3,*)
c
c Local declarations:
      DIMENSION  r0box(3),boxctr(3)
      DIMENSION  origlm(3), t(3)
c
      DATA  three,four / 3.d0,4.d0 /
c  Exponent for Gaussian neutralizing charge:
c   Should be small enough to be well integrable on any reasonable
c   grid but large enough to be contained within any viable unit cell.
      DATA  alfion / 0.30d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
      esion = zero
      nptr = n1r*n2r*n3r
c
c  We have a charged system, with a local point charge,
c  we wish to neutralized for input into the fft Poisson solver
c
      if( elchrg .eq. zero ) RETURN
c
      qgaus = -elchrg
c
c Set up some constants:
      pi = four*ATAN( one )
      fourpi = four*pi
      rtalf = SQRT( alfion )
      aldpi = alfion/pi
c
c Offset origin by enough to cancel the dipole, if requested:
c
c  First, set center of neutralizing charge at center of cell:
      do  j=1,3
        origlm(j) = rchrg(j)
      enddo
c
      dipsq = zero
      do  idim=1,3
        dipsq = dipsq + dipole(idim)**2
      enddo
c     Use a threshhold to skip dipole due to numerical noise
      if( ion_opt.gt.0 .and. dipsq.gt. 4.d-4 )then
c       Note you have to be careful of qgaus sign
        do  idim=1,3
          origlm(idim) = origlm(idim) + dipole(idim) / qgaus
        enddo
      endif
c
c  Branch to chosen method for Poisson solver for charge
c
      if( ion_opt.eq.0 )then
c
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c             Jellium neutralizing charge, flat potential
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       I do not know how to explicitly evaluate the jellium self-energy
c       integral (need Madelung term, Ewald summations, etc), but I have
c       an empirically derived expression that converges to the correct
c       result.  Since this branch is superceded by the more rigorous
c       gaussian local moment countercharge scheme, no better is needed.
c       Someone who cares more about this can refine this in the future.
c
        call CROSS( rprim(1,2),rprim(1,3), t )
        vol = ABS( rprim(1,1)*t(1) + rprim(2,1)*t(2) + rprim(3,1)*t(3) )
        arg = 9*pi / (two*vol)
c
        essph = 0.8d0*qgaus*qgaus*arg**(one/three)
        esion = 0.925368d0 * essph
        sigma = qgaus/vol
c       Make potential so that grid integral matches self-energy term:
        vsigma = two*esion / qgaus
c
        esions = zero
        do  nr=1,nptr
          rhion(nr) = sigma
          vion(nr) = vsigma
          esions = esions + rhion(nr)*vion(nr)
        enddo
        esiongr = nptr*half*weight*sigma*vsigma
c
        esions = half*weight*esions
        rhoion = nptr*weight*sigma
        if( lstout .gt. 0 )then
          call FLGETIWR( IWR )
          write(IWR,'(a,4i10)')  'IONPOT: n1r,n2r,n3r=',n1r,n2r,n3r
          write(IWR,'(a,f20.8)') 'IONPOT: flat charge =',rhoion
          write(IWR,'(a,f20.8)') 'IONPOT: flat esion  =',esion
          write(IWR,'(a,f20.8)') 'IONPOT: flat esions =',esions
          write(IWR,'(a,f20.8)') 'IONPOT: flat esiongr=',esiongr
        endif
c
      elseif( IABS( ion_opt ) .eq. 1 )then
c
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Flat background charge, sphere-shaped potential
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Our unit cell is a parallelpiped.  To simplify the arithmetic,
c       assume the potential follows the functional form of a sphere.
c       Further, for computing the grid potential, I'm going to assume
c       that functional form for the potential is the functional form
c       for the inside of the sphere, even for points lying outside the
c       sphere radius (rather than the Z/r proper outside a sphere).
c       This is a cheaper calculation (no branching) and I can get away
c       with it because in a local Gaussian basis, the orbitals have no
c       amplitude there, which, further, explains why I do not care
c       about discontinuities in the potential at the boundaries.
c
        call CROSS( rprim(1,2),rprim(1,3), t )
        vol = ABS( rprim(1,1)*t(1) + rprim(2,1)*t(2) + rprim(3,1)*t(3) )
        arg = 9*pi / (two*vol)
c
c  For the self-energy correction, I use the functional form for a
c  sphere, and scale with an empirically derived factor to get the
c  proper energy correction for a cube unit cell.
        essph = 0.8d0*qgaus*qgaus*arg**(one/three)
        esion = 0.925368 * essph
c  Get sphere radius matching unit cell volume, density of background:
        rad = ( (three*vol) / (fourpi) ) ** (one/three)
        radsq = rad*rad
        sigma = qgaus/vol
c
c  Compute potential (Ry) based on distance from sphere origin:
        fac = fourpi*sigma
c
        irpot = 0
        do 150 ibox=1,nbox
c
c  Get basic box info
          call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $     ira0,irb0,irc0, idra,idrb,idrc )
          nrbox = idra*idrb*idrc
c
          iv0 = irpot
          irpot = irpot + nrbox
c
c  Generate coordinates
c
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          r0box(1) = r0box(1) - origlm(1)
          r0box(2) = r0box(2) - origlm(2)
          r0box(3) = r0box(3) - origlm(3)
          call MKRBOX( hh,r0box, idra,idrb,idrc, r )
c
          do 140 nr=1,nrbox
            rsq = r(1,nr)**2 + r(2,nr)**2 + r(3,nr)**2
c
c  Flat neutralizing charge, and associated potential (sphere):
            rhion(iv0+nr) = sigma
            vion(iv0+nr) = fac*( radsq - rsq/three )
  140     continue
c
  150   continue
c
      elseif( IABS( ion_opt ) .eq. 2 )then
c
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               LMCC - 0D point charge - Gaussians
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
CC Following  check doesn't work anymore with new conventions.
Cc  Find position of Gaussian LMCC origin on the grid:
C        do id=1,3
C          r1(id) = origlm(id)
C        enddo
Cc
C        call PROJLV( r1, hh(1,1), hh(1,2), hh(1,3), t )
C        do 160 id=1,3
C          if( t(id).lt.zero .and. id.gt.ndim )then
C            call STOPXERR( 'LMCC-ion/Gaussian ctr outside cell' )
C          endif
C  160   continue
c
c  Go ahead and compute single Gaussian LMCC:
c
        call IONGPOT( ndim, esion, 1,qgaus,origlm,
     $   rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $   rhion, vion, r )
c
      else
c
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Not a valid charge treatment option
c       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        call STOPXERR( 'ionopt-x/bad ion treatment option/IONPOT' )
c
      endif
c
c  Send out new position of charge:
      do  j=1,3
        rchrg(j) = origlm(j)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
