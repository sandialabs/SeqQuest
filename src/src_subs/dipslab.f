c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DIPSLAB
c
c
      subroutine DIPSLAB( ion_opt,ndim,
     $ dipole,rdipole,  esion,vacplus,
     $ rprim, n1r,n2r,n3r,  hh, weight, ndbox,nbox,
     $ rhlmcc, vlmcc, r )
c---------------------------------------------------------------
c Purpose: Compute countercharge density and potential to treat
c          slab dipole moment
c
c Written: Peter A. Schultz, 30-March-1999 (for v2.32)
c
c Revision history:
c   1Jul11-PAS/2.62: adapted for consistency w/ double-layer
c  15Feb11-PAS/2.62: naming and cosmetics, slab dipole vacuum potl
c   8Mar02-PAS/2.52: clean unused variables/literal constants d0
c---------------------------------------------------------------
c
c  Notes:
c   To compute potential due to a slab dipole, construct a pair
c   of one-dimensional Gaussians planes, of opposite sign, separated
c   by the distance 2*dzlmcc. dzlmcc is specified in this routine.
c   Dipole potential is computed here, remove the lmcc dipole density
c   from the delta-density, leaving the fft density dipole-free.
c
c   Implements method set out in: P.A.Schultz, PRB 60, 1551 (1999).
c
c   Key assumption: z-axis is surface normal
c
c   Should be merged with SLABCHRG after code is vetted.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input arrays:
      DIMENSION  dipole(3), rdipole(3)
c  unit cell vectors, and grid generating basis vectors
      DIMENSION  rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3)
c Output arrays:
c  lmcc slab dipole density, and associated electrostatic potential:
      DIMENSION  rhlmcc(*),vlmcc(*)
c
c Scratch arrays:
      DIMENSION  r(3,*)
c
c Local declarations:
      LOGICAL    ingaus
      DIMENSION  r0box(3),boxctr(3)
      DIMENSION  origlm(3), rvec3(3), gg(3,3)
      DIMENSION  nro(3),nrx(3)
c
      DATA  four / 4.d0 /
c
c   Parameters for slab lmcc:
      DATA  alflmcc / 0.30d0 /, dzlmcc / 1.5d0 /
c      alflmcc = Gaussian exponent for dipole charges
c      dzlmcc  = z-offset from dipole center for lmcc gaussian planes
c   These are adjustable parameters, constrained by the need to be ...
c    ... small enough to be well integrable on any reasonable grid
c    ... but large enough to be contained within any viable unit cell
c
c >>>> EXECUTABLE CODE:
c
      esion = zero
      vacplus = zero
c
c Routine is for neutral slab, with planar dipole
c
      if( ndim.ne.2 .or. ion_opt.ne.2 ) RETURN
c
c Set up some constants 
      fourpi = four*pi
      rtalf = SQRT( alflmcc )
      alfdpi = alflmcc / pi
c
c Compute slab normal vector
      call CROSS( rprim(1,1),rprim(1,2), rvec3 )
      area = SQRT( rvec3(1)**2 + rvec3(2)**2 + rvec3(3)**2 )
      vol = rvec3(1)*rprim(1,3) + rvec3(2)*rprim(2,3)
     $    + rvec3(3)*rprim(3,3)
      if( vol .lt. zero )then
        rvec3(1) = -rvec3(1)
        rvec3(2) = -rvec3(2)
        rvec3(3) = -rvec3(3)
      endif
      rvec3(1) = rvec3(1)/area
      rvec3(2) = rvec3(2)/area
      rvec3(3) = rvec3(3)/area
c
c Project dipoles and grid spacing onto slab normal (z-vector)
      z0dipole = rvec3(1)*rdipole(1) + rvec3(2)*rdipole(2)
     $         + rvec3(3)*rdipole(3)
      zhh = rvec3(1)*hh(1,3) + rvec3(2)*hh(2,3) + rvec3(3)*hh(3,3)
c
      dipolez = rvec3(1)*dipole(1) + rvec3(2)*dipole(2)
     $        + rvec3(3)*dipole(3)
      dipzmag = ABS( dipolez )
c
      if( dipzmag .lt. 2.d-2 ) RETURN
c
c----------------------------------------------------------------------
c  Slab system has significant dipole.
c  Construct lmcc density that matches dipole normal to the surface,
c  and compute associate potential from this density.
c----------------------------------------------------------------------
c
c Compute net charge needed for each Gaussian to cancel dipole
      qgaus = dipzmag * half/dzlmcc
c
c Get range of lmcc Gaussians on the grid:
      rngalf = SQRT( cutgrid/alflmcc )
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
      if( ion_opt .eq. 2 )then
c
c      ***** Slab system *****
c
       pole = one
 2100   continue
c       Start with negative-z pole of dipole, then switch to z(+)
        pole = -pole
        qgaus = -qgaus
c
c       First, offset gaussian pole from center of neutralizing dipole:
        zgaus = - z0dipole - pole*(dipolez/dipzmag)*dzlmcc
c
        nro(3) = (-zgaus/zhh) + half + one
c
        ic1 = nro(3) - nrx(3)
        icn = nro(3) + nrx(3)
        if( ic1.lt.1 ) ic1 = 1
        if( icn.gt.n3r ) icn = n3r
c
        gcoef = qgaus* SQRT( alfdpi ) / area
c       Coefficients for gaussian and error function terms
        vgcoef = -gcoef*fourpi / alflmcc
        vecoef = -gcoef*fourpi / SQRT( alfdpi )
c
        origlm(1) = zero
        origlm(2) = zero
        origlm(3) = zgaus
c
        esiongr = zero
        irpot = 0
        do 2250 ibox=1,nbox
c
c         Get basic box info
          call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $     ira0,irb0,irc0, idra,idrb,idrc )
          nrbox = idra*idrb*idrc
c
          iv0 = irpot
          irpot = irpot + nrbox
c
          ingaus = .true.
          if( (irc0+idrc).lt.ic1 .or. icn.le.irc0 ) ingaus = .false.
c
c         Generate coordinates of points within box
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          r0box(1) = r0box(1) + origlm(1)
          r0box(2) = r0box(2) + origlm(2)
          r0box(3) = r0box(3) + origlm(3)
          call MKRBOX( hh,r0box, idra,idrb,idrc, r )
c
          if( ingaus )then
c           This box is within the gaussian range
c
            rhcoef = gcoef
            do  nr=1,nrbox
              zval = r(3,nr)
              zmag =  ABS( zval )
              zsq = zmag*zmag
              rhadd = rhcoef* EXP( -alflmcc*zsq )
              vadd = vgcoef* EXP( -alflmcc*zsq )
     $         + vecoef*zmag*DERF( rtalf*zmag )
              rhlmcc(iv0+nr) = rhlmcc(iv0+nr) + rhadd
              vlmcc(iv0+nr) = vlmcc(iv0+nr) + vadd
              esiongr = esiongr + vlmcc(iv0+nr)*rhlmcc(iv0+nr)
            enddo
          else
c           Outside gaussian range, asymptotic potential
c
            do  nr=1,nrbox
              zval = r(3,nr)
              zmag =  ABS( zval )
              rhadd = zero
              vadd = vecoef*zmag
              vlmcc(iv0+nr) = vlmcc(iv0+nr) + vadd
              esiongr = esiongr + vlmcc(iv0+nr)*rhlmcc(iv0+nr)
            enddo
          endif
c
c         Close loop over boxes
 2250   continue
c
c       Compute self-energy of dipole on grid:
        esiongr = weight*esiongr/two
c
        esion = esiongr
c
       if( pole.lt.zero ) goto 2100
c
c      Compute asymptotic vacuum potential (due to dipole)
       vacplus = -vecoef*two*dzlmcc*(dipolez/dipzmag)
C       write(IWR,*) 'Vacuum(-z),vacuum(+z)[Ry]=', -vacplus, vacplus
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
