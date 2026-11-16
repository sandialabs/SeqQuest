c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLABCHRG
c
c
      subroutine SLABCHRG( ion_opt,ndim,
     $ elchrg, origws, zchrglm,zchrgcc,
     $ dipole,rdipole,  esion,vacplus,
     $ rprim, n1r,n2r,n3r,  hh, weight, ndbox,nbox,
     $ rhlmcc, vlmcc, r )
c---------------------------------------------------------------
c Purpose: Apply above-surface slab counter-charge plane to
c          balance charge, compute dipole potential from dipole
c          planes and remove the surface charge from the density.
c
c Written: Peter A. Schultz, 22-February-2011 (for 2.62)
c          (Initial version adapted from DIPSLAB---nuetral slab)
c
c Revision history:
c  14Sep11-PAS/2.62: finally debugged and consistent with DIPSLA
c---------------------------------------------------------------
c
c  Notes:
c   Construct a pair of one-dimensional Gaussian planes:
c   one at the centroid (zchrglm) of the net charge on the slab
c   and a counter-charge plane at specified point (zchrgcc) in vacuum.
c   Compute the full LMCC dipole potential here.
c   Compute the net charge plane to neutralize the slab net charge.
c
c   Extends method set out in: P.A.Schultz, PRB 60, 1551 (1999),
c   which describes the construction of a dipole potential
c
c   Key assumption: z-axis is surface normal
c
c   Makes no checks if choices of counterslab positions are sensible
c
c   Should be merged with DIPSLAB after code is vetted.
c
c  NB: ASSUMPTION INTO CODE THAT Z_VECTOER IS SLAB NORMAL
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
      DIMENSION  origws(3)
c
      DIMENSION  zchrglm(3), zchrgcc(3)
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
cdef  Not ready for prime time yet, but getting ready ...
      PARAMETER  (nlmd=2)
      DIMENSION  ilmq(nlmd), qlmcc(nlmd), zlmcc(nlmd)
c
      DATA  four / 4.d0 /
c
c   Parameters for slab lmcc:
      DATA  alflmcc / 0.30d0 /, dzlmcc / 1.5d0 /
      DATA  zlmaxdif / 10.0d0 /
c      alflmcc = Gaussian exponent for dipole charges
c      dzlmcc  = z-offset from dipole center for lmcc gaussian planes
c      zlmaxdif = largest jump in slab charge plane
c   These are adjustable parameters, constrained by the need to be ...
c    ... small enough to be well integrable on any reasonable grid
c    ... but large enough to be contained within any viable unit cell
c
c >>>> EXECUTABLE CODE:
c
      esion = zero
      vacplus = zero
c
c Routine is for charged slab, with planar countercharge
c
cxxx:start
      call FLGETIWR( IWR )
      write(IWR,*) 'DBG/SLABCHRG: ====== elchrg,ion_opt=',elchrg,ion_opt
  913 format(1x,'  SLABQ: ',a,6f12.6)
cxxx:end
      if( elchrg .eq. zero ) RETURN
      if( ndim.ne.2 .or. IABS( ion_opt ) .ne. 3 ) RETURN
c
      qgaus = -elchrg
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
cxxx:start
      write(IWR,913) 'dipole=',dipole
      write(IWR,913) 'dipolez(revc3.dipole)=',dipolez
cxxx:end
      dipzmag = ABS( dipolez )
c
c Determine z-location of charge in slab, for local moment charge
      ilmq(1) = 1
      zgauslm = origws(3) + dipolez / qgaus
cxxx:start
      write(IWR,913) 'old zlm zchrglm=',(zchrglm(i),i=1,3)
      write(IWR,913) 'new zlm zgauslm=',zgauslm
cxxx:end
      zlmdif = zgauslm - zchrglm(3)
      if( ABS( zlmdif ) .gt. zlmaxdif )then
c       Regulate this to make it more stable
        zgauslm = zchrglm(3) + zlmaxdif*( zlmdif/ABS(zlmdif) )
        write(IWR,913) 'regulated zlm zgauslm=',zgauslm
      endif
      zchrglm(1) = zero
      zchrglm(2) = zero
      zchrglm(3) = zgauslm
c Load z-location of countercharge plane (xxx worry about offset)
      ilmq(2) = 0
      zgauscc = zchrgcc(3)
c
c----------------------------------------------------------------------
c  Slab system has a charge.
c  Construct one lmcc density plane for slab charge, with potential
c  Construct potential (no density) from second counterdensity plane in vacuum
c----------------------------------------------------------------------
c
c Hijack vetted DIPSLAB machinery to better ensure code works
      z0dipole = ( zgauslm + zgauscc ) / two
      dzlmcc = zgauslm - z0dipole
cxxx:start dev debug diagnostics
      call FLGETIWR( IWR )
      write(IWR,913) 'z_lm,z_cc=',zgauslm,zgauscc
      write(IWR,913) 'z0_lm,dz_lm=',z0dipole,dzlmcc
      write(IWR,913) '  qgaus=',qgaus
      write(IWR,913) '  dipolez,dipzmag=',dipolez,dipzmag
      rhlmsum = 0.d0
      nptslm = 0
cxxx:end
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
      if( IABS( ion_opt ) .eq. 3 )then
c
c      ***** Slab system *****
c
       ilmgaus = 0
       pole = one
c      Do ign weidity to align internals of loop with dipslab
       qgaus = -qgaus
       zgauspole = zgauslm
 2100   continue
c       Start with negative-z pole of dipole, then switch to z(+)
        ilmgaus = ilmgaus + 1
        pole = -pole
        qgaus = -qgaus
c
c       First, offset gaussian pole from center of neutralizing dipole:
C        zgaus = - z0dipole - pole*(dipolez/dipzmag)*dzlmcc
        zgaus = -zgauspole
        zgauspole = zgauscc
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
cxxx:start
      rhlmtot = 0.d0
      velmtot = 0.d0
      nboxin = 0
cxxx:end
        esiongr = zero
        esioncc0 = zero
        esioncc = zero
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
c           Turn off lm density for off-slab countercharge
cxxx:
       nboxin = nboxin + 1
            rhcoef = gcoef
            if( ilmq(ilmgaus) .eq. 0 ) rhcoef = zero
            do  nr=1,nrbox
              zval = r(3,nr)
              zmag =  ABS( zval )
              zsq = zmag*zmag
              expaz2 =  EXP( -alflmcc*zsq )
cxxx:
C              rhadd = rhcoef* EXP( -alflmcc*zsq )
              rhadd = rhcoef* expaz2 
C              vadd = vgcoef* EXP( -alflmcc*zsq )
              vadd = vgcoef* expaz2
     $         + vecoef*zmag*DERF( rtalf*zmag )
cxxx:
       rhlmtot = rhlmtot + rhadd
       velmtot = velmtot + vadd
       esioncc0 = esioncc0 + rhlmcc(iv0+nr)*vadd
       esioncc = esioncc + (gcoef*expaz2)*vlmcc(iv0+nr)
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
       esioncc0 = esioncc0 + rhlmcc(iv0+nr)*vadd
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
cxxx
        esioncc0 = weight*esioncc0 / two
        esioncc = weight*esioncc / two
c
        esion = esiongr
cxxx:
      write(IWR,913) '  q_gaus,z_gaus= ',qgaus,zgaus
      write(IWR,*) '   DBG/SLABQ: rhlmtot,vlmtot,nboxin',
     $ rhlmtot*weight,velmtot*weight, nboxin
CC      write(IWR,913) 'DBG: aft, vgeof,vecoef=',vgcoef,vecoef
c
       if( pole.lt.zero ) goto 2100
c
c      Compute asymptotic vacuum potential (due to dipole)
       vacplus = -vecoef*two*dzlmcc*(dipolez/dipzmag)
c      Add dble layer energy into (minus sign due to attractive dipole).
cxxx: this sign is problematic ....
       esion = esion + esioncc0
cxxx:
       write(IWR,913) 'dipole vacuum offset vacplus=',vacplus
       write(IWR,913) 'esiongr,esioncc,esion(Ry)=',
     $  esiongr,esioncc0,esion
C       write(IWR,*) 'Vacuum(-z),vacuum(+z)[Ry]=', -vacplus, vacplus
c
      endif
cxxx:start
C      call GRIDAVE( IWR, n1r,n2r,n3r, 1, rhlmcc, ave, r, 'RHLMCC' )
C      call GRIDAVE( IWR, n1r,n2r,n3r, 1, vlmcc, ave, r, 'VLMCC' )
      write(IWR,*) 'DBG/SLABCHRG =========== end slabchrg ==========='
cxxx:end
c
c    That's all Folks!
c
      RETURN
      END
