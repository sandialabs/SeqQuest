c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKGRHO
c
c
      subroutine MKGRHO( idrhfl0,idrhofl,   nptr,
     $ grho,g2rho,gpggp,  drho,ddrho )
c---------------------------------------------------------------
c Purpose: compute, on grid, density gradients needed for gga
c          calculation from density derivatives on grid.
c
c Written: Peter A. Schultz, 11-December-2000, for v2.46 (gga)
c
c Revision history:
c  20Jul01-PAS/2.49: spin-polarized dft. mv idrhofl rewind up
c---------------------------------------------------------------
c
c Notes:
c  We need three quantities:
c   (1) grho  = | grad[rho] |
c   (2) g2rho = grad**2 [rho]
c   (3) gpggp = grad[rho] DOT grad[ | grad[rho] | ]
c  On disk(s) we have all first and second derivatives of rho:
c    d/dx,d/dy,d/dz
c    d2/dxdx, d2/dydy, d2/dzdz, d2/dxdy, d2/dydz, d2/dzdx
c  Need to accumulate:
c    grho  = SQRT( (dp/dx)**2 + (dp/dy)**2 + (dp/dz)**2 )
c    g2rho = d2p/dxdx + d2p/dydy + d2p/dzdz
c    gpggp = SUM[ du*dv*duv ] / grho
c  For vectorization, I accumulate du*dv*duv without regard
c  for grho=0, and take care of this at the end where I do the
c  final zero-checking and division by grho.
c  I need only a single second derivative field in memory
c  at any one time (plus a scratch work field to get the
c  reference + scf second derivatives), and accumulate
c  needed sums and discard. I compute grho late, so that
c  I can use it as temp scratch for this purpose.
c
c Variables:
c    idrhfl0 - file with derivatives from reference density
c    idrhofl - file with derivatives from delta (scf) density
c    nptr    - size of grid field
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, half = 0.5d0, two = 2.d0 )
c output arrays: (gga gradients)
      DIMENSION  grho(nptr),g2rho(nptr),gpggp(nptr)
      DIMENSION  wkrd(*)
c scratch arrays: (1st derivs, 2nd deriv)
      DIMENSION  drho(nptr,3)
      DIMENSION  ddrho(nptr)
c local declarations:
c     This defines threshhold for what is considered a "zero" gradient
      DATA  vsmall / 1.d-8 /
c
c >>>> EXECUTABLE CODE:
c
c Clear out gradient terms:
      call MKZERO( nptr, grho )
      call MKZERO( nptr, g2rho )
      call MKZERO( nptr, gpggp )
c
c Retrieve first derivatives of reference atom density:
      do  iu=1,3
        call READBIG( idrhfl0, nptr, drho(1,iu) )
      enddo
c
c Retrieve first derivatives of scf delta-density ...
      do  iu=1,3
        call READBIG( idrhofl, nptr, grho )
c        ... and combine to get total density derivatives
        do  nr=1,nptr
          drho(nr,iu) = drho(nr,iu) + grho(nr)
        enddo
      enddo
c
c Retrieve double (XX,YY,ZZ) derivs of reference and scf density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhofl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr) + grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          g2rho(nr) = g2rho(nr) + ddrho(nr)
          gpggp(nr) = gpggp(nr) + drho(nr,iu)*drho(nr,iu)*ddrho(nr)
        enddo
      enddo
c
c Retrieve 2nd (XY,YZ,ZX) derivs of reference and scf density ...
      do  iu=1,3
        iv = iu + 1
        if( iv.eq.4 ) iv = 1
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhofl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr) + grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          gpggp(nr) = gpggp(nr) + two*drho(nr,iu)*drho(nr,iv)*ddrho(nr)
        enddo
      enddo
c
c Assemble | grad[rho] | term in grho(), now that scratch work done:
      do  nr=1,nptr
        grho(nr) = SQRT( drho(nr,1)**2 + drho(nr,2)**2 + drho(nr,3)**2 )
      enddo
c
c And finish the gpggp term, with conditional division by grho:
      do  nr=1,nptr
        if( grho(nr).gt.vsmall )then
          gpggp(nr) = gpggp(nr) / grho(nr)
        else
          gpggp(nr) = zero
        endif
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKGRHO0
c
c
      entry MKGRHO0( idrhfl0,   nptr,
     $ grho,g2rho,gpggp,  drho,ddrho )
c---------------------------------------------------------------
c Purpose: compute, on grid, density gradients needed for gga
c          calculation from density derivatives on grid for
c          reference density.
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
c Clear out gradient terms:
      call MKZERO( nptr, grho )
      call MKZERO( nptr, g2rho )
      call MKZERO( nptr, gpggp )
c
c Retrieve first derivatives of reference atom density:
      do  iu=1,3
        call READBIG( idrhfl0, nptr, drho(1,iu) )
      enddo
c
c Retrieve double (XX,YY,ZZ) derivs of reference atom density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, ddrho )
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          g2rho(nr) = g2rho(nr) + ddrho(nr)
          gpggp(nr) = gpggp(nr) + drho(nr,iu)*drho(nr,iu)*ddrho(nr)
        enddo
      enddo
c
c Retrieve 2nd (XY,YZ,ZX) derivs of reference atom density ...
      do  iu=1,3
        iv = iu + 1
        if( iv.eq.4 ) iv = 1
        call READBIG( idrhfl0, nptr, ddrho )
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          gpggp(nr) = gpggp(nr) + two*drho(nr,iu)*drho(nr,iv)*ddrho(nr)
        enddo
      enddo
c
c Assemble | grad[rho] | term in grho(), now that scratch work done:
      do  nr=1,nptr
        grho(nr) = SQRT( drho(nr,1)**2 + drho(nr,2)**2 + drho(nr,3)**2 )
      enddo
c
c And finish the gpggp term, with conditional division by grho:
      do  nr=1,nptr
        if( grho(nr).gt.vsmall )then
          gpggp(nr) = gpggp(nr) / grho(nr)
        else
          gpggp(nr) = zero
        endif
      enddo
c
c    That's all Folks!
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKGRHOA
c
c
      entry MKGRHOA( idrhfl0,idrhofl,idrhsfl,   nptr,
     $ grho,g2rho,gpggp,  drho,ddrho, wkrd )
c---------------------------------------------------------------
c Purpose: compute, on grid, up-spin density gradients for gga
c          calculation from density derivatives on grid.
c
c Written: Peter A. Schultz, 11-August-2001, for v2.49 (gga)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Notes:
c  We need three quantities:
c   (1) grho  = | grad[rho] |
c   (2) g2rho = grad**2 [rho]
c   (3) gpggp = grad[rho] DOT grad[ | grad[rho] | ]
c  On disk(s) we have all first and second derivatives of rho:
c    d/dx,d/dy,d/dz
c    d2/dxdx, d2/dydy, d2/dzdz, d2/dxdy, d2/dydz, d2/dzdx
c  Need to accumulate:
c    grho  = SQRT( (dp/dx)**2 + (dp/dy)**2 + (dp/dz)**2 )
c    g2rho = d2p/dxdx + d2p/dydy + d2p/dzdz
c    gpggp = SUM[ du*dv*duv ] / grho
c  For vectorization, I accumulate du*dv*duv without regard
c  for grho=0, and take care of this at the end where I do the
c  final zero-checking and division by grho.
c  I need only a single second derivative field in memory
c  at any one time (plus a scratch work field to get the
c  reference + scf second derivatives), and accumulate
c  needed sums and discard. I compute grho late, so that
c  I can use it as temp scratch for this purpose.
c
c Variables:
c    idrhfl0 - file with derivatives from reference density
c    idrhofl - file with derivatives from total delta (scf) density
c    idrhsfl - file with derivatives from up delta (scf) density
c    nptr    - size of grid field
c
c >>>> EXECUTABLE CODE:
c
c Clear out gradient terms:
      call MKZERO( nptr, grho )
      call MKZERO( nptr, g2rho )
      call MKZERO( nptr, gpggp )
c
c Retrieve first derivatives of reference density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, drho(1,iu) )
c        ... and divide by two to get up-spin part:
        call DSCAL( nptr, half, drho(1,iu),1 )
      enddo
c
c Retrieve first derivatives of total scf delta-density ...
      do  iu=1,3
        call READBIG( idrhsfl, nptr, grho )
        call READBIG( idrhofl, nptr, wkrd )
c        ... accumulate dn-spin gradient:
        do  nr=1,nptr
          g2rho(nr) = g2rho(nr) +
     $     ( drho(nr,iu) + wkrd(nr) - grho(nr) )**2
        enddo
c        ... and combine to get total up-spin first derivatives
        do  nr=1,nptr
          drho(nr,iu) = drho(nr,iu) + grho(nr)
        enddo
      enddo
c
c Complete *dn*-spin gradient, and reinitialize borrowed array:
      do  nr=1,nptr
        wkrd(nr) = SQRT( g2rho(nr) )
      enddo
      call MKZERO( nptr, g2rho )
c
c Retrieve double (XX,YY,ZZ) derivs of reference and scf density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhsfl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr)/two + grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          g2rho(nr) = g2rho(nr) + ddrho(nr)
          gpggp(nr) = gpggp(nr) + drho(nr,iu)*drho(nr,iu)*ddrho(nr)
        enddo
      enddo
c
c Retrieve 2nd (XY,YZ,ZX) derivs of reference and scf density ...
      do  iu=1,3
        iv = iu + 1
        if( iv.eq.4 ) iv = 1
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhsfl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr)/two + grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          gpggp(nr) = gpggp(nr) + two*drho(nr,iu)*drho(nr,iv)*ddrho(nr)
        enddo
      enddo
c
c Assemble | grad[rho] | term in grho(), now that scratch work done:
      do  nr=1,nptr
        grho(nr) = SQRT( drho(nr,1)**2 + drho(nr,2)**2 + drho(nr,3)**2 )
      enddo
c
c And finish the gpggp term, with conditional division by grho:
      do  nr=1,nptr
        if( grho(nr).gt.vsmall )then
          gpggp(nr) = gpggp(nr) / grho(nr)
        else
          gpggp(nr) = zero
        endif
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKGRHOB
c
c
      entry MKGRHOB( idrhfl0,idrhofl,idrhsfl,  nptr,
     $ grho,g2rho,gpggp,  drho,ddrho, wkrd )
c---------------------------------------------------------------
c Purpose: compute, on grid, dn-spin density gradients for gga
c          calculation from density derivatives on grid.
c
c Written: Peter A. Schultz, 11-August-2001, for v2.49 (gga)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Notes:
c  We need three quantities:
c   (1) grho  = | grad[rho] |
c   (2) g2rho = grad**2 [rho]
c   (3) gpggp = grad[rho] DOT grad[ | grad[rho] | ]
c  On disk(s) we have all first and second derivatives of rho:
c    d/dx,d/dy,d/dz
c    d2/dxdx, d2/dydy, d2/dzdz, d2/dxdy, d2/dydz, d2/dzdx
c  Need to accumulate:
c    grho  = SQRT( (dp/dx)**2 + (dp/dy)**2 + (dp/dz)**2 )
c    g2rho = d2p/dxdx + d2p/dydy + d2p/dzdz
c    gpggp = SUM[ du*dv*duv ] / grho
c  For vectorization, I accumulate du*dv*duv without regard
c  for grho=0, and take care of this at the end where I do the
c  final zero-checking and division by grho.
c  I need only a single second derivative field in memory
c  at any one time (plus a scratch work field to get the
c  reference + scf second derivatives), and accumulate
c  needed sums and discard. I compute grho late, so that
c  I can use it as temp scratch for this purpose.
c
c Variables:
c    idrhfl0 - file with derivatives from reference density
c    idrhofl - file with derivatives from total delta (scf) density
c    idrhsfl - file with derivatives from up delta (scf) density
c    nptr    - size of grid field
c
c >>>> EXECUTABLE CODE:
c
c Clear out gradient terms:
      call MKZERO( nptr, grho )
      call MKZERO( nptr, g2rho )
      call MKZERO( nptr, gpggp )
c
c Retrieve first derivatives of reference density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, drho(1,iu) )
c        ... and divide by two to get dn-spin part:
        call DSCAL( nptr, half, drho(1,iu),1 )
      enddo
c
c Retrieve first derivatives of ...
      do  iu=1,3
c        ... total delta-density:
        call READBIG( idrhofl, nptr, wkrd )
c        ... and up-spin delta-density:
        call READBIG( idrhsfl, nptr, grho )
c        ... and combine to get total dn-spin first derivatives:
        do  nr=1,nptr
          drho(nr,iu) = drho(nr,iu) + wkrd(nr) - grho(nr)
        enddo
      enddo
c
c Retrieve double (XX,YY,ZZ) derivs of reference and scf density ...
      do  iu=1,3
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhofl, nptr, wkrd )
        call READBIG( idrhsfl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr)/two + wkrd(nr) - grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          g2rho(nr) = g2rho(nr) + ddrho(nr)
          gpggp(nr) = gpggp(nr) + drho(nr,iu)*drho(nr,iu)*ddrho(nr)
        enddo
      enddo
c
c Retrieve 2nd (XY,YZ,ZX) derivs of reference and scf density ...
      do  iu=1,3
        iv = iu + 1
        if( iv.eq.4 ) iv = 1
        call READBIG( idrhfl0, nptr, ddrho )
        call READBIG( idrhofl, nptr, wkrd )
        call READBIG( idrhsfl, nptr, grho )
c        ... combine
        do  nr=1,nptr
          ddrho(nr) = ddrho(nr)/two + wkrd(nr) - grho(nr)
        enddo
c        ... and accumulate gga gradient terms:
        do  nr=1,nptr
          gpggp(nr) = gpggp(nr) + two*drho(nr,iu)*drho(nr,iv)*ddrho(nr)
        enddo
      enddo
c
c Assemble | grad[rho] | term in grho(), now that scratch work done:
      do  nr=1,nptr
        grho(nr) = SQRT( drho(nr,1)**2 + drho(nr,2)**2 + drho(nr,3)**2 )
      enddo
c
c And finish the gpggp term, with conditional division by grho:
      do  nr=1,nptr
        if( grho(nr).gt.vsmall )then
          gpggp(nr) = gpggp(nr) / grho(nr)
        else
          gpggp(nr) = zero
        endif
      enddo
c
c    That's all Folks!
c
      RETURN
      END
