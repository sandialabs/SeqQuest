c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTSTR2
c
c
      subroutine DFTSTR2( idrhfl0,idrhofl,idrhsfl,
     $ ndim, strlap, nptr,weight,
     $ dxc2, ddrho0,dddrho )
c---------------------------------------------------------------
c Purpose: compute, on grid, integrals needed for non-spin gga
c          stress calculation from 2nd derivatives.
c
c Written: Peter A. Schultz, 11-June-2003, for 2.55 (BLYP)
c
c Revision history:
c  10Oct03-PAS/2.57: accommodate spin-dependent stress for GGA
c---------------------------------------------------------------
c
c Notes:
c   Evaluate the grid integral for non-spin GGA stress evaulation
c   coming from the Laplacian term (e.g. BLYP):
c
c     stress(i,j) = INT( dxc2(r) * d2_rh(r)/(dx_i*dx_j) )
c
c  On disk(s) we have all first and second derivatives of rho:
c    d/dx,d/dy,d/dz
c    d2/dxdx, d2/dydy, d2/dzdz, d2/dxdy, d2/dydz, d2/dzdx
c
c Variables:
c    idrhfl0 - file with derivatives from reference density
c    idrhofl - file with derivatives from delta (scf) density
c    nptr    - size of grid field
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0 )
c
c Input declarations:
      DIMENSION  dxc2(nptr)
c Output stuff:
      DIMENSION  strlap(3,3)
c Scratch:
      DIMENSION  ddrho0(nptr), ddrho(nptr)
c
c >>>> EXECUTABLE CODE:
c
      if( ndim .eq. 0 ) RETURN
c
c  Diagonal terms
c
      do  iu=1,3
        if( iu .le. ndim )then
c         Retrieve double (XX,YY,ZZ) derivs of reference ...
          call READBIG( idrhfl0, nptr, ddrho0 )
c          ... and scf delta-density
          call READBIG( idrhofl, nptr, ddrho )
c          ... combine
          do  nr=1,nptr
            ddrho(nr) = ddrho(nr) + ddrho0(nr)
          enddo
c          ... and evaluate integral
          strii = zero
          do  nr=1,nptr
            strii = strii + dxc2(nr)*ddrho(nr)
          enddo
          strlap(iu,iu) = weight*strii
        else
c         Skip thru to position files at beginning of off-diags
          call READSKP( idrhfl0, nptr )
          call READSKP( idrhofl, nptr )
        endif
      enddo
c
c  Off-diagonal terms
c
      do  iu=1,3
        iv = iu + 1
        if( iv .eq. 4 ) iv = 1
        if( MAX(iu,iv) .le. ndim )then
c         Retrieve 2nd (XY,YZ,ZX) derivs of reference ...
          call READBIG( idrhfl0, nptr, ddrho0 )
c          ... and scf delta-density
          call READBIG( idrhofl, nptr, ddrho )
c          ... combine
          do  nr=1,nptr
            ddrho(nr) = ddrho(nr) + ddrho0(nr)
          enddo
c          ... and evaluate integral
          strij = zero
          do  nr=1,nptr
            strij = strij + dxc2(nr)*ddrho(nr)
          enddo
          strlap(iu,iv) = weight*strij
          strlap(iv,iu) = weight*strij
C        else
C          call READSKP( idrhfl0, nptr )
C          call READSKP( idrhofl, nptr )
        endif
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
