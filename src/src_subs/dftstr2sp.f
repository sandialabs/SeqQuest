c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTSTR2
c
c
      subroutine DFTSTR2SP( idrhfl0,idrhofl,idrhsfl,
     $     ndim, strlap, nptr,weight,
     $     dxc2up,  dxc2dn,   ddrho0,    ddrhoa,    ddrhob )
c---------------------------------------------------------------
c Purpose: compute, on grid, integrals needed for spin gga
c          stress calculation from 2nd derivatives.
c
c Written: Peter A. Schultz, 10-October-2003, for 2.57
c
c Revision history:
c  10Oct03-PAS/2.57: modeled off DFTSTR2
c---------------------------------------------------------------
c
c Notes:
c   Evaluate the grid integrals for spin-GGA stress evaluation
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
c    idrhsfl - file with derivatives from up-spin del-density
c    nptr    - size of grid field
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, half=0.5d0 )
c
c Input declarations:
      DIMENSION  dxc2up(nptr), dxc2dn(nptr)
c Output stuff:
      DIMENSION  strlap(3,3)
c Scratch:
      DIMENSION  ddrho0(nptr), ddrhoa(nptr), ddrhob(nptr)
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
          call READBIG( idrhofl, nptr, ddrhob )
c          ... and spin-up del-density
          call READBIG( idrhsfl, nptr, ddrhoa )
c          ... construct dd-del-dn = dd-rho - dd-del-up
          do  nr=1,nptr
            ddrhob(nr) = ddrhob(nr) - ddrhoa(nr)
          enddo
c          ... construct dd-up&dn = dd-del-up&dn + dd-rho0/2
          do  nr=1,nptr
            ddrhoa(nr) = ddrhoa(nr) + half*ddrho0(nr)
            ddrhob(nr) = ddrhob(nr) + half*ddrho0(nr)
          enddo
c          ... and evaluate integrals
          striia = zero
          striib = zero
          do  nr=1,nptr
            striia = striia + dxc2up(nr)*ddrhoa(nr)
            striib = striib + dxc2dn(nr)*ddrhob(nr)
          enddo
          strlap(iu,iu) = weight * ( striia + striib )
        else
c         Skip thru to position files at beginning of off-diags
          call READSKP( idrhfl0, nptr )
          call READSKP( idrhofl, nptr )
          call READSKP( idrhsfl, nptr )
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
          call READBIG( idrhofl, nptr, ddrhob )
c          ... and spin-up del-density
          call READBIG( idrhsfl, nptr, ddrhoa )
          do  nr=1,nptr
            ddrhob(nr) = ddrhob(nr) - ddrhoa(nr)
          enddo
          do  nr=1,nptr
            ddrhoa(nr) = ddrhoa(nr) + half*ddrho0(nr)
            ddrhob(nr) = ddrhob(nr) + half*ddrho0(nr)
          enddo
          strija = zero
          strijb = zero
          do  nr=1,nptr
            strija = strija + dxc2up(nr)*ddrhoa(nr)
            strijb = strijb + dxc2dn(nr)*ddrhob(nr)
          enddo
          strlap(iu,iv) = weight * ( strija + strijb )
          strlap(iv,iu) = weight * ( strija + strijb )
C        else
C          call READSKP( idrhfl0, nptr )
C          call READSKP( idrhofl, nptr )
C          call READSKP( idrhsfl, nptr )
        endif
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
