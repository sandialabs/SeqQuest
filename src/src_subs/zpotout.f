c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ZPOTOUT
c
c
      subroutine ZPOTOUT( igridfl,IWRT,
     $ iatmfmt, ntyp,natm, natmnm, itypa, atmnm,typnm, ratm,
     $ nspin, hh, n1r,n2r,n3r,nptr,
     $ wkz, delrho, delves )
c---------------------------------------------------------------
c Purpose: output plane(z)-averaged potential and density
c
c Written: P.A. Schultz, 2011, for qDFT-flDFT coupling prep in 2.62
c          12-january-2011
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c On input:
c   igridfl > 0 - read density and potential from disk file
c           else, use density and potential as passed
c   IWRT    > 0 - write analysis to file unit = IWRT
c           else, do not write analysis
c   hh(3,3)     = basis vectors for real-space grid
c   n1r,n2r,n3r = linear dimensions of grid
c   nptr        = total number of grid points (=n1r*n2r*n3r)
c   wkz()    = scratch workspace, of at least wkz(n3r) size
c   delrho() = delta density on grid, implicit dims(n1r,n2r,n3r)
c   delves() = electrostatic potential, implicit dims(n1r,n2r,n3r)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (zero=0.d0, one=1.d0)
c
c Input arrays:
c   ... atomic positions
      DIMENSION  itypa(*)
      CHARACTER*(*)  typnm, atmnm
      DIMENSION  typnm(*), atmnm(*)
      DIMENSION  ratm(3,natm)
c   ... potentials
      DIMENSION  hh(3,3)
      DIMENSION  delrho(nptr), delves(nptr)
c
c Scratch arrays:
      DIMENSION  wkz(*)
c
c Local declarations:
      DIMENSION  ave(2)
c
c >>>> EXECUTABLE CODE:
c
      if( IWRT .le. 0 ) RETURN
c
c  If input not loaded, get the density and potential from disk file
c
      if( igridfl .ge. 0 )then
c       Retrieve the density and potential from the grid file
        REWIND( unit=igridfl )
c
c       Retrieve converged grid delta density from igridfl
c
c       Read delta-rho (total, or up-spin if spin-polarized)
        call READBIG( igridfl, nptr, delrho )
        if( nspin .eq. 2 )then
c         Read dn-spin delta-rho and sum into up-delrho to get total delta-rho
          call READBIG( igridfl, nptr, delves )
          call DAXPY( nptr, one, delves,1, delrho,1 )
        endif
c
c       Read the delta-electrostatic potential
        call READBIG( igridfl, nptr, delves )
      endif
c
c  Do the analysis and send it to output
c
      if( IWRT .ge. 0 )then
c       Output the density to the specified file:
c
        write(IWRT,*) '### ATOM POSITIONS: name, type, x,y,z (Bohr) ###'
        call CONFIGWR( IWRT   , iatmfmt, natmnm, ntyp,typnm,
     $   natm, atmnm, itypa, ratm )
c
        write(IWRT,*) '##### Z-AVERAGE DENSITY and POTENTIAL #####'
        write(IWRT,'(a/f10.7,a)') 'spacing of points',hh(3,3),' Bohr'
        write(IWRT,'(a/i6)') 'number of points in z-direction',n3r
        call GRIDAVE( IWRT,
     $   n1r,n2r,n3r, 1,delrho,ave, wkz, 'del-density (e/bohr^3)')
        call GRIDAVE( IWRT,
     $   n1r,n2r,n3r, 1,delves,ave, wkz, 'electrostatic del-Ves (Ry)')
      endif
c
      write(IWRT,*) '##### END Z-AVERAGE DENSITY and POTENTIAL #####'
c
c    That's all Folks!
c
      RETURN
      END
