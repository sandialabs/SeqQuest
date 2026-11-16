c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BULKPOL
c
c
      subroutine BULKPOL( IWR, ndim, elchrg,dielec, rprim, engypol )
c---------------------------------------------------------------
c Purpose: generate bulk polarization energy in response to charge
c          in cell.
c
c Written: Peter A. Schultz, 2-October-2003, for v2.57
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Thanks to Art Edwards (AFRL) for finding the method used here.
c
c Notes:
c  The model used here for the bulk polarization in response to
c  a local charged defect is the Jost model:
c     W. Jost, J. Chem. Phys 1, 466 (1933).
c  where the polarization energy is approximated by the polarization
c  response of an infinite bulk dielectric material outside a
c  spherical hollow of the same volume as the simulation supercell.
c  The formula for the bulk polarization response (in Hartree):
c    E(Jost) = ( 1 - 1/dielec ) * ( Q^2 / 2*R(Jost) )
c  where:
c    E = output bulk polarization energy (Hartree).
c    dielec = the static dielectric constant of the material
c    Q = net charge in the simulation cell (au)
c    R = is the radius of a sphere that gives the same volume as
c        the parallelpiped supercell (bohr)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, three=3.d0 )
c
c Input:
c   rprim - supercell vectors Ry/bohr)
      DIMENSION  rprim(3,*)
c
c Local:
      DIMENSION  vtmp(3)
      DATA  toev / 13.605805d0 /
      DATA  pi / 3.14159265358979d0 /
      DATA  fourthrd / 1.333333333333333d0 /
c
c >>>> EXECUTABLE CODE:
c
      engypol = zero
c If not bulk, with a net charge, and with a dielectric, skip ...
      if( ndim.ne.3 .or. ABS( elchrg ).eq.zero .or.
     $    dielec.le.zero ) RETURN
c
c Obtain the volume of the supercell
      call CROSS( rprim(1,2),rprim(1,3), vtmp )
      volum = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $           + rprim(3,1)*vtmp(3) )
c
c Compute the radius of a sphere with the same volume
      r3 = volum / ( fourthrd*pi )
      rln = LOG( r3 ) / three
      rjost = EXP( rln )
c
c The Jost model bulk polarization energy around spherical cavity
      ejost = ( one - one/dielec ) * elchrg**2 / (two*rjost)
c
c Convert to Rydberg
      engypol = two*ejost
c
c Tell the user about it
 9020 format(1x,a,3f16.6)
      if( IWR .gt. 0 )then
        write(IWR,9020) 'BULKPOL: bulk polarization response'
        write(IWR,9020) ' Jost model: Epol = (1-1/e)*q^2 / (2*Rjost)'
        write(IWR,9020) ' Dielectric constant, e=', dielec
        write(IWR,9020) ' cell volume(bohr^3), Rjost(bohr)=',
     $                    volum,  rjost
        engyev = engypol*toev
        write(IWR,9020) ' Bulk polarization ENERGY (Ry)=',engypol
        write(IWR,9020) ' Bulk polarization ENERGY (eV)=',engyev
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
