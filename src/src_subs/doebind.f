c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DOEBIND
c
c
      subroutine DOEBIND( IWR, engytotl, engyref, engybind,
     $ toev, natm, itypa, atengy )
c---------------------------------------------------------------
c Purpose: compute binding energy for system wrt reference atom
c
c Written: Peter A. Schultz, 3-October-2001, for v2.50
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  itypa(*)
      DIMENSION  atengy(natm)
c
c local declarations:
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      engyref = zero
      do  iatm=1,natm
        ityp = itypa(iatm)
        engyref = engyref + atengy(ityp)
      enddo
c
      engybind = engytotl - engyref
      if( engyref .ne. zero )then
        engybev = toev*engybind
        write(IWR,'(1x,a,f20.10,a,f20.4,a)')
     $   'Binding energy =',engybind,' Ry =',engybev,' eV'
      endif
c
c    That's all Folks!
c
      RETURN
      END
