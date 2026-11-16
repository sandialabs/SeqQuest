c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FDEFECT
c
c
      subroutine FDEFECT( IWR, ifdefct, natm, frc,fdefct,
     $ nsyma,rmatsym,naofsym )
c---------------------------------------------------------------
c Purpose: symmetrize force defect (non-zero net drift)
c
c Written: Peter A. Schultz, 18-May-2001, for v2.47
c
c Revision history:
c  18May01-PAS/2.47: extract from main routine
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  frc(3,natm)
      DIMENSION  fdefct(3)
      DIMENSION  rmatsym(*)
      DIMENSION  naofsym(natm,*)
c local declarations:
      DIMENSION  vec0(3)
      DATA  vec0 / 0.d0,0.d0,0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( ifdefct.gt.0 .and. ifdefct.le.natm )then
c
c  Check that symmetry leaves defect atom fixed,
c  If not, then it is not symmetry preserving to install
c  a force defect, and we must skip this correction
        do  isym=1,nsyma
          jatm = naofsym(ifdefct,isym)
          if( jatm.ne.ifdefct ) RETURN
        enddo
c
c  Atom is fixed by symmetry, install force defect:
C        call SYMVEC( 1,fdefct, vec0, nsyma, rmatsym )
        write(IWR,'(a,3f12.6)') 'Force f-defect ->',fdefct
        frc(1,ifdefct) = frc(1,ifdefct) - fdefct(1)
        frc(2,ifdefct) = frc(2,ifdefct) - fdefct(2)
        frc(3,ifdefct) = frc(3,ifdefct) - fdefct(3)
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
