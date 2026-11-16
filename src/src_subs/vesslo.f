c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VESSLO
c
c
      subroutine VESSLO( IWR, ndim, n1r,n2r,n3r,nptr, weight,
     $ ws1,ws2,ws3,ftarray,
     $ rhoslo, espot, gvecsq )
c---------------------------------------------------------------
c Purpose: compute real-space es potential from slow density
c
c Written: P.A. Schultz, based on PJF's old routines
c
c Revision history:
c  26Jan11-PAS/2.62: tweaks for output (for solvation coupling)
c  30Mar99-PAS/2.32: additional diagnostics, change b.c.
c   9Dec98-PAS/2.29: promote source to explicit double precision
c  22Oct98-PAS/2.27: boundary condition options
c  23Jun98-PAS/2.23: extend graded output controls
c---------------------------------------------------------------
c
cdev262: bulk es pot offset
c To-do: should offset bulk es potential rather than by zero,
c   as is currently being done for ibndopt = 1.
c   The fft library does this automatically, but best to make
c   this explicit here rather than count on outside function.
c
c Atomic units: distance in bohr, energy in Rydberg
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input arrays: density(unchanged), G-vectors(destroyed)
      DIMENSION  rhoslo(nptr), gvecsq(nptr)
c Output arrays: electrostatic potential
      DIMENSION  espot(2,nptr)
c
c Scratch arrays:
      DIMENSION  ws1(*),ws2(*),ws3(*), ftarray(2,*)
c
c Local declarations:
      DIMENSION  ave(2)
c
c Set default for type of boundary condition
c  For ibndopt = 0, use fft potential unshifted
c              = 1, use z-edge potential as (non-bulk) vacuum zero
c              = 2, use average pot as zero
      DATA  ibndopt / 1 /
      DATA  zero, eight / 0.d0, 8.d0 /
      DATA  lstout / 0 /
      DATA  IWRT / 0 /
c
c >>>> EXECUTABLE CODE:
c
      eightpi = eight*pi
c
      if( lstout.gt.3 ) IWRT = IWR
      if( lstout.gt.3 )
     $ call GRIDAVE(IWR,n1r,n2r,n3r,1,rhoslo,ave,ftarray,'slow density')
c
c Load (real*8) density into (complex*16) workspace for fft
      call MKZERO( 2*nptr, espot )
      call DCOPY( nptr, rhoslo,1, espot,2 )
c
c  Fourier transform density:
c
      call CFFT3D( espot, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c
c  Set rho(g=0), ie. net charge, to zero:
c
      indexg0 = (n1r+2)/2 + n1r*(n2r/2) + n1r*n2r*(n3r/2)
      defect = weight*DBLE(nptr)*espot(1,indexg0)
c
      if( lstout.gt.0 .and. IWR.ge.0 )then
        if( lstout.gt.2 )
     $  write(IWR,900) indexg0, espot(1,indexg0),espot(2,indexg0)
  900   format(/1x,'indexg0=',i12,'; rho(g=0)=',1p2d18.8)
c
        write(IWR,910) 'delta-charge defect (electrons)=',defect
      endif
  910 format(1x,a,1p2d18.8)
c
c  Compute Coulomb potential (screening eliminated)
c
      espot(1,indexg0) = zero
      espot(2,indexg0) = zero
      ave1 = zero
      ave2 = zero
c Set G=0 to small number to prevent divide-by-zero error in loop
      gvecsq(indexg0) = 1.d-15
      do  ng=1,nptr
        espot(1,ng) = eightpi*espot(1,ng) / gvecsq(ng)
        espot(2,ng) = eightpi*espot(2,ng) / gvecsq(ng)
      enddo
c
c  Fourier transform es potential back to coordinate space
c
      call CFFT3DB( espot, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c
      call GRIDAVE( IWRT, n1r,n2r,n3r,2,espot,ave,ftarray, 'slow Ves' )
c
      if( lstout.ge.1 .and. IWR.ge.0 )then
c       Print out boundary values:
        write(IWR,910) 'global average fft es-potential(Ry)=', ave(1)
        if( ndim .lt. 3)
     $  write(IWR,910) 'z-edge average fft es-potential(Ry)=',
     $                  ftarray(1,1)
      endif
c
c  Shift electrostatic potential (e.g. to rescale vacuum to zero)
c
      if( ibndopt.gt.0 )then
c
        if( ibndopt.eq.1 )then
          if( ndim.lt.3 )then
c           Use z=zmax planar edge average as vacuum zero
            n3vac = 1
C            n3vac = n3r
C            if( (n3r/2).eq.((n3r+1)/2) ) n3vac = 1
            ave1 = ftarray(1,n3vac)
            ave2 = ftarray(2,n3vac)
          else
cpas        Bulk 3D, no boundary, use raw electrostatic potential
cpas         (dev262 - should be altered to global-average-offset)
            ave1 = zero
            ave2 = zero
          endif
        else
c         Offset potential by global average, new global average is zero
          ave1 = ave(1)
          ave2 = ave(2)
        endif
c
        do  nr=1,nptr
          espot(1,nr) = espot(1,nr) - ave1
          espot(2,nr) = espot(2,nr) - ave2
        enddo
c
c       Leave record of vacuum shift, will be very important for tracking b.c.
        if( lstout.ge.1 .and. IWR.ge.0 )then
          write(IWR,*)
          write(IWR,910) 'Vacuum shift for fft es-potential(Ry)=', ave1
        endif
c
      endif
c
c  Compact the es pot array to keep only the real part
C      iq = 2
C      id = -1
C      ir = 0
C      do  nr=1,nptr
C        ir = ir + iq - 1
C        iq = iq + id
C        id = -id
C        espot(iq,ir) = espot(1,nr)
C      enddo
      call DCOPY( nptr, espot,2, gvecsq,1 )
      call DCOPY( nptr, gvecsq,1, espot,1 )
c
c    That's all Folks!
c
      RETURN
      END
