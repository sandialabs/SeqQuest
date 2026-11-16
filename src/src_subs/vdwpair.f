c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWCROSS
c
c
      subroutine VDWCROSS( ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
c---------------------------------------------------------------
c Purpose: interpolated pairwise vdW terms from atomic terms
c Written: P.A. Schultz, 8-February-2013, for 2.65
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Inputs:
      DIMENSION  c6vdw(ntyp),r0vdw(ntyp)
c Outputs:
      DIMENSION  c6ijvdw(ntyp,ntyp), r0ijvdw(ntyp,ntyp)
c
      do 200 ityp=1,ntyp
        c6i = c6vdw(ityp)
        r0i = r0vdw(ityp)
        do 100 jtyp=ityp,ntyp
          c6j = c6vdw(jtyp)
          r0j = r0vdw(jtyp)
          call VDWIJMIX( c6i,c6j,c6ij, r0i,r0j,r0ij )
          c6ijvdw(ityp,jtyp) = c6ij
          r0ijvdw(ityp,jtyp) = r0ij
          c6ijvdw(jtyp,ityp) = c6ij
          r0ijvdw(jtyp,ityp) = r0ij
  100   continue
  200 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWIJMIX
c
c
      subroutine VDWIJMIX( c6i,c6j,c6ij, r0i,r0j,r0ij )
c---------------------------------------------------------------
c Purpose: interpolated pairwise vdW terms from atomic terms
c Written: P.A. Schultz, 8-February-2013, for 2.65
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      REAL*8       zero,one,two
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0 )
c Inputs:
      REAL*8   c6i,c6j, r0i,r0j
c Outputs:
      REAL*8   c6ij, r0ij
Cc
C      INTEGER  i_dftd2,i_ulg
C      DATA     i_dftd2,i_ulg / 1,2 /
c
      c6ij = zero
      r0ij = one
      if( c6i.gt.zero .and. c6j.gt.zero )then
c       Both atoms have valid vdW forms
c
c       All vdW (so far) use root norm for c6 coeffs:
        c6ij = SQRT( c6i*c6j )
c
        if( i_vdw .eq. i_dftd2 )then
c         The pair interaction r0 is sum of atomic radii
          r0ij = r0i + r0j
        elseif( i_vdw .eq. i_ulg )then
c         The r0 is root norm of atomic r0
          r0ij = SQRT( r0i*r0j )
        else
c         Invalid vdW option
          r0ij = SQRT( r0i*r0j )
          call STOPXERR( 'vdwij: invalid i_vdw' )
        endif
c
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWPAIR
c
c
      subroutine VDWPAIR( e, de, r, c6, r0 )
c---------------------------------------------------------------
c Purpose: compute energy and gradient of vdW pair potential
c Written: P.A. Schultz, 6-February-2013, for 2.65
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      REAL*8       zero, one, two, six
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, six=6.d0 )
c
c Output:
      REAL*8   e, de
c Input:
      REAL*8   r, c6, r0
c
c Local:
      REAL*8  r5,r6,r7,expfac,f,df
c
c >>>> EXECUTABLE CODE
c
      d_vdw = vdw_d
      e = zero
      de = zero
      if( r .lt. 0.1 .or. i_vdw .eq. 0 ) RETURN
c
      r5 = r**5
      r6 = r5*r
      r7 = r6*r
c
      if( i_vdw .eq. i_dftd2 )then
c       Grimme06 DFT-D2
        expfac = EXP( -d_vdw*( r/r0 - one ) )
        f = one / (one + expfac)
        e = -(c6/r6) * f
c
        df = (d_vdw/r0)*expfac * f**2
        de = -(c6/r6) * df + 6.d0*c6/r7 * f
c
      elseif( i_vdw .eq. i_ulg )then
c       Kim/Choi/WAG ULG12
        e = -c6 / (r6 + d_vdw*r0**6)
        de = c6*(6.d0*r5) / (r6 + d_vdw*r0**6)**2
c
      else
c       Not a valid vdw option
        call STOPXERR( 'vdwpair: bad vdw option' )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWSTAMP
c
c
      subroutine VDWSTAMP( IWR, ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
c---------------------------------------------------------------
c Purpose: interpolated pairwise vdW terms from atomic terms
c Written: P.A. Schultz, 8-February-2013, for 2.65
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Inputs:
      DIMENSION  c6vdw(ntyp),r0vdw(ntyp)
c Outputs:
      DIMENSION  c6ijvdw(ntyp,ntyp), r0ijvdw(ntyp,ntyp)
c
      if( IWR.le.0 ) RETURN
c
      mtyp = ntyp
 9100 format(2x,a)
 9101 format(4x,1p8d12.4)
c
      write(IWR,9100) ' ***** Recap of vdW parameters *****'
c
      call VDWGET_S( sfac )
      call VDWGET_D( dfac )
      call VDWGET_RC( rcfac )
      write(IWR,'(3x,a,3f12.6)') 'sfactor,dfactor,rc=',sfac,dfac,rcfac
c
      write(IWR,9100) ' C6atom (Ry):'
      write(IWR,9101) (c6vdw(ityp),ityp=1,mtyp)
c
      write(IWR,9100) ' R0atom (Bohr):'
      write(IWR,9101) (r0vdw(ityp),ityp=1,mtyp)
c
      write(IWR,9100) ' C6pair matrix[type,type] (Ry):'
      do  ityp=1,mtyp
        write(IWR,9101) (c6ijvdw(ityp,jtyp),jtyp=1,mtyp)
      enddo
c
      write(IWR,9100) ' R0pair matrix[type,type] (Bohr):'
      do  ityp=1,mtyp
        write(IWR,9101) (r0ijvdw(ityp,jtyp),jtyp=1,mtyp)
      enddo
c
      write(IWR,9100) ' <<<<< End vdW parameter recap'
c
      RETURN
      END
