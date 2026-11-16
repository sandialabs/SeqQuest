c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c vdW control module
c Purpose: Manages data about the nature of vdW corrections
c Written: Peter A. Schultz, Febraury-2013, for 2.65
c
c Module data:
c  i_vdw, i_<vdW> = index of activated vdW, and indices of available vdW corrections
c  vdw_s          = s-factor, for global (fcnal-dependent) scaling of vdW correction
c  vdw_d          = d-factor, parameter defining the long-to-short interpolation
c  vdw_rc         = R-cutoff, the range at which vdW correction is cutoff
c                   (NB: Quest short-circuits this to be the longest overlap interaction)
c Routines:
c  VDWSET - set the kind of (vdW) force field desired)
c  i-fcn IVDWTYPE() - return type (number) of ff
c                     NB: this is not clean module design, used for internal debugging
c  L-fcn DO_VDW()   - returns boolean for if vdW corrections are on or off
c  VDWGETNAME       - returns name of activated vdW correction
c  VDWINIT(iflag)   - initializes module (for iflag<0), sets up global vdW params
c  VDWSET_S/VDWSET_S   - get/set vdW s-factor (global fcnal scaling)
c  VDWSET_D/VDWSET_D   - get/set vdW d-factor (interpolation param)
c  VDWSET_RC/VDWSET_RC - get/set vdW R-cutoff (vdW, which Quest short-circuits
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWSET
c
c
      subroutine VDWSET( vdwlbl )
c---------------------------------------------------------------
c Purpose: set up flags for vdW correction
c Written: Peter A. Schultz, February-2013, for v2.65 (vdW)
c---------------------------------------------------------------
c
c  Sets i_vdw (communicated via common)
c
c  Default vdw is OFF
c
c  To add new vdw correction methods:
c    Extend COMMON /VDWCTRL / with new i_{labels}
c    Add new unique integer values for i_{labels" below
c    Add check for {label} in if-block here
c    Add code in vdw* routines to implement the new vdW types
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
c Input declarations:
      CHARACTER*(*)  vdwlbl
c Local declarations:
      CHARACTER*6  vdwtyp
c
c >>>> EXECUTABLE CODE:
c
c The following are vdW id tags.  They must be distinct from one another.
      i_dftd2 = 1
      i_ulg   = 2
c      i_newvdw = 3 i.e. extend existing sequence
c
c Now go looking for a matching vdw type:
c
      vdwtyp = vdwlbl(1:6)
c
      if(     vdwtyp .eq. 'OFF   ' .or. vdwtyp .eq. 'off   '
     $   .or. vdwtyp .eq. 'NONE  ' .or. vdwtyp .eq. 'none  ' )then
c
        i_vdw = 0
c
      elseif( vdwtyp .eq. 'DFT-D2' .or. vdwtyp .eq. 'GRIMME'
     $   .or. vdwtyp .eq. 'DFTD2 ' .or. vdwtyp .eq. 'Grimme'
     $   .or. vdwtyp .eq. 'D2    '                           )then
c
        i_vdw = i_dftd2
c
      elseif( vdwtyp .eq. 'ULG   ' .or. vdwtyp .eq. 'ulg   '
     $   .or. vdwtyp .eq. 'ULG12 ' .or. vdwtyp .eq. 'ulg12 ' )then
c
        i_vdw = i_ulg
c
      else
        i_vdw = -1
        call STOPXERR( 'this vdw option is not implemented' )
      endif
c
      iflag = i_vdw
      call VDWINIT( iflag )
c
c    That's all Folks!
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IVDWTYPE
c
      INTEGER FUNCTION IVDWTYPE( )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      IVDWTYPE = i_vdw
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DO_VDW
c
      LOGICAL FUNCTION DO_VDW()
c
c DO NOT PUT (DEBUG) WRITE INTO THIS FUNCTION, CODE FAILS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      if( i_vdw .gt. 0 )then
        DO_VDW = .true.
      else
        DO_VDW = .false.
      endif
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWGETNAME
c
      subroutine VDWGETNAME( vdwtyp )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
c Local declarations:
      CHARACTER*6  vdwtyp
c
      if( i_vdw .eq. 0 )then
c
        vdwtyp = 'NONE  '
c
      elseif( i_vdw .eq. i_dftd2 )then
c
        vdwtyp = 'DFT-D2'
c
      elseif( i_vdw .eq. i_ulg )then
c
        vdwtyp = 'ULG   '
c
      else
c
        vdwtyp = 'UNKNOW'
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWINIT
c
      subroutine VDWINIT( iflag )
c---------------------------------------------------------------
c Purpose: set up global parameters for vdW corrections
c---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
c Simple initialize:
      if( iflag.lt.0 ) i_vdw = 0
      if( i_vdw .eq. 0 ) RETURN
c
c Check that we have a valid initialization (SET->INIT)
      if( i_vdw .ne. iflag ) call STOPXERR( 'bad VDWINIT call' )
c
      vdw_rc = 60.d0
      if( i_vdw .eq. i_dftd2 )then
c       Set s-factor (0.75/pbe, 1.2/blyp) & d-factor for DFT-D2
        vdw_s = 0.75
        vdw_d = 20.
      elseif( i_vdw .eq. i_ulg )then
c       Set s-factor (0.7012 assumes PBE) & d-factor (aka b_{lg}) for ULG
        vdw_s = 0.7012
        vdw_d = 0.6966
      else
        call STOPXERR( 'VDWINIT: unknown i_vdw' )
      endif
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWGET_S
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWSET_S
c
      subroutine VDWSET_S( svdw )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      vdw_s = svdw
c
      RETURN
c
      entry VDWGET_S( svdw )
c
      if( i_vdw.gt.0 ) svdw = vdw_s
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWGET_D
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWSET_D
c
      subroutine VDWSET_D( dvdw )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      vdw_d = dvdw
c
      RETURN
c
      entry VDWGET_D( dvdw )
c
      if( i_vdw.gt.0 ) dvdw = vdw_d
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWGET_RC
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWSET_RC
c
      subroutine VDWSET_RC( rcvdw )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      vdw_s = rcvdw
c
      RETURN
c
      entry VDWGET_RC( rcvdw )
c
      if( i_vdw.gt.0 ) rcvdw = vdw_rc
c
      RETURN
      END
