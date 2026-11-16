c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGET
c
c
c---------------------------------------------------------------
c Purpose: access information about the nature of the dft
c
c Written: Peter A. Schultz, 7-August-2001, for v2.49 (spin)
c          based on a suggestion by Richard P. Muller (Caltech)
c
c Revision history:
c  13Nov12-JEM/2.65: HSE flags added
c  24May06-AEM/2.60: AM05 added
c  11Jun03-PAS/2.55: add query about Laplacian
c   7Aug01-PAS/2.49: adapted from implied Muller routines
c---------------------------------------------------------------
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGETGGA
c
c
      subroutine DFTGETGGA( do_gga )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /DFTTYPS/ i_fcnal, i_gga,i_spin,
     $ i_capz,i_pw91,i_pbe,i_blyp,i_am05,i_hse,i_hse2
c
c output declaration:
      LOGICAL  do_gga, do_spin, do_lap
      CHARACTER*6  fcnal
c
c >>>> EXECUTABLE CODE:
c
      if( i_gga .eq. 0 )then
        do_gga = .false.
      elseif( i_gga .eq. 1 )then
        do_gga = .true.
      else
        call STOPXERR( 'dftget - bad record of gga' )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGETLAP
c
c
      entry DFTGETLAP( do_lap )
c
      do_lap = .false.
c Right now, the only laplacian-based functional is BLYP:
      if( i_fcnal .eq. i_blyp )then
        do_lap = .true.
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGETSPIN
c
c
      entry DFTGETSPIN( do_spin )
c
      if( i_spin .eq. 0 )then
        do_spin = .false.
      elseif( i_spin .eq. 1 )then
        do_spin = .true.
      else
        call STOPXERR( 'dftget - bad record of spin' )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGETTYP
c
c
      entry DFTGETTYP( fcnal )
c
      if( i_fcnal .eq. 0 )then
c
        fcnal = 'NONE  '
c
      elseif( i_fcnal .eq. i_capz )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'CAPZ  '
        else
          fcnal = 'CAPZSP'
        endif
c
      elseif( i_fcnal .eq. i_pw91 )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'PW91  '
        else
          fcnal = 'PW91SP'
        endif
c
      elseif( i_fcnal .eq. i_pbe )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'PBE   '
        else
          fcnal = 'PBE-SP'
        endif
c
      elseif( i_fcnal .eq. i_blyp )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'BLYP  '
        else
          fcnal = 'BLYPSP'
        endif
c
      elseif( i_fcnal .eq. i_am05 )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'AM05  '
        else
          fcnal = 'AM05SP'      
        endif  
c
      elseif( i_fcnal .eq. i_hse )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'HSE   '
        else
          fcnal = 'HSE-SP'      
        endif  
c
      elseif( i_fcnal .eq. i_hse2 )then
c
        if( i_spin .eq. 0 )then
          fcnal = 'HSE2  '
        else
          fcnal = 'HSE2SP'      
        endif  
c
      else
c
        fcnal = 'UNKNOW'
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
