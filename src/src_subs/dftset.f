c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTSET
c
c
      subroutine DFTSET( dftlbl )
c---------------------------------------------------------------
c Purpose: set up for flags for dft of particular flavor
c
c Written: Peter A. Schultz, December-2000, for v2.46 (gga)
c
c Revision history:
c  13Nov12-JEM/2.65: HSE added
C  17Spt08-RJM/2.62: added AM05 spin 
c  24May06-AEM/2.60: AM05 added      
c  25Jun03-AEM/2.56: BLYPSP added
c  10Jun03-AEM/2.55: BLYP added
c   7Aug01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs, writes
c  25Mar01-PAS/2.46a: make gga default be PBE
c---------------------------------------------------------------
c
c  Set i_fcnal (communicated via common) for a given dftlbl,
c  and set whether gga/spin with i_gga/i_spin
c
c  Default lda is capz=ceperley/alder param'd by perdew/zunger
c  Default gga is pbe
c
c  To add new functionals, extend COMMON/DFTTYPS/ with new
c  i_{labels}, add new unique integer values below, add
c  check for label in if-block here, and add code in
c  routines "dft*gga"/"dft*lda" for new gga/lda functionals.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /DFTTYPS/ i_fcnal, i_gga,i_spin,
     $ i_capz,i_pw91,i_pbe,i_blyp,i_am05,i_hse,i_hse2
c
c input declarations:
      CHARACTER*(*)  dftlbl
c local declarations:
      CHARACTER*6  fcnal
c
c >>>> EXECUTABLE CODE:
c
c The following are fcnal id tags.  They must be distinct
c from one another.
      i_capz = 1
      i_pw91 = 2
      i_pbe  = 3
      i_blyp = 4
      i_am05 = 5
      i_hse  = 6
      i_hse2 = 7
c      i_newdft = 6 i.e. extend existing sequence
c
c Now go looking for a matching fcnal:
c
      fcnal = dftlbl(1:6)
c
      if(     fcnal .eq. 'capz  ' .or. fcnal .eq. 'CAPZ  '
     $   .or. fcnal .eq. 'lda   ' .or. fcnal .eq. 'LDA   ' )then
c
        i_fcnal = i_capz
        i_gga = 0
        i_spin = 0
c
      elseif( fcnal .eq. 'capzsp' .or. fcnal .eq. 'CAPZSP'
     $   .or. fcnal .eq. 'lsd   ' .or. fcnal .eq. 'LSD   '
     $   .or. fcnal .eq. 'lda-sp' .or. fcnal .eq. 'LDA-SP' )then
c
        i_fcnal = i_capz
        i_gga = 0
        i_spin = 1
c
      elseif( fcnal .eq. 'pw91  ' .or. fcnal .eq. 'PW91  ' )then
c
        i_fcnal = i_pw91
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'pw91sp' .or. fcnal .eq. 'PW91SP' )then
c
        i_fcnal = i_pw91
        i_gga = 1
        i_spin = 1
c
      elseif( fcnal .eq. 'pbe   ' .or. fcnal .eq. 'PBE   '
     $   .or. fcnal .eq. 'gga   ' .or. fcnal .eq. 'GGA   ' )then
c
        i_fcnal = i_pbe
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'pbe-sp' .or. fcnal .eq. 'PBE-SP'
     $   .or. fcnal .eq. 'gga-sp' .or. fcnal .eq. 'GGA-SP' )then
c
        i_fcnal = i_pbe
        i_gga = 1
        i_spin = 1
c
      elseif( fcnal .eq. 'blyp  ' .or. fcnal .eq. 'BLYP  ' )then
c
        i_fcnal = i_blyp
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'blypsp' .or. fcnal .eq. 'BLYPSP' )then
c
        i_fcnal = i_blyp
        i_gga = 1
        i_spin = 1
c
      elseif( fcnal .eq. 'am05  ' .or. fcnal .eq. 'AM05  ' )then
c
        i_fcnal = i_am05
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'am05sp' .or.  fcnal .eq. 'AM05SP' )then
c
        i_fcnal = i_am05
        i_gga = 1
        i_spin = 1
c
      elseif( fcnal .eq. 'hse   ' .or. fcnal .eq. 'HSE   ' )then
c
        i_fcnal = i_hse
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'hse-sp' .or. fcnal .eq. 'HSE-SP' )then
c
        i_fcnal = i_hse
        i_gga = 1
        i_spin = 1
c
      elseif( fcnal .eq. 'hse2  ' .or. fcnal .eq. 'HSE2  ' )then
c
        i_fcnal = i_hse2
        i_gga = 1
        i_spin = 0
c
      elseif( fcnal .eq. 'hse2sp' .or. fcnal .eq. 'HSE2SP' )then
c
        i_fcnal = i_hse2
        i_gga = 1
        i_spin = 1
c
      else
c
        call STOPXERR( 'no fcnal/this dft flavor not available' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
