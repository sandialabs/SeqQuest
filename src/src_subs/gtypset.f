c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GTYPSET
c
c
      subroutine GTYPSET( glabel )
c---------------------------------------------------------------
c Purpose: set up flags for geometry update method
c
c Written: Peter A. Schultz, 30-October-2001, for v2.50 (dmd)
c
c Revision history:
c  21Jun07-PAS/2.60: bfgs disabled
c   7May04-PAS/2.58: add explicit quantum MD flag
c  19Dec02-PAS/2.55: add LBFGS option
c---------------------------------------------------------------
c
c  Setting ig_type (communicated via common) for a given glabel
c Relaxation schemes:
c  ig_steep = simple steepest descent
c  ig_broy  = modified Broyden (D.D. Johnson).
c  ig_md0   = simple dynamics
c  ig_dmd   = "damped dynamics" a la Jonsson - component-wise
c  ig_dmdat = "damped dynamics" a la Jonsson - atom-by-atom
c  ig_asd   = accelerated steepest descent (modified damped dyn.)
c  ig_lbfgs = limited memory BFGS algorithm
c  ig_cg    = conjugate gradient (not yet coded)
c  ig_qmd   = quantum molecular dynamics (RPM)
c
c  To add new methods, extend COMMON/GEOMTYPS/ with new
c  ig_{labels}, add new unique integer values below, add
c  check for label in if-block here, and add code in
c  routine "GRELAX" to implement the new scheme.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /GEOMTYPS/ ig_type, ig_steep,ig_broy,ig_md0,
     $ ig_dmd,ig_dmdat,ig_asd,ig_lbfgs,ig_cg,ig_qmd
c
c input declarations:
      CHARACTER*(*)  glabel
c local declarations:
      CHARACTER*6  gmeth
c
c >>>> EXECUTABLE CODE:
c
c The following are geometry method id tags.  They must be distinct
c from one another.
      ig_steep =  1
      ig_broy  =  2
      ig_md0   =  3
      ig_dmd   =  4
      ig_dmdat =  5
      ig_asd   =  6
      ig_lbfgs =  7
      ig_cg    =  8
      ig_qmd   =  9
c      ig_gmeth = 10 i.e. extend existing sequence
c
c Now go looking for a matching method:
c
      gmeth = glabel(1:6)
c
      if(     gmeth .eq. 'steepe' .or. gmeth .eq. 'STEEPE' )then
c       Steepest descent
        ig_type = ig_steep
c
      elseif( gmeth .eq. 'broyde' .or. gmeth .eq. 'BROYDE' )then
c       Modified Broyden scheme (D.D. Johnson)
        ig_type = ig_broy
c
      elseif( gmeth .eq. 'mdsimp' .or. gmeth .eq. 'MDSIMP' )then
c       Simple dynamics outside of MD engine
        ig_type = ig_md0
c
      elseif( gmeth .eq. 'dynami' .or. gmeth .eq. 'DYNAMI' )then
c       Simple dynamics
        ig_type = ig_qmd
c
      elseif( gmeth .eq. 'damped' .or. gmeth .eq. 'DAMPED' )then
c       "damped dynamics" a la Jonsson (component-wise)
        ig_type = ig_dmd
c
      elseif( gmeth .eq. 'dmd-at' .or. gmeth .eq. 'DMD-AT' )then
c       "damped dynamics" atom-wise
        ig_type = ig_dmdat
c
      elseif( gmeth .eq. 'asd   ' .or. gmeth .eq. 'ASD   ' )then
c       Accelerated Steepest Descent (PAS-modified damped dyn.)
        ig_type = ig_asd
c
      elseif( gmeth .eq. 'lbfgs ' .or. gmeth .eq. 'LBFGS ' .or.
     $        gmeth .eq. 'bfgs  ' .or. gmeth .eq. 'BFGS  ' )then
c       Limited memory BFGS scheme
        ig_type = ig_lbfgs
c       Write out warning that BFGS installation not clean:
        call FLGETIWR( IWR )
C        call BFGS_WARN( IWR )
        write(IWR,*) 'BFGS disabled - choose different minimization'
        call STOPXERR( 'BFGS disable' )
c
      elseif( gmeth .eq. 'conjug' .or. gmeth .eq. 'CONJUG' .or.
     $        gmeth .eq. 'cg    ' .or. gmeth .eq. 'CG    ' )then
c       Conjugate gradient method
        ig_type = ig_cg
        call STOPXERR( 'conjugate gradient scheme not coded yet' )
c
      elseif( gmeth(1:4) .eq. 'gtyp' .or. gmeth(1:4) .eq. 'GTYP' )then
c       *** development only *** Set igtype directly
        read(gmeth(5:6),'(i2)')  ig_type
c
      else
        call STOPXERR( 'gtypset - unknown geometry update method' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
