c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKCUTS
c
c
      subroutine MKCUTS( convgf, convsl, convii )
c---------------------------------------------------------------
c Purpose: construct cutoffs from convergence parameters
c
c Written: Peter A. Schultz, 13-August-1998, for v2.24
c
c Revision history:
c  21Jun01-PAS/2.48: replaces STOPs
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  I take the input convergence parameters, translate them into
c  exponent cutoffs, and communicate those cutoffs via common
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      cutgwfn = -LOG(convgf)
      cutslow = cutgwfn - LOG(convsl)
      cutii   = cutgwfn - LOG(convii)
      if( cutgwfn.lt.zero .or. cutslow.lt.zero .or.
     $ cutii.lt.zero ) call STOPXERR( 'bad-cuts/bad conv params' )
c
c    That's all, Folks!
c
      RETURN
      END
