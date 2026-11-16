c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> YMM
c
c
      subroutine YMM( y,sth,cph,sph,m )
c---------------------------------------------------------------
c Revision history:
c  11Mar02-PAS/2.52: last literal reals to d0
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  calculates y-sub-m,m spherical harmonic
c     sth=sin(theta)
c     cph=cos(phi)
c     sph=sin(phi)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( half = 0.5d0 )
c
      DOUBLE COMPLEX    y
c
      DIMENSION  ynorm(12)
      DATA  ynorm /
     $   .28209479177388d+00,
     $   .69098829894267d+00,
     $   .15450968080928d+01,
     $   .33377905890623d+01,
     $   .70805230791198d+01,
     $   .14852230510108d+02,
     $   .30917383269124d+02,
     $   .64005064137030d+02,
     $   .13194982000571d+03,
     $   .27113110604413d+03,
     $   .55565341895088d+03,
     $   .11362831437960d+04/
c
c >>>> EXECUTABLE CODE:
c
      y = ynorm(m+1)
      if( m.ne.0 ) y = y*( -half*DCMPLX(cph,sph)*sth )**m
c
c    That's all Folks!
c
      RETURN
      END
