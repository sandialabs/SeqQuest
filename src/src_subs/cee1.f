c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CEE1
c
c
      FUNCTION CEE1( a,b,c, d,e,f )
c---------------------------------------------------------------
c Revision history:
c  24Sep08-PAS/2.62: cosmetic to removed compiler warning
c   7Mar02-PAS/2.52: explicitly promote constants to dp (d0)
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      INTEGER  a,b,c, d,e,f
c local declarations:
      DIMENSION  fact(21)
      DATA  fact /
     $  1.d0,1.d0,2.d0,6.d0,24.d0,120.d0,720.d0,5040.d0,40320.d0,
     $  362880.d0,3628800.d0,39916800.d0,479001600.d0,6227020800.d0,
     $  87178291200.d0,1307674368000.d0,20922789888000.d0,
     $  3.55687428096d+14,6.402373705728d+15,
     $  1.216451004088d+17,2.432902008177d+18 /
      DATA  zero,one / 0.d0,1.d0 /
      LOGICAL  ODD
c statement function:
      ODD(int) = (AND(int,1).eq.1)
c
c >>>> EXECUTABLE CODE:
c
      p = fact(a+b-c+1)
      q = fact(a-b+c+1)
      r = fact(b+c-a+1)
      s = fact(a+b+c+2)
      t = fact(a+d+1)
      u = fact(a-d+1)
      v = fact(b+e+1)
      w = fact(b-e+1)
      x = fact(c+f+1)
      y = fact(c-f+1)
      p = SQRT( p*q*r/s ) * SQRT( t*u*v*w*x*y )
      if( ODD(c-f) ) p = -p
      ceesum = zero
      is = -1
      q = -one
c
  100 continue
      is = is + 1
      q = -q
      r = fact(is+1)
      i = a + b - c - is
      if( i.lt.0 ) goto 990
      s = fact(i+1)
      i = a - d - is
      if( i.lt.0 ) goto 990
      t = fact(i+1)
      i = b + e - is
      if( i.lt.0 ) goto 990
      u = fact(i+1)
      i = -b + c + d + is
      if( i.lt.0 ) goto 100
      v = fact(i+1)
      i = -a + c - e + is
      if( i.lt.0 ) goto 100
      w = fact(i+1)
      ceesum = ceesum + p*q/(r*s*t*u*v*w)
      cee1 = ceesum
      goto 100
c
  990 continue
      cee1 = ceesum
c
c    That's all Folks!
c
      RETURN
      END
