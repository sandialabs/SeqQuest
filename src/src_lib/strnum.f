c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRNUM
c
c
      subroutine STRNUM( nstr, str, num, maxnum )
c---------------------------------------------------------------
c Purpose: put integer into string of appropriate length
c
c Written: Peter A. Schultz, 11-July-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Integer "num" to be written must be greater than or equal to zero.
c Integer "num" is written to character string of minimum length
c     needed to write largest possible input integer "maxnum"
c     with leading zeroes as needed.
c If error (num<0.or.maxnum<num), return null, i.e. 0-len, string
c If num=maxnum=0, return null string
c
c E.g. STRNUM( nstr, str, 13, 1000 ) returns str='0013' (nstr=4)
c      STRNUM( nstr, str,  0,   13 ) returns str='00'   (nstr=2)
c      STRNUM( nstr, str,  0,    0 ) returns nstr=0 (zero-len str)
c      STRNUM( nstr, str,113,  113 ) returns str='113'  (nstr=3)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output: character string with number
      CHARACTER*(*)  str
c
c Local declarations
      CHARACTER*(1)  chzero
      DATA           chzero / '0' /
c
c >>>> EXECUTABLE CODE:
c
      nstr = 0
c
c Check that number is within its allowed bounds:
      if( num.lt.0 .or. num.gt.maxnum .or. maxnum.eq.0 ) RETURN
c
c Determine length of full string:
      istr = 0
      mx = 1
   10 if( mx.gt.maxnum ) goto 19
      istr = istr + 1
      str(istr:istr) = chzero
      mx = mx*10
      goto 10
   19 nstr = istr
c
c Write in characters for this number:
      n0 = num
      istr = nstr
   20 nx = n0 / 10
      id = n0 - 10*nx
      write( str(istr:istr), '(i1)' ) id
      n0 = nx
      istr = istr - 1
      if( n0.gt.0 ) goto 20
c
      RETURN
c
c    That's all Folks!
c
      END
