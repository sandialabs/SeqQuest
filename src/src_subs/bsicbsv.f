c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BSICBSV
c
c
      subroutine BSICBSV( nrs,nra, b, xx, s1,s2 )
c---------------------------------------------------------------
c Revision history:
c   8Aug00-PAS/2.45: renamed from "chbsv" to "bsicbsv"
c  22Oct99-PAS/2.40: local 2X opt/2-unroll k-loop and (nra<->2)
c  14Jun99-PAS/2.40: a simple optimization
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  ai65(),ai55() - coefs. for 6.5 and 5.5 order bessel fcns.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (n=14)
c
      DIMENSION  b(nrs,2), xx(nra),  s1(nra,2),s2(nra,2)
c local declarations:
      DIMENSION  ai55(14),ai65(14)
      DATA  (ai55(i),i=1,14) /
     2     9.84033426637529d-05,     2.99702005100095d-05,
     3     9.19213342363118d-06,     1.02971223258331d-06,
     4     1.74888326874802d-07,     1.56403099137419d-08,
     5     1.88861632972884d-09,     1.41056070635215d-10,
     6     1.33792649824125d-11,     8.59457809289520d-13,
     7     6.75553647465783d-14,     3.81127987453573d-15,
     8     2.56705246725694d-16,     1.29201110336002d-17/
      DATA  (ai65(i),i=1,14) /
     2     7.31335376759418d-06,     1.94255942832937d-06,
     3     5.83433173657527d-07,     5.89076413815039d-08,
     4     9.74426182172852d-09,     8.00554251040564d-10,
     5     9.40404192453049d-11,     6.53461776118609d-12,
     6     6.03109914162418d-13,     3.63752560527624d-14,
     7     2.78439618029305d-15,     1.48513584937351d-16,
     8     9.75140047077775d-18,     4.66539056854547d-19/
c
      DATA  zero,half / 0.0d0,0.5d0 /
c
c >>>> EXECUTABLE CODE:
c
      a1k = ai55(n)
      a2k = ai65(n)
      do  nr=1,nra
        s1(nr,1) = zero
        s1(nr,2) = a1k
        s2(nr,1) = zero
        s2(nr,2) = a2k
      enddo
c
      do 150 k=n-1,3,-2
c
        a1k = ai55(k)
        a2k = ai65(k)
        do  nr=1,nra
          s1(nr,1) = xx(nr)*s1(nr,2) - s1(nr,1) + a1k
          s2(nr,1) = xx(nr)*s2(nr,2) - s2(nr,1) + a2k
        enddo
c
        a1k = ai55(k-1)
        a2k = ai65(k-1)
        do  nr=1,nra
          s1(nr,2) = xx(nr)*s1(nr,1) - s1(nr,2) + a1k
          s2(nr,2) = xx(nr)*s2(nr,1) - s2(nr,2) + a2k
        enddo
c
  150 continue
c
      a1k = ai55(1)
      a2k = ai65(1)
      do  nr=1,nra
        b(nr,1) = half*xx(nr)*s1(nr,2) - s1(nr,1) + a1k
        b(nr,2) = half*xx(nr)*s2(nr,2) - s2(nr,1) + a2k
      enddo
c
c    That's all Folks!
c
      RETURN
      END
