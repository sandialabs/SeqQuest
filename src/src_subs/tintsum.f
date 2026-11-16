c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TINTSUM
c
c
      subroutine TINTSUM( telem, p,li,lj,ci,cj,ali,alj,alsum,pmag )
c---------------------------------------------------------------
c Purpose: compute quintuple sum of symmetry coefficients
c          times the radial integrals for the kinetic matrix
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      DIMENSION  p(3), ci(10),cj(10)
c local declarations:
      EQUIVALENCE  (nq(1),nq1),(nqmn(1),nqmn1),(psq(1),psq1)
      COMMON  /nqs/   nq1,nq2,nq3
      COMMON  /nqmns/ nqmn1,nqmn2,nqmn3
      COMMON  /psqs/  psq1,psq2,psq3
      DIMENSION  nq(3),nqmn(3),psq(3)
c
      LOGICAL    EVEN
      DIMENSION  trad(5,3),lang(10)
      DATA  lang / 0, 1,1,1, 2,2,2,2,2,2 /
c
      DIMENSION  ns(9),np(10,3),bpq(5,5),nqmin(5)
      DATA  ns / 1, 2,3,4, 5,6,7,9,10 /
      DATA  np / 0, 1,0,0, 1,0,1,2,0,0,
     $           0, 0,1,0, 1,1,0,0,2,0,
     $           0, 0,0,1, 0,1,1,0,0,2 /
      DATA  bpq / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $            0.d0, 1.d0, 1.d0, 0.d0, 0.d0,
     $            0.d0, 0.d0, 1.d0, 3.d0, 3.d0,
     $            0.d0, 0.d0, 0.d0, 1.d0, 6.d0,
     $            0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
      DATA  nqmin / 1,2,2,3,3 /
c
      DATA  zero,one,two,three / 0.d0,1.d0,2.d0,3.d0 /
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      trc = two*rtpi*alj*alj/(alsum*SQRT(alsum))
      cl = (three*ali/alj + two - pmag*pmag/(two*alsum))/two
C      do  lp1=1,3
      do  lp1=1,lang(lj)+1
        cliq = cl
        fiq = trc
C        do  iqs=1,lang(li)+lang(lj)+1
        do  iqs=1,5
          fiq = fiq/(alsum+alsum)
          cliq = cliq - one
          trad(iqs,lp1) = fiq*cliq
        enddo
        cl = cl + alsum/alj
      enddo
      do  j=1,3
        psq(j) = p(j)*p(j)
      enddo
c
      nsimx = ns(li)
      nsjmx = ns(lj)
c
      telem = zero
      nsjmn = lj
      if( lj.eq.9 ) nsjmn = 8
c
      do 102 nsi=1,nsimx
        if( ci(nsi).eq.zero ) goto 102
c
        do 101 nsj=nsjmn,nsjmx
          if( cj(nsj).eq.zero ) goto 101
c
          pr1 = ci(nsi)*cj(nsj)
          lp1 = lang(nsj) + 1
          do 20 j=1,3
            nq(j) = np(nsi,j) + np(nsj,j) + 1
            nqmn(j) = nqmin(nq(j))
            if( EVEN(nq(j)) ) pr1 = pr1*p(j)
   20     continue
c     cdir$ novector
          do 50 iq1=nqmn1,nq1
            pr12 = pr1
            do 40 iq2=nqmn2,nq2
              pr123 = pr12
              do 30 iq3=nqmn3,nq3
                telem = telem + trad(iq1+iq2+iq3-2,lp1)*
     $           pr123*bpq(nq1,iq1)*bpq(nq2,iq2)*bpq(nq3,iq3)
                pr123 = pr123*psq3
   30         continue
              pr12 = pr12*psq2
   40       continue
            pr1 = pr1*psq1
   50     continue
c     cdir$ vector
  101   continue
  102 continue
c
      RETURN
      END
