c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLANG
c
c
      subroutine VLANG( vij1,vint1,vij2,vint2,vij3,vint3,vij4,vint4,
     $ p,li,lj,ci,cj )
c---------------------------------------------------------------
c Purpose: compute quintuple sum of symmetry coefficients times
c          the radial integrals in "venint"
c
c Revision history:
c   6Mar00-PAS/2.43: rename "engaus" to "vlang"
c  18Nov99-PAS/2.42: cosmetic cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  21Jun89-PAS: (1) EVEN fcn (2) fold over nsi,nsj sums
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  vint1(*),vint2(*),vint3(*),vint4(*)
      DIMENSION  p(3),ci(10),cj(10)
c local declarations:
      EQUIVALENCE  (nq(1),nq1),(nqmn(1),nqmn1),(psq(1),psq1)
      COMMON  /nqs/   nq1,nq2,nq3
      COMMON  /nqmns/ nqmn1,nqmn2,nqmn3
      COMMON  /psqs/  psq1,psq2,psq3
      DIMENSION  nq(3),nqmn(3),psq(3)
      DIMENSION  nsx(2)
      EQUIVALENCE  (nsx(1),nsi), (nsx(2),nsj)
      LOGICAL    EVEN
c
      DIMENSION  t(5)
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
      DATA  zero / 0.d0 /
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      nsimx = ns(li)
      nsjmx = ns(lj)
      if( nsimx.le.nsjmx )then
        ki = 1
        kj = 2
      else
        ki = 2
        kj = 1
        itmp = nsimx
        nsimx = nsjmx
        nsjmx = itmp
      endif
c
      do  j=1,3
        psq(j) = p(j)*p(j)
      enddo
c
      do  i=1,5
        t(i) = zero
      enddo
c
      vij1 = zero
      vij2 = zero
      vij3 = zero
      vij4 = zero
      do 102 nsi=1,nsimx
        do 101 nsj=nsi,nsjmx
c
          if( nsj.ne.nsi .and. nsj.le.nsimx )then
            pr1 = ci(nsi)*cj(nsj) + ci(nsj)*cj(nsi)
          else
            pr1 = ci(nsx(ki))*cj(nsx(kj))
          endif
          if( pr1.eq.zero ) goto 101
c
          do  j=1,3
            nq(j) = np(nsi,j) + np(nsj,j) + 1
            nqmn(j) = nqmin(nq(j))
            if( EVEN(nq(j)) ) pr1 = pr1*p(j)
          enddo
c
c     cdir$ novector
          do 50 iq1=nqmn1,nq1
            pr12 = pr1*bpq(nq1,iq1)
            do 40 iq2=nqmn2,nq2
              pr123 = pr12*bpq(nq2,iq2)
              do 30 iq3=nqmn3,nq3
                t(iq1+iq2+iq3-2) = t(iq1+iq2+iq3-2) + pr123*bpq(nq3,iq3)
                pr123 = pr123*psq3
   30         continue
              pr12 = pr12*psq2
   40       continue
            pr1 = pr1*psq1
   50     continue
c     cdir$ vector
c
  101   continue
  102 continue
c
      do  i=1,5
        vij1 = vij1 + t(i)*vint1(i)
        vij2 = vij2 + t(i)*vint2(i)
        vij3 = vij3 + t(i)*vint3(i)
        vij4 = vij4 + t(i)*vint4(i)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
