c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONECTR
c
c
      subroutine ONECTR( vij1,vradsr1,vij2,vradsr2,vij3,vradsr3,
     $ ipoint,na,lmx1ctr, alsum,p,li,lj,ci,cj )
c---------------------------------------------------------------
c Purpose: generalization of engaus for use with "vofn"
c
c Written: P.J. Feibelman
c
c Revision history:
c  24Sep08-PAS/2.62: change to get rid of compiler warning
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  15May89-PAS: pack vradsr, optimize code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vij1(na),vij2(na),vij3(na)
c input arrays:
      DIMENSION  vradsr1(70,na),vradsr2(70,na),vradsr3(70,na)
      DIMENSION  p(3),ci(10),cj(10)
c local declarations:
      DIMENSION  potcoef(5,5,3,3),npot(3),ntimes(3),bfac(70)
      EQUIVALENCE  (nq(1),nq1),(nqmn(1),nqmn1),(psq(1),psq1)
      COMMON  /nqs/   nq1,nq2,nq3
      COMMON  /nqmns/ nqmn1,nqmn2,nqmn3
      COMMON  /psqs/  psq1,psq2,psq3
      DIMENSION  nq(3),nqmn(3),psq(3)
      DIMENSION  nsx(2)
      EQUIVALENCE  (nsx(1),nsi),(nsx(2),nsj)
      LOGICAL  EVEN
c
      DIMENSION  number(3),lang(10),istart(7)
      DATA  number / 1,4,10 /
      DATA  lang / 0, 1,1,1, 2,2,2,2,2,2 /
      DATA  istart / 1,9,10,18,19,27,28 /
c
      EQUIVALENCE  (bpq,beq)
      DIMENSION  beq(49)
      DIMENSION  ns(9),np(10,3),bpq(7,7),nqmin(7)
      DATA  ns / 1, 2,3,4, 5,6,7,9,10 /
      DATA  np / 0, 1,0,0, 1,0,1,2,0,0,
     $           0, 0,1,0, 1,1,0,0,2,0,
     $           0, 0,0,1, 0,1,1,0,0,2 /
      DATA  bpq / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $            0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $            0.d0, 0.d0, 1.d0, 3.d0, 3.d0, 0.d0, 0.d0,
     $            0.d0, 0.d0, 0.d0, 1.d0, 6.d0,15.d0,15.d0,
     $            0.d0, 0.d0, 0.d0, 0.d0, 1.d0,10.d0,45.d0,
     $            0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0,15.d0,
     $            0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
      DATA  nqmin / 1,2,2,3,3,4,4 /
c
      DATA  zero,one / 0.d0,1.d0 /
c statement function:
cvax  EVEN(int) = ((int.and.1).eq.0)
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
c  Set potential coefficients ...
c   ... zeroth derivative of potential:
      do  j=1,3
        do  npower=1,5
          potcoef(npower,1,1,j) = one
        enddo
      enddo
      if( lmx1ctr.gt.0 )then
c  ... first derivative of potential:
        alsumt2 = alsum + alsum
        do  j=1,3
          do  npow=1,5
            potcoef(npow,1,2,j) =  DBLE(npow-1)
            potcoef(npow,2,2,j) = p(j)
            potcoef(npow,3,2,j) = -alsumt2
          enddo
        enddo
        if( lmx1ctr.gt.1 )then
c  ... second derivative of potential:
          alsumt4 = alsumt2 + alsumt2
          do  j=1,3
            pj = p(j)
            do  npow=1,5
              potcoef(npow,1,3,j) =  DBLE( MAX( (npow-1)*(npow-2), 0 ) )
              potcoef(npow,2,3,j) = (pj+pj)* DBLE(npow-1)
              potcoef(npow,3,3,j) = pj*pj -  DBLE(npow+npow-1)*alsumt2
              potcoef(npow,4,3,j) = -alsumt4*pj
              potcoef(npow,5,3,j) = alsumt4*alsum
            enddo
          enddo
        endif
      endif
c
      nv = lang(li) + lang(lj) + lmx1ctr + 1
      nqv = nv*number(lmx1ctr+1)
      do  iqv=1,nqv
        bfac(iqv) = zero
      enddo
      p1 = p(1)
      p2 = p(2)
      p3 = p(3)
      psq1 = p1*p1
      psq2 = p2*p2
      psq3 = p3*p3
c
      if( li.le.lj )then
        ki = 1
        kj = 2
        nsimx = ns(li)
        nsjmx = ns(lj)
      else
        ki = 2
        kj = 1
        nsimx = ns(lj)
        nsjmx = ns(li)
      endif
c
      do 102 nsi=1,nsimx
        do 101 nsj=nsi,nsjmx
          if( nsj.ne.nsi .and. nsj.le.nsimx )then
            pr0 = ci(nsi)*cj(nsj) + ci(nsj)*cj(nsi)
          else
            pr0 = ci(nsx(ki))*cj(nsx(kj))
          endif
          if( pr0.eq.zero ) goto 101
c
          npow1 = np(nsi,1) + np(nsj,1) + 1
          npow2 = np(nsi,2) + np(nsj,2) + 1
          npow3 = np(nsi,3) + np(nsj,3) + 1
          iv = -2
          do 90 lpot=1,number(lmx1ctr+1)
            do  j=1,3
              npj = np(lpot,j)
              npot(j) = npj
              ntimes(j) = npj + npj + 1
            enddo
c
            do 80 nt1=1,ntimes(1)
              nq1 = npow1 - 1 - npot(1) + nt1
              if( nq1.lt.1 ) goto 80
              potc1 = potcoef(npow1,nt1,np(lpot,1)+1,1)
              if( potc1.eq.zero ) goto 80
              nqmn1 = nqmin(nq1)
              ist1 = istart(nq1)
              if( EVEN(nq1) ) potc1 = p1*potc1
c
              do 70 nt2=1,ntimes(2)
                nq2 = npow2 - 1 - npot(2) + nt2
                if( nq2.lt.1 ) goto 70
                potc12 = potcoef(npow2,nt2,np(lpot,2)+1,2)
                if( potc12.eq.zero ) goto 70
                potc12 = potc12*potc1
                nqmn2 = nqmin(nq2)
                ist2 = istart(nq2)
                if( EVEN(nq2) ) potc12 = p2*potc12
c
                do 60 nt3=1,ntimes(3)
                  nq3 = npow3 - 1 - npot(3) + nt3
                  if( nq3.lt.1 ) goto 60
                  potc123 = potcoef(npow3,nt3,np(lpot,3)+1,3)
                  if( potc123.eq.zero ) goto 60
                  potc123 = potc123*potc12
                  nqmn3 = nqmin(nq3)
                  ist3 = istart(nq3)
                  if( EVEN(nq3) ) potc123 = p3*potc123
c
                  pr1 = pr0*potc123
                  ix1 = ist1
c     cdir$ novector
                  do 50 iq1=nqmn1,nq1
                    pr12 = pr1
                    ix2 = ist2
                    do 40 iq2=nqmn2,nq2
                      pr123 = pr12
                      do 30 iq3=nqmn3,nq3
                        bfac(iv+iq1+iq2+iq3) = bfac(iv+iq1+iq2+iq3) +
     $                   pr123*beq(ix1)*beq(ix2)*bpq(nq3,iq3)
                        pr123 = pr123*psq3
   30                 continue
                      pr12 = pr12*psq2
                      ix2 = ix2 + 7
   40               continue
                    pr1 = pr1*psq1
                    ix1 = ix1 + 7
   50             continue
c     cdir$ vector
   60           continue
   70         continue
   80       continue
c
            iv = iv + lang(li) + lang(lj) + lang(lpot) + 1
   90     continue
c
  101   continue
  102 continue
c
c      nqv = iv + 2
      do  ia=1,ipoint+1
        sum1 = zero
        sum2 = zero
        sum3 = zero
        do  iqv=1,nqv
          sum1 = sum1 + bfac(iqv)*vradsr1(iqv,ia)
          sum2 = sum2 + bfac(iqv)*vradsr2(iqv,ia)
          sum3 = sum3 + bfac(iqv)*vradsr3(iqv,ia)
        enddo
        vij1(ia) = sum1
        vij2(ia) = sum2
        vij3(ia) = sum3
      enddo
c
c    That's all Folks!
c
      RETURN
      END
