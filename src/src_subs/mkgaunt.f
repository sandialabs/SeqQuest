c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKGAUNT
c
c
      subroutine MKGAUNT( gntcomb,gaunt )
c---------------------------------------------------------------
c Purpose: compute linear combinations of gaunt coefficients
c          needed because of the use of real spherical harmonics
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   8Aug94-PAS/2.11: routine still a mess even after makeover
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
      DIMENSION  gaunt(9,25,3)
c local declarations:
      DOUBLE COMPLEX    zone,fre,fim,fac1,fc1,fac2,fc2
      DIMENSION  l(9),m(9),irl(9)
      DATA  l   / 0, 1,1,1, 2,2,2,2,2 /
      DATA  m   / 0, 1,1,0, 2,1,1,2,0 /
      DATA  irl / 1, 1,2,1, 2,2,1,1,1 /
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
      LOGICAL  ODD
c statement function:
      ODD(int) = (AND(int,1).eq.1)
c
c create gaunt coefficients:
c
      do  ll1=1,3
        l1 = ll1 - 1
        do  mm1=1,2*l1+1
          m1 = mm1 - 1 - l1
          lm1 = l1**2 + mm1
          do  ll2=1,5
            l2 = ll2 - 1
            do  mm2=1,2*l2+1
              m2 = mm2 - 1 - l2
              lm2 = l2**2+mm2
              do 10 kk=1,3
                k = kk - 1
                if( IABS(l1-l2).gt.k .or. (l1+l2).lt.k .or.
     $              IABS(m1-m2).gt.k ) goto 10
                gaunt(lm1,lm2,kk) = CEE( l1,m1,l2,m2,k )
   10         continue
            enddo
          enddo
        enddo
      enddo
c
c clear array of stored values:
c
      do  lmp1=1,9
       do  ll2=1,5
        do  mm2=1,5
         do  i2=1,2
          do  km=1,9
           gntcomb(lmp1,ll2,mm2,i2,km) = (0.,0.)
          enddo
         enddo
        enddo
       enddo
      enddo
c
      do 500 lmp1=1,9
      l1 = l(lmp1)
      m1 = m(lmp1)
      i1 = irl(lmp1)
c
      do 430 mm2=1,5
      m2 = mm2 - 1
      do 420 ll2=mm2,5
      l2 = ll2 - 1
      do 410 i2=1,2
c
      do 300 k1=1,3
        k = k1 - 1
c
c  triangle inequality
        if( k.gt.(l1+l2) .or. k.lt.IABS(l1-l2) ) goto 300
c
        zone = DCMPLX( one, zero )
        rt2 = one/SQRT(two)
        fre = DCMPLX( rt2, zero )
        fim = DCMPLX( zero, rt2 )
c
        do 200 im1=1,2
          if( m1.eq.0 )then
            if( im1.eq.2 ) goto 200
            fac1 = zone
            mc1 = 0
          else
            if( i1.eq.2 )then
              fc1 = fim
            else
              fc1 = fre
            endif
            if( im1.eq.2 )then
              fac1 = fc1
              if( ODD(i1-1) ) fac1 = -fac1
              mc1 = -m1
            else
              fac1 = fc1
              if( ODD(m1) ) fac1 = -fac1
              mc1 = m1
            endif
          endif
c
          do 100 im2=1,2
            if( m2.eq.0 )then
              if( im2.eq.2 ) goto 100
              fac2 = zone
              mc2 = 0
            else
              if( i2.eq.2 )then
                fc2 = fim
              else
                fc2 = fre
              endif
              if( im2.eq.2 )then
                fac2 = fc2
                if( ODD(m2+i2-1) ) fac2 = - fac2
                mc2 = -m2
              else
                fac2 = fc2
                mc2 = m2
              endif
            endif
            msubk = mc1 + mc2
            if( IABS(msubk).gt.k ) goto 100
c
c      compute m-sub-k index
c
            km = k**2 + k + msubk + 1
            gntcomb(lmp1,ll2,mm2,i2,km) = gntcomb(lmp1,ll2,mm2,i2,km) +
     $       fac1*fac2* gaunt( l1**2+l1-mc1+1, l2**2+l2+mc2+1, k+1 )
c
  100     continue
  200   continue
  300 continue
c
  410 continue
  420 continue
  430 continue
  500 continue
c
      RETURN
      END
