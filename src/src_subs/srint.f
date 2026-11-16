c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SRINT
c
c
      subroutine SRINT( vradsr1,v1ctr1,vradsr2,v1ctr2,vradsr3,v1ctr3,
     $ alsum,pmag,exparg, lmx1ctr,ltop,  na,nuc,neq,ieqn,
     $ nrs,nr0,nrd, radmsh,rwght,
     $ besstor,radfac,w1,w2,w3,w4,w5,w6 )
c---------------------------------------------------------------
c Written: P.J. Feibelman
c
c Revision history:
c   7Dec99-PAS/2.42: cosmetic loop cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c    Aug94-PAS/2.11: reopt, clean, rationalize
c   2Feb93-PAS: put in cutoff control
c  14May89-PAS: install optimized bsiall/opt code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vradsr1(7*10,neq+1)
      DIMENSION  vradsr2(7*10,neq+1)
      DIMENSION  vradsr3(7*10,neq+1)
c input arrays:
      DIMENSION  v1ctr1(nrd,10,na)
      DIMENSION  v1ctr2(nrd,10,na)
      DIMENSION  v1ctr3(nrd,10,na)
      DIMENSION  radmsh(nrs),rwght(nrs), ieqn(*)
c scratch arrays:
      DIMENSION  besstor(nrs,7),radfac(nrs,7),
     $ w1(nrs),w2(nrs),w3(nrs),w4(nrs),w5(nrs),w6(nrs)
c local declarations:
      DIMENSION  number(3),lang(10)
      DATA  number / 1,4,10 /
      DATA  lang / 0, 1,1,1, 2,2,2,2,2,2 /
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      llmax = ltop + lmx1ctr
c
      if( pmag.gt.zero )then
        do  nr=1,nrs
          r = radmsh(nr)
          totarg = (-alsum*r+pmag)*r - exparg
          radfac(nr,1) = rwght(nr)*EXP(totarg)
        enddo
      else
        do  nr=1,nrs
          radfac(nr,1) = rwght(nr)*EXP(-alsum*radmsh(nr)**2-exparg)
        enddo
      endif
c
c    Compute bessel functions on radial mesh
c
      if( pmag.gt.zero )then
c
        call BSIALL( nrs, llmax, pmag, besstor, radmsh,
     $   w1,w2,w3,w4,w5,w6 )
c
        if( llmax.gt.1 )then
          denom = one/pmag
          do 40 ll=2,llmax
            do  nr=1,nrs
              radfac(nr,ll) = radfac(nr,ll-1)*radmsh(nr)*denom
              radfac(nr,ll-1) = radfac(nr,ll-1)*besstor(nr,ll-1)
            enddo
   40     continue
        endif
        do  nr=1,nrs
          radfac(nr,llmax) = radfac(nr,llmax)*besstor(nr,llmax)
        enddo
c
      else
c
        if( llmax.gt.1 )then
          denom = one
          do 140 ll=2,llmax
            denom = one/ DBLE(2*ll-1)
            do  nr=1,nrs
              radfac(nr,ll) = radfac(nr,ll-1)*radmsh(nr)**2*denom
            enddo
  140     continue
        endif
c
      endif
c
c    Compute integrals
c
      llp0 = 0
      do 240 lpot=1,number(lmx1ctr+1)
c
        llp = llp0
        do 180 ll=1,ltop+lang(lpot)
          llp = llp + 1
          sum1 = zero
          sum2 = zero
          sum3 = zero
          do  nr=1,nrs
            sum1 = sum1 + v1ctr1(nr0+nr,lpot,nuc)*radfac(nr,ll)
            sum2 = sum2 + v1ctr2(nr0+nr,lpot,nuc)*radfac(nr,ll)
            sum3 = sum3 + v1ctr3(nr0+nr,lpot,nuc)*radfac(nr,ll)
          enddo
          vradsr1(llp,1) = vradsr1(llp,1) + sum1
          vradsr2(llp,1) = vradsr2(llp,1) + sum2
          vradsr3(llp,1) = vradsr3(llp,1) + sum3
  180   continue
c
c    Obtain vradsr for equivalent atom triples, if any
c
        if( neq.gt.0 )then
c
          do 190 ieq=1,neq
            jeq = ieqn(ieq)
c
            llp = llp0
            do 185 ll=1,ltop+lang(lpot)
              llp = llp + 1
              sum1 = zero
              sum2 = zero
              sum3 = zero
              do  nr=1,nrs
                sum1 = sum1 + v1ctr1(nr0+nr,lpot,jeq)*radfac(nr,ll)
                sum2 = sum2 + v1ctr2(nr0+nr,lpot,jeq)*radfac(nr,ll)
                sum3 = sum3 + v1ctr3(nr0+nr,lpot,jeq)*radfac(nr,ll)
              enddo
              vradsr1(llp,ieq+1) = vradsr1(llp,ieq+1) + sum1
              vradsr2(llp,ieq+1) = vradsr2(llp,ieq+1) + sum2
              vradsr3(llp,ieq+1) = vradsr3(llp,ieq+1) + sum3
  185       continue
c
  190     continue
c
        endif
c
        llp0 = llp
  240 continue
c
      RETURN
      END
