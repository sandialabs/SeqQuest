c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FTORB
c
c
      subroutine FTORB( orbint, norb,nk, natm,ntyp,nshld,nald,
     $ itypa, numshl,lshel,nala,ala,cala, norba, vk )
c---------------------------------------------------------------
c Purpose: compute exact fourier transform of each orbital
c
c Written: P.J. Feibelman, revised by P.A. Schultz
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c   7Mar02-PAS/2.52: ang constants extracted
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c output array:
      DIMENSION  orbint(2,norb,nk)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  norba(ntyp)
      DIMENSION  vk(3,nk)
c
c local declarations:
      DIMENSION  oint(5)
c ph() contains i**l
      DIMENSION  ph(2,0:3)
      DATA  ph / 1.d0,0.d0, 0.d0,1.d0, -1.d0,0.d0, 0.d0,-1.d0 /
      DIMENSION  vnorm(5)
c
      DATA  zero,one,two, four / 0.d0,1.d0,2.d0, 4.d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
      do 200 k=1,nk
        xk = vk(1,k)
        yk = vk(2,k)
        zk = vk(3,k)
        vksq = xk*xk + yk*yk + zk*zk
        vksqd4 = vksq / four
c
        do 100 ityp=1,ntyp
c
          iorbat0 = 0
          do 80 ishl=1,numshl(ityp)
            lshl = lshel(ishl,ityp)
            norbsh = 2*lshl + 1
            nali = nala(ishl,ityp)
c
            sumg = zero
            do 30 ial=1,nali
              ali = ala(ial,ishl,ityp)
c 07dec93-PAS: following replaces dangerous construction
c      cali = (pi/ali)**1.5*cala(ial,ishl,ityp)* EXP(-0.25*vksq/ali)
              alfac = (pi/ali)*SQRT( pi/ali )
              addg = alfac*cala(ial,ishl,ityp)*EXP( -vksqd4/ali )
              if( lshl.gt.0 )then
                do  l=1,lshl
                  addg = addg/(ali+ali)
                enddo
              endif
              sumg = sumg + addg
   30       continue
c
c install angular factor
c
            goto (40,41,42), lshl+1
c     s-fcn
   40       fac = sumg/vnorm(1)
            oint(1) = fac
            goto 50
c     p-fcn
   41       fac = sumg/vnorm(2)
            oint(1) = -xk*fac
            oint(2) = -yk*fac
            oint(3) = -zk*fac
            goto 50
c     d-fcn
   42       fac = sumg/vnorm(3)
            oint(1) = xk*yk*fac
            oint(2) = yk*zk*fac
            oint(3) = zk*xk*fac
            fac = sumg/vnorm(4)
            oint(4) = (xk*xk-yk*yk)*fac
            fac = sumg/vnorm(5)
            oint(5) = (two*zk*zk-xk*xk-yk*yk)*fac
c
   50       continue
c
c search for atoms of this type, multiply by i**l, deposit values:
c
            i0 = 0
            do 60 iatm=1,natm
              if( itypa(iatm).eq.ityp )then
                i = i0 + iorbat0
                do  ish=1,norbsh
                  i = i + 1
                  orbint(1,i,k) = orbint(1,i,k) + oint(ish)*ph(1,lshl)
                  orbint(2,i,k) = orbint(2,i,k) + oint(ish)*ph(2,lshl)
                enddo
              endif
              i0 = i0 + norba(itypa(iatm))
   60       continue
            iorbat0 = iorbat0 + norbsh
c
   80     continue
c
  100   continue
c
  200 continue
c
      if( lstout.ge.1 )then
        call FLGETIWR( IWR )
        mk = nk
        if( lstout.eq.1 ) mk = 1
        kn = 0
        write(IWR,*) 'FTORB: Int(orb)='
  300   continue
        if( kn.lt.mk )then
          k1 = kn + 1
          kn = k1 + 2
          if( kn.gt.mk ) kn = mk
          write(IWR,9300) k1,k1+1,k1+2
 9300     format(1x,'i',2x, 3(1x, 9x,'#k=',i2,10x) )
          do  i=1,norb
            write(IWR,9310) i,(orbint(1,i,k),orbint(2,i,k),k=k1,kn)
 9310       format(1x,i3,   3(1x, 1p2d12.4) )
          enddo
          goto 300
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
