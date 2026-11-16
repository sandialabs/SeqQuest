c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WF1COFR
c
c
      subroutine WF1COFR( orb,grorb,grgrorb,
     $ nuc,jatm, lj,jtyp,jshl,nalj,
     $ natm,ntyp,nshld,nald, nlat, ala,cala, ratm,rlat,
     $ g,gg,rji,f,gauss )
c---------------------------------------------------------------
c Purpose: for each basis function, compute value, gradient, and
c          second gradient at nucleus of all atoms but its own.
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c  26Aug94-PAS/2.11: cleaned up, underflow protection
c---------------------------------------------------------------
c
c Output used to construct charge density as sum of 1-ctr stuff
c containing all short-range near-nucleus rapid variation
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      DIMENSION  orb(nlat),grorb(3,nlat),grgrorb(3,3,nlat)
      DIMENSION  ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  ratm(3,natm),rlat(3,nlat)
c scratch declarations:
      DIMENSION  g(3,nlat),gg(3,3,nlat),rji(3,nlat),
     $           f(nlat),gauss(nlat)
c
c local declarations:
      DIMENSION  vnorm(9)
      DATA  zero,half,one,two, four / 0.d0,0.5d0,1.d0,2.d0, 4.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutexp
      do  jlat=1,nlat
        do  jr=1,3
          do  ir=1,3
            grgrorb(ir,jr,jlat) = zero
            gg(ir,jr,jlat) = zero
          enddo
          grorb(jr,jlat) = zero
          g(jr,jlat) = zero
        enddo
        orb(jlat) = zero
      enddo
c
c      lattice vector loop
c
      do  jr=1,3
        do  jlat=1,nlat
          rji(jr,jlat) = ratm(jr,nuc) - ratm(jr,jatm) - rlat(jr,jlat)
        enddo
      enddo
c
      goto (10,20,30,40,50,60,70,80,90) lj
c
   10 do  jlat=1,nlat
        f(jlat) = one
      enddo
      goto 200
c
   20 do  jlat=1,nlat
        f(jlat) = rji(1,jlat)
        g(1,jlat) = one
      enddo
      goto 200
c
   30 do  jlat=1,nlat
        f(jlat) = rji(2,jlat)
        g(2,jlat) = one
      enddo
      goto 200
c
   40 do  jlat=1,nlat
        f(jlat) = rji(3,jlat)
        g(3,jlat) = one
      enddo
      goto 200
c
   50 do  jlat=1,nlat
        f(jlat) = rji(1,jlat)*rji(2,jlat)
        g(1,jlat) = rji(2,jlat)
        g(2,jlat) = rji(1,jlat)
        gg(1,2,jlat) = one
        gg(2,1,jlat) = one
      enddo
      goto 200
c
   60 do  jlat=1,nlat
        f(jlat) = rji(2,jlat)*rji(3,jlat)
        g(2,jlat) = rji(3,jlat)
        g(3,jlat) = rji(2,jlat)
        gg(2,3,jlat) = one
        gg(3,2,jlat) = one
      enddo
      goto 200
c
   70 do  jlat=1,nlat
        f(jlat) = rji(3,jlat)*rji(1,jlat)
        g(1,jlat) = rji(3,jlat)
        g(3,jlat) = rji(1,jlat)
        gg(3,1,jlat) = one
        gg(1,3,jlat) = one
      enddo
      goto 200
c
   80 do  jlat=1,nlat
        f(jlat) = rji(1,jlat)**2 - rji(2,jlat)**2
        g(1,jlat) =  two*rji(1,jlat)
        g(2,jlat) = -two*rji(2,jlat)
        gg(1,1,jlat) =  two
        gg(2,2,jlat) = -two
      enddo
      goto 200
c
   90 do  jlat=1,nlat
        f(jlat) = two*rji(3,jlat)**2 - rji(1,jlat)**2 - rji(2,jlat)**2
        g(1,jlat) = -two*rji(1,jlat)
        g(2,jlat) = -two*rji(2,jlat)
        g(3,jlat) = four*rji(3,jlat)
        gg(1,1,jlat) = -two
        gg(2,2,jlat) = -two
        gg(3,3,jlat) = four
      enddo
      goto 200
c
  200 continue
c
c    Alpha loop:
c
      vn = vnorm(lj)
      do 400 jal=1,nalj
        alj = ala(jal,jshl,jtyp)
        calj = cala(jal,jshl,jtyp)/vn
c
c  evaluate value of orbital j at site nuc
        rjimx = argmx/alj
        do  jlat=1,nlat
          rjisq = rji(1,jlat)**2 + rji(2,jlat)**2 + rji(3,jlat)**2
          if( rjisq .gt. rjimx )then
            gauss(jlat) = zero
          else
            gauss(jlat) = calj*EXP( -alj*rjisq )
            orb(jlat) = orb(jlat) + gauss(jlat)*f(jlat)
          endif
        enddo
c
c  evaluate gradients of orbital j at site nuc
        do 321 n=1,3
          do  jlat=1,nlat
            grorb(n,jlat) = grorb(n,jlat) + gauss(jlat)*
     $       ( g(n,jlat) - two*alj*rji(n,jlat)*f(jlat) )
            grgrorb(n,n,jlat) = grgrorb(n,n,jlat) - gauss(jlat)*
     $       two*alj*f(jlat)
          enddo
  321   continue
c
        do 332 n=1,3
          do 331 nn=1,3
            do  jlat=1,nlat
              grgrorb(nn,n,jlat) = grgrorb(nn,n,jlat) + gauss(jlat)*
     $         ( gg(nn,n,jlat)
     $         - two*alj*(rji(n,jlat)*g(nn,jlat)+rji(nn,jlat)*g(n,jlat))
     $         + four*alj**2*rji(n,jlat)*rji(nn,jlat)*f(jlat) )
            enddo
  331     continue
  332   continue
c
  400 continue
c
c    That's all Folks!
c
      RETURN
      END
