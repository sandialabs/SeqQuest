c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> XCONMSH
c
c
      subroutine XCONMSH( natm,ntyp,nshld, nlat, n1r,n2r,n3r,nptr,
     $ itypa, numshl,lshel, ratm,hh,rlat, mngr3,ngr,
     $ nfityp,almnv1c,alxcmin,lmx1ctr,
     $ alfxc,nalfxc,cv,ce,nalfxcd,
     $ vxcmsh,excmsh,expfac,r )
c---------------------------------------------------------------
c Purpose: Initializes xc potential and energy density
c          on mesh using fit coefficients from v1ctrij
c
c Written: P.J. Feibelman
c
c Revision history:
c   1Mar00-PAS/2.43: bugfix for non-xcfit
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   1May97-PAS/2.17: naming conventions, new grid coord rtn
c  11Nov94-PAS/2.10: cleanup/first round optimization
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm), numshl(nshld),lshel(nshld,ntyp)
      DIMENSION  ratm(3,natm),hh(3,3),rlat(3,nlat)
      DIMENSION  mngr3(*),ngr(*)
      DIMENSION  nfityp(ntyp),almnv1c(ntyp),alxcmin(ntyp),lmx1ctr(ntyp)
      DIMENSION  nalfxc(ntyp),alfxc(nalfxcd,ntyp)
      DIMENSION  cv(10,nalfxcd,natm),ce(10,nalfxcd,natm)
c output arrays:
      DIMENSION  vxcmsh(nptr),excmsh(nptr)
c scratch arrays:
      DIMENSION  expfac(nptr),r(3,nptr)
c local declarations:
      DIMENSION  r0(3),rnuc(3),rlast(3),rshift(3),number(3)
      DATA  number / 1,4,10 /
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      r0(1) = zero
      r0(2) = zero
      r0(3) = zero
      call MKRBOX( hh,r0, n1r,n2r,n3r, r )
c
c    Loop on atoms:
c
      iorb = 1
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
        nrat = ngr(iorb)
        imin3 = mngr3(iorb)
        do  ishl=1,numshl(ityp)
          if( ngr(iorb).gt.nrat )then
            nrat = ngr(iorb)
            imin3 = mngr3(iorb)
          endif
          iorb = iorb + 2*lshel(ishl,ityp) + 1
        enddo
        if( nfityp(ityp).eq.0 ) goto 1000
c
        nr0 = n1r*n2r*(imin3-1)
c  Check against 1-c cutoff AND xc-fit cutoff:
        alxmin = alxcmin(ityp)
        alxmin = MAX(alxcmin(ityp),almnv1c(ityp))
        numlmx = number(lmx1ctr(ityp)+1)
c
c    Loop on fit alphas for this atom type:
c
        do  ial=1,nalfxc(ityp)
          ial1 = ial
          alxc = alfxc(ial,ityp)
          if( alxc.ge.alxmin ) goto 22
        enddo
c  No gaussians reach cutoff, skip center
        goto 1000
c
   22   continue
        rlast(1) = zero
        rlast(2) = zero
        rlast(3) = zero
c
        do 900 ilat=1,nlat
c
          do  ir=1,3
            rnuc(ir) = ratm(ir,iatm) + rlat(ir,ilat)
          enddo
          do  ir=1,3
            rshift(ir) = rnuc(ir) - rlast(ir)
          enddo
          do  ir=1,3
            rlast(ir) = rnuc(ir)
          enddo
c  Shift coordinates:
          do  nr=nr0+1,nr0+nrat
            r(1,nr) = r(1,nr) - rshift(1)
            r(2,nr) = r(2,nr) - rshift(2)
            r(3,nr) = r(3,nr) - rshift(3)
          enddo
c
c    Loop over fit alphas:
c
          do 200 ial=ial1,nalfxc(ityp)
c
c    Compute exponentials when argument is small enough.
c
            do  nr=nr0+1,nr0+nrat
              rsq = r(1,nr)**2 + r(2,nr)**2 + r(3,nr)**2
              arg = alxc*rsq
              if( arg.gt.cutgrid )then
                expfac(nr) = zero
              else
                expfac(nr) = EXP(-arg)
              endif
            enddo
c
c    Compute counter-potential using this information ...
c
            if( numlmx.ge.1 )then
c         ... s-like xc potential
              do  nr=nr0+1,nr0+nrat
                vxcmsh(nr) = vxcmsh(nr) - expfac(nr)*cv(1,ial,iatm)
                excmsh(nr) = excmsh(nr) - expfac(nr)*ce(1,ial,iatm)
              enddo
            endif
c
            if( numlmx.ge.4 )then
c         ... p-like xc potential
              do  nr=nr0+1,nr0+nrat
                x = r(1,nr)
                y = r(2,nr)
                z = r(3,nr)
                vxcmsh(nr) = vxcmsh(nr) - expfac(nr)*( cv(1,ial,iatm)
     $           + x*cv(2,ial,iatm)+y*cv(3,ial,iatm)+z*cv(4,ial,iatm) )
                excmsh(nr) = excmsh(nr) - expfac(nr)*( ce(1,ial,iatm)
     $           + x*ce(2,ial,iatm)+y*ce(3,ial,iatm)+z*ce(4,ial,iatm) )
              enddo
            endif
c
            if( numlmx.ge.10 )then
c         ... d-like xc potential
              do  nr=nr0+1,nr0+nrat
                x = r(1,nr)
                y = r(2,nr)
                z = r(3,nr)
                vxcmsh(nr) = vxcmsh(nr) - expfac(nr)*( cv(1,ial,iatm)
     $           +x*(cv(2,ial,iatm)+y*cv(5,ial,iatm)+x*cv(8,ial,iatm) )
     $           +y*(cv(3,ial,iatm)+z*cv(6,ial,iatm)+y*cv(9,ial,iatm) )
     $           +z*(cv(4,ial,iatm)+x*cv(7,ial,iatm)+z*cv(10,ial,iatm)))
                excmsh(nr) = excmsh(nr) - expfac(nr)*( ce(1,ial,iatm)
     $           +x*(ce(2,ial,iatm)+y*ce(5,ial,iatm)+x*ce(8,ial,iatm) )
     $           +y*(ce(3,ial,iatm)+z*ce(6,ial,iatm)+y*ce(9,ial,iatm) )
     $           +z*(ce(4,ial,iatm)+x*ce(7,ial,iatm)+z*ce(10,ial,iatm)))
              enddo
            endif
c
  200     continue
c
  900   continue
c  Restore coordinates:
          do  nr=nr0+1,nr0+nrat
            r(1,nr) = r(1,nr) - rlast(1)
            r(2,nr) = r(2,nr) - rlast(2)
            r(3,nr) = r(3,nr) - rlast(3)
          enddo
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
