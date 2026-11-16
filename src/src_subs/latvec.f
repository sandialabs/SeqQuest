c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LATVEC
c
c
      subroutine LATVEC( ndim, nlatd,nlat,nlat1c,nlat2c,
     $ natm,ntyp,nshld,nald, itypa, numshl,nala,ala,
     $ ratm,rprim, rlat, rl1c,rl2c )
c---------------------------------------------------------------
c Purpose: computes lattice vectors so hamiltonian is converged
c     to degree specified by "cutgwfn", given longest gaussians
c
c Written: Peter A. Schultz, 8-November-1994, based on PJF's
c  code, which in turn was based on appelbaum/hamann lcao code
c
c Revision history:
c  14Jun05-PAS/2.59: File units, with a fix
c  21Jun01-PAS/2.48: replace STOPs
c  27Mar00-PAS/2.43: separate out calculations of Bloch phases
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  24Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp)
      DIMENSION  ratm(3,natm),rprim(3,3)
c output array:
      DIMENSION  rlat(3,*)
c scratch arrays:
      DIMENSION  rl1c(*), rl2c(*)
c local declarations:
      DIMENSION  g(3,3),rl(3),rcm(3),ruc(3), nnp(3),rt(3),rt1c(3)
      LOGICAL    sorted
      DATA  zero,half / 0.d0,0.5d0 /
      DATA  lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      almin = ala(1,1,1)
      do  ityp=1,ntyp
        do  ishl=1,numshl(ityp)
          do  ial=1,nala(ishl,ityp)
            if( ala(ial,ishl,ityp).lt.almin ) almin = ala(ial,ishl,ityp)
          enddo
        enddo
      enddo
c
c Compute linear independent length of each primitive lattice vector
c   (consider the effect of reduced dimensionality)
c
      call CROSS( rprim(1,1),rprim(1,2), g(1,3) )
      call CROSS( rprim(1,3),rprim(1,1), g(1,2) )
      call CROSS( rprim(1,2),rprim(1,3), g(1,1) )
      if( ndim.lt.3 )then
        call CROSS( g(1,3),rprim(1,1), g(1,2) )
        call CROSS( rprim(1,2),g(1,3), g(1,1) )
        if( ndim.lt.2 ) call CROSS( g(1,2),g(1,3), g(1,1) )
      endif
      do  idim=1,3
        gmag = SQRT(g(1,idim)**2+g(2,idim)**2+g(3,idim)**2)
        g(1,idim) = g(1,idim)/gmag
        g(2,idim) = g(2,idim)/gmag
        g(3,idim) = g(3,idim)/gmag
        rl(idim) = g(1,idim)*rprim(1,idim) + g(2,idim)*rprim(2,idim)
     $           + g(3,idim)*rprim(3,idim)
        rl(idim) = ABS(rl(idim))
      enddo
c
c Get center-of-"existence" of center coordinates:
      atmn = DBLE(natm)
      do  j=1,3
        ruc(j) = half*( rprim(j,1) + rprim(j,2) + rprim(j,3) )
        rcm(j) = zero
        do  iatm=1,natm
          rcm(j) = rcm(j) + ratm(j,iatm)
        enddo
        rcm(j) = rcm(j)/atmn
      enddo
c Get radius of a sphere/circle/line sure to contain all atoms:
      radatm = zero
      do  iatm=1,natm
        rt(1) = rcm(1) - ratm(1,iatm)
        rt(2) = rcm(2) - ratm(2,iatm)
        rt(3) = rcm(3) - ratm(3,iatm)
        do  idim=ndim+1,3
          cg = g(1,idim)*rt(1) + g(2,idim)*rt(2) + g(3,idim)*rt(3)
          rt(1) = rt(1) - cg*g(1,idim)
          rt(2) = rt(2) - cg*g(2,idim)
          rt(3) = rt(3) - cg*g(3,idim)
        enddo
        rmag = SQRT( rt(1)**2 + rt(2)**2 + rt(3)**2 )
        if( rmag.gt.radatm ) radatm = rmag
      enddo
c Bias center-of-mass:
      do  idim=1,ndim
        do  j=1,3
          rcm(j) = rcm(j) - 0.001*rprim(j,idim)
        enddo
      enddo
c
c find longest diagonal within unit cell
      raduc = zero
      if( ndim.eq.1 )then
        raduc = half*SQRT(rprim(1,1)**2+rprim(2,1)**2+rprim(3,1)**2)
      elseif( ndim.eq.2 )then
        rlp = SQRT( (rprim(1,1)+rprim(1,2))**2 +
     $   (rprim(2,1)+rprim(2,2))**2 + (rprim(3,1)+rprim(3,2))**2 )
        rlm = SQRT( (rprim(1,1)-rprim(1,2))**2 +
     $   (rprim(2,1)-rprim(2,2))**2 + (rprim(3,1)-rprim(3,2))**2 )
        raduc = half* MAX( rlp, rlm )
      elseif( ndim.eq.3 )then
        rt(1) = half*( rprim(1,1) + rprim(1,2) + rprim(1,3) )
        rt(2) = half*( rprim(2,1) + rprim(2,2) + rprim(2,3) )
        rt(3) = half*( rprim(3,1) + rprim(3,2) + rprim(3,3) )
        sqpp = zero
        sqpm = zero
        sqmp = zero
        sqmm = zero
        do  j=1,3
          sqpp = sqpp +   rt(j)**2
          sqpm = sqpm + ( rt(j) - rprim(j,3) )**2
          sqmp = sqmp + ( rt(j) - rprim(j,2) )**2
          sqmm = sqmm + ( rt(j) - rprim(j,2) - rprim(j,3) )**2
        enddo
        raduc = MAX( sqpp, sqpm, sqmp, sqmm )
        raduc = SQRT(raduc)
      endif
c
      rang1c = SQRT( cutgrid / (almin+almin) )
      rang2c = 2*SQRT( cutslow / (almin+almin) )
      r1cut = rang1c + radatm + raduc
      r2cut = rang2c + radatm + radatm
      l1 = 0
      l2 = 0
      l3 = 0
      goto (70,71,72,73), ndim+1
   73 l3 = (rang2c/rl(3)) + 1
   72 l2 = (rang2c/rl(2)) + 1
   71 l1 = (rang2c/rl(1)) + 1
   70 continue
c
c  Form trial linear combinations of primitive vectors
c
      ilat = 0
      do 103 i3=-l3,l3
        nnp(3) = i3
        do 102 i2=-l2,l2
          nnp(2) = i2
          do 101 i1=-l1,l1
            nnp(1) = i1
c
            do  j=1,3
              rt(j) = zero
            enddo
            do  j=1,3
              do  i=1,ndim
                rt(j) = rt(j) + rprim(j,i)* DBLE(nnp(i))
              enddo
            enddo
c bias results by center-of-mass for 1c range
            do  j=1,3
              rt1c(j) = rt(j) + rcm(j) - ruc(j)
            enddo
c
            r2c = SQRT(rt(1)**2+rt(2)**2+rt(3)**2)
            r1c = SQRT(rt1c(1)**2+rt1c(2)**2+rt1c(3)**2)
c this will force primary unit cell to the front
            if( i1.eq.0 .and. i2.eq.0 .and. i3.eq.0 ) r1c = -1.
            if( r1c.gt.r1cut .and. r2c.gt.r2cut ) goto 101
c
            ilat = ilat + 1
            if( ilat.gt.nlatd )then
c too many vectors, find out how many total before stopping
              call FLGETIERR( IERRFL )
              write(IERRFL,9100)  ilat,r2c,r2cut
 9100         format(' r2c(',i4,')=',f13.6,' while r2cut=',f13.6)
            else
c "rl1c" holds 1c rlat magnitudes, "rl2c" 2c mags:
              rl1c(ilat) = r1c
              rl2c(ilat) = r2c
              rlat(1,ilat) = rt(1)
              rlat(2,ilat) = rt(2)
              rlat(3,ilat) = rt(3)
            endif
c
  101     continue
  102   continue
  103 continue
      if( ilat.gt.nlatd )then
        write(IERRFL,*)'LATVEC: nlat= ',ilat,' exceeds max nlatd=',nlatd
        call STOPXERR( 'nlat-dim/exceeded max num of lattice vectors' )
      endif
      nlat = ilat
      nlat2c = nlat
c
c  Order lattice vectors wrt magnitude
c
      nlat1c = 0
  200 sorted = .true.
      do 220 ilat=nlat1c+1,nlat-1
        if( rl1c(ilat).le.rl1c(ilat+1) ) goto 220
c
        x = rl1c(ilat)
        rl1c(ilat) = rl1c(ilat+1)
        rl1c(ilat+1) = x
        x = rl2c(ilat)
        rl2c(ilat) = rl2c(ilat+1)
        rl2c(ilat+1) = x
        do  j=1,3
          rt(j) = rlat(j,ilat)
          rlat(j,ilat) = rlat(j,ilat+1)
          rlat(j,ilat+1) = rt(j)
        enddo
        sorted = .false.
  220 continue
      if( .not. sorted ) goto 200
c
      if( nlat1c.eq.0 )then
c       sort again, this time over 2c ranges
        do  ilat=1,nlat
          if( rl1c(ilat).le.r1cut ) nlat1c = ilat
          x = rl1c(ilat)
          rl1c(ilat) = rl2c(ilat)
          rl2c(ilat) = x
        enddo
        goto 200
      endif
c
      if( lstout.gt.0 )then
        call FLGETIWR( IWR )
        write(IWR,*) ' '
        write(IWR,*) 'LATVEC: nlat1c,nlat2c,nlat =',nlat1c,nlat2c,nlat
        if( lstout.gt.1 )then
          nlatprt = MIN( 160, nlat )
          write(IWR,9260)
 9260     format(2x,'ilat',10x,'rlat length',23x,'( lattice vector )')
 9261     format(1x,  i4,  2x,      2f12.6,5x,'('  ,3(f12.6,1x),    ')')
          do  ilat=1,nlatprt
            write(IWR,9261)  ilat, rl1c(ilat),rl2c(ilat),
     $       (rlat(i,ilat),i=1,3)
          enddo
          if( nlat.gt.nlatprt )then
            write(IWR,9261)  nlat, rl1c(ilat),rl2c(ilat),
     $       (rlat(i,ilat),i=1,3)
          endif
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
