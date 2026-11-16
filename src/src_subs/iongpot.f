c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IONGPOT
c
c
      subroutine IONGPOT( ndim, esion, nglm,qlmcc,rchrg,
     $ rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $ rhion, vion, r )
c---------------------------------------------------------------
c Purpose: Compute densities and potentials on grid from local
c          moment countercharge (LMCC) expressed as an array
c          of spherical gaussians (cluster or periodic).
c
c Written: Peter A. Schultz, 29-March-1999 (for v2.32)
c
c Revision history:
c  30Jun11-PAS/2.62: cosmetics
c   8Mar02-PAS/2.52: clean unused variables/literal constants d0
c---------------------------------------------------------------
c
c notes:
c   nglm    = number of Gaussians in LMCC array
c   qlmcc() = charge in each Gaussian
c   rchrg() = position of each Gaussian
c   esion   = self-energy of Gaussian LMCC array
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  rchrg(3,*), qlmcc(*)
c   grid generating primitive vector, and unit cell vectors
      DIMENSION  rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3)
c input/output arrays:
      DIMENSION  rhion(*),vion(*)
c scratch arrays:
      DIMENSION  r(3,*)
c
c local declarations:
      LOGICAL    ingaus
      DIMENSION  r0box(3),boxctr(3)
      DIMENSION  origlm(3), t(3)
      DIMENSION  gg(3,3), r1(3),r2(3),r3(3),rl(3,27)
      DIMENSION  nro(3),nrx(3)
c
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  three,four / 3.d0,4.d0 /
c  Exponent for Gaussian neutralizing charge:
      DATA  alfion / 0.30d0 /
c   Should be small enough to be well integrable on any reasonable
c   grid but large enough to be contained within any viable unit cell.
c
c >>>> EXECUTABLE CODE:
c
      esion = zero
c
c  Get range of this Gaussian on the grid:
      rngalf = SQRT( cutgrid/alfion )
      pi = four*ATAN( one )
      rtalf = SQRT( alfion )
      aldpi = alfion/pi
c
      do 1000 iglm=1,nglm
c
c Specify charge and position of this Gaussian:
        qgaus = qlmcc(iglm)
        do  j=1,3
          origlm(j) = rchrg(j,iglm)
        enddo
c
        gcoef = qgaus*aldpi* SQRT( aldpi )
        vcoef = two*qgaus
        vcoef0 = two*vcoef * SQRT( aldpi )
c
        if( ndim .eq. 0 )then
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           We have a finite cluster/molecular system
c  I can use the raw parallelpiped unit cell as my local volume
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Find position of coordinate origin on the grid:
          do  j=1,3
            r1(j) = origlm(j)
          enddo
c
          call PROJLV( r1, hh(1,1), hh(1,2), hh(1,3), t )
          do  id=1,3
            nro(id) = ( t(id) + half + one )
          enddo
c
          call CROSS( hh(1,2),hh(1,3), gg(1,1) )
          call CROSS( hh(1,3),hh(1,1), gg(1,2) )
          call CROSS( hh(1,1),hh(1,2), gg(1,3) )
          do  j=1,3
            gmag = SQRT( gg(1,j)**2 + gg(2,j)**2 + gg(3,j)**2 )
            gg(1,j) = gg(1,j)/gmag
            gg(2,j) = gg(2,j)/gmag
            gg(3,j) = gg(3,j)/gmag
          enddo
c
          do  id=1,3
            rhh = ABS( gg(1,id)*hh(1,id) + gg(2,id)*hh(2,id)
     $               + gg(3,id)*hh(3,id) )
            rnx = rngalf/rhh + half
            nrx(id) = rnx
          enddo
c
          ia1 = nro(1) - nrx(1)
          ian = nro(1) + nrx(1)
          ib1 = nro(2) - nrx(2)
          ibn = nro(2) + nrx(2)
          ic1 = nro(3) - nrx(3)
          icn = nro(3) + nrx(3)
          if( ia1.lt.1 ) ia1 = 1
          if( ian.gt.n1r ) ian = n1r
          if( ib1.lt.1 ) ib1 = 1
          if( ibn.gt.n2r ) ibn = n2r
          if( ic1.lt.1 ) ic1 = 1
          if( icn.gt.n3r ) icn = n3r
c
          esiongr = zero
          irpot = 0
          do 200 ibox=1,nbox
c
c  Get basic box info
            call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $       ira0,irb0,irc0, idra,idrb,idrc )
            nrbox = idra*idrb*idrc
c
            iv0 = irpot
            irpot = irpot + nrbox
c
            ingaus = .true.
            if( (ira0+idra).lt.ia1 .or. ian.le.ira0 .or.
     $          (irb0+idrb).lt.ib1 .or. ibn.le.irb0 .or.
     $          (irc0+idrc).lt.ic1 .or. icn.le.irc0 ) ingaus = .false.
c
c  Generate coordinates
c
            call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
            r0box(1) = r0box(1) - origlm(1)
            r0box(2) = r0box(2) - origlm(2)
            r0box(3) = r0box(3) - origlm(3)
            call MKRBOX( hh,r0box, idra,idrb,idrc, r )
c
            if( ingaus )then
c             This box is within the gaussian range
              do 160 nr=1,nrbox
                rsq = r(1,nr)**2 + r(2,nr)**2 + r(3,nr)**2
                rmag = SQRT( rsq )
c
                rhion(iv0+nr) = rhion(iv0+nr) + gcoef*EXP( -alfion*rsq )
                if( rmag .gt. 1.d-4 )then
                  vion(iv0+nr) = vion(iv0+nr)
     $             + vcoef*DERF( rtalf*rmag ) / rmag
                else
                  vion(iv0+nr) = vion(iv0+nr)
     $             + vcoef0*(one-alfion*rsq/three)
                endif
                esiongr = esiongr + vion(iv0+nr)*rhion(iv0+nr)
  160         continue
            else
c             Should replace with asymptotic limit between:
              do 180 nr=1,nrbox
                rsq = r(1,nr)**2 + r(2,nr)**2 + r(3,nr)**2
                rmag = SQRT( rsq )
                vion(iv0+nr) = vion(iv0+nr) + vcoef / rmag
                esiongr = esiongr + vion(iv0+nr)*rhion(iv0+nr)
  180         continue
            endif
c
  200     continue
c
c Finalize grid self-energy term for background charge:
          esiongr = weight*esiongr/two
c
c Compute analytic self-energy term for background charge:
          esiong = qgaus*qgaus* SQRT( two*alfion/pi )
c Since put and take is on grid, use grid integral as self-energy:
c  (for reasonable grids and alpha, should not make a difference)
          esion = esiongr
c
        else
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     We have a system that is periodic in at least one dimension
c   Must use Wigner-Seitz cell as my local volume (in periodic dims)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Note this code would work for ndim=0 cluster as well
c
          nla = 0
          nlb = 0
          nlc = 0
          if( ndim.gt.0 ) nla = 1
          if( ndim.gt.1 ) nlb = 1
          if( ndim.gt.2 ) nlc = 1
c use a nearest neighbor search to find longest range for which you do
c not need to check for residence within wigner-seitz cell
          rngws = 9.9d9
          il = 0
          do 343 kc=-nlc,nlc
            r3(1) = kc*rprim(1,3)
            r3(2) = kc*rprim(2,3)
            r3(3) = kc*rprim(3,3)
            do 342 kb=-nlb,nlb
              r2(1) = r3(1) + kb*rprim(1,2)
              r2(2) = r3(2) + kb*rprim(2,2)
              r2(3) = r3(3) + kb*rprim(3,2)
              do 341 ka=-nla,nla
                il = il + 1
                rl(1,il) = r2(1) + ka*rprim(1,1)
                rl(2,il) = r2(2) + ka*rprim(2,1)
                rl(3,il) = r2(3) + ka*rprim(3,1)
                rsqlat = rl(1,il)**2 + rl(2,il)**2 + rl(3,il)**2
                if( ka.eq.0 .and. kb.eq.0 .and. kc.eq.0 )then
                  ilo = il
                elseif( rsqlat.lt.rngws )then
                  rngws = rsqlat
                endif
  341         continue
  342       continue
  343     continue
          rngws = half*SQRT( rngws )
          rsqws = rngws**2
c
c  Put lattice vector zero to front of lattice vector list:
          nlnn = il
          rl(1,ilo) = rl(1,1)
          rl(2,ilo) = rl(2,1)
          rl(3,ilo) = rl(3,1)
          rl(1,1) = zero
          rl(2,1) = zero
          rl(3,1) = zero
c
          esiongr = zero
          irpot = 0
          do 400 ibox=1,nbox
c
c  Get basic box info
            call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $       ira0,irb0,irc0, idra,idrb,idrc )
            nrbox = idra*idrb*idrc
c
            iv0 = irpot
            irpot = irpot + nrbox
c
c  Generate coordinates
c
            call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
            r0box(1) = r0box(1) - origlm(1)
            r0box(2) = r0box(2) - origlm(2)
            r0box(3) = r0box(3) - origlm(3)
            call MKRBOX( hh,r0box, idra,idrb,idrc, r )
c
            ilmn = 1
c
            do 380 nr=1,nrbox
              rsq = (r(1,nr)+rl(1,ilmn))**2 + (r(2,nr)+rl(2,ilmn))**2
     $            + (r(3,nr)+rl(3,ilmn))**2
c
              if( rsq.gt.rsqws .and. nlnn.gt.1 )then
c               Check if some periodic lattice image of charge is closer
                do 360 il=1,nlnn
                  rsql = (r(1,nr)+rl(1,il))**2 + (r(2,nr)+rl(2,il))**2
     $                 + (r(3,nr)+rl(3,il))**2
                  if( rsql.lt.rsq )then
                    rsq = rsql
                    ilmn = il
                  endif
  360           continue
              endif
              rmag = SQRT( rsq )
c
              if( rmag.lt.rngalf )then
c               This point is within the gaussian countercharge range
                rhion(iv0+nr) = rhion(iv0+nr) + gcoef*EXP( -alfion*rsq )
                if( rmag .gt. 1.d-4 )then
                  vion(iv0+nr) = vion(iv0+nr)
     $             + vcoef*DERF( rtalf*rmag ) / rmag
                else
                  vion(iv0+nr) = vion(iv0+nr)
     $             + vcoef0*(one-alfion*rsq/three)
                endif
              else
c               Outside countercharge range, but still have potential
                vion(iv0+nr) = vion(iv0+nr) + vcoef / rmag
              endif
              esiongr = esiongr + vion(iv0+nr)*rhion(iv0+nr)
c
  380       continue
c
  400     continue
c
          esiongr = weight*esiongr/two
c
c Compute analytic self-energy term for background charge:
          esiong = qgaus*qgaus* SQRT( two*alfion/pi )
c Since put and take is on grid, use grid integral as self-energy:
c  (for reasonable grids and alpha, should not make a difference)
          esion = esiongr
c
        endif
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
