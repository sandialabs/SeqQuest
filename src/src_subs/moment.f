c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MOMENT
c
c
      subroutine MOMENT( IWR, elchrg,dipole,rdipole, ndim,
     $ origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $ rho, r )
c---------------------------------------------------------------
c Purpose: compute the dipole moment of a charge distribution
c          known on the grid
c
c Written: Peter A. Schultz, 28-October-1998, for v2.27 (ions)
c
c Revision history:
c  14Jun00-PAS/2.45: deal with grid defect charge and dipoles
c  12Jun00-PAS/2.44: bias latvec/ws to make slab match bulk
c   9Mar99-PAS/2.32: boxed wigner-seitz periodic supercell dipole
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Nov98-PAS/2.27: compute and send position of cluster dipole
c---------------------------------------------------------------
c
c Input variables:
c  origws() = center of the wigner-seitz (WS) cell
c  rprim()  = primitive lattice vectors
c  hh()     = primitive vectors of grid
c  n1r,n2r,n3r = dimension of grid in three dimensions
c  ndbox()  = defines box-wise breakup of grid
c  rho()    = input electron density on grid, sorted into boxes
c
c Output variables:
c  dipole() = computed dipole moment
c  rdipole() = computed dipole position
c
c Notes:
c  This code has a zero bias, i.e., points on/near the boundary
c  contribute only to the cell with that point on the zero side
c  of the cell.  However, this "point" represents a *volume*,
c  which should be shared among the cells it overlaps.  Here, I
c  assume that the supercell approximation is sound, and that
c  the delta density at this boundary is zero.  This should be
c  dealt with properly - this ideal is never truly achieved
c  in charged systems as the density polarizes due to the charge.
c
c  Ouch.  Because of grid-inaccuracies, the density on the grid
c  does not exactly integrate to the correct number of electrons.
c  Need to address this here explicitly somehow as calculation
c  is affected by this: computed dipole depends on frame of
c  reference with net q.  The ad hoc solution here is to take the
c  defect charge off its charge centroid and compute dipole from
c  that.  Increasing vacuum region and changing reference point
c  then yields the same computed dipole.  14Jun00-PAS/2.44
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  origws(3)
c   grid generating primitive vector, and unit cell vectors
      DIMENSION  rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3)
      DIMENSION  rho(*)
c output arrays:
      DIMENSION  dipole(3), rdipole(3)
c scratch arrays:
      DIMENSION  r(3,*)
c local declarations:
      DIMENSION  ctrpos(3),ctrneg(3)
      DIMENSION  r0box(3),boxctr(3)
      DIMENSION  orig(3), r1(3),r2(3),r3(3),rl(3,27)
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c vsmall is diff between R**2's considered to be equal.
      DATA  vsmall / 0.9d-6 /
      DATA  lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      nla = 1
      nlb = 1
      nlc = 1
      if( ndim.gt.0 ) nla = 3
      if( ndim.gt.1 ) nlb = 3
      if( ndim.gt.2 ) nlc = 3
c use a nearest neighbor search to find longest range for which you do
c not need to check for residence within wigner-seitz cell
      rngws = 9.9d9
      il = 0
c
      isc = 1
      kc = 0
      do 343 ilc=1,nlc
        r3(1) = kc*rprim(1,3)
        r3(2) = kc*rprim(2,3)
        r3(3) = kc*rprim(3,3)
        kc = kc + isc*ilc
        isc = -isc
c
        isb = 1
        kb = 0
        do 342 ilb=1,nlb
          r2(1) = r3(1) + kb*rprim(1,2)
          r2(2) = r3(2) + kb*rprim(2,2)
          r2(3) = r3(3) + kb*rprim(3,2)
          kb = kb + isb*ilb
          isb = -isb
c
          isa = 1
          ka = 0
          do 341 ila=1,nla
            il = il + 1
            rl(1,il) = r2(1) + ka*rprim(1,1)
            rl(2,il) = r2(2) + ka*rprim(2,1)
            rl(3,il) = r2(3) + ka*rprim(3,1)
            ka = ka + isa*ila
            isa = -isa
            rsqlat = rl(1,il)**2 + rl(2,il)**2 + rl(3,il)**2
            if( il.ne.1 .and. rsqlat.lt.rngws )then
              rngws = rsqlat
            endif
  341     continue
  342   continue
  343 continue
      rngws = half*SQRT( rngws )
      rsqws = rngws**2
c
c  Note: origin is at front of lattice vector list:
      nlnn = il
      rl(1,1) = zero
      rl(2,1) = zero
      rl(3,1) = zero
c
      sumrho = zero
      sumdif = zero
      sumpos = zero
      sumneg = zero
      do  j=1,3
        dipole(j) = zero
        ctrpos(j) = zero
        ctrneg(j) = zero
c Center calculation about center of WS cell:
        orig(j) = -origws(j)
      enddo
c
      irpot = 0
      do 380 ibox=1,nbox
c
c  Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
        iv0 = irpot
        irpot = irpot + nrbox
c
c  Generate coordinates
c
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
        r0box(1) = r0box(1) + orig(1)
        r0box(2) = r0box(2) + orig(2)
        r0box(3) = r0box(3) + orig(3)
        call MKRBOX( hh,r0box, idra,idrb,idrc, r )
c
        if( ndim.eq.0 )then
c         Cluster: WS cell same as unit cell
          do 350 nr=1,nrbox
            r1(1) = r(1,nr)
            r1(2) = r(2,nr)
            r1(3) = r(3,nr)
c
            rh = rho(iv0+nr)
            sumrho = sumrho + rh
            dipole(1) = dipole(1) + rh*r1(1)
            dipole(2) = dipole(2) + rh*r1(2)
            dipole(3) = dipole(3) + rh*r1(3)
c
            sumdif = sumdif + ABS(rh)
            if( rh.lt.zero )then
              sumneg = sumneg - rh
              ctrneg(1) = ctrneg(1) - rh*r1(1)
              ctrneg(2) = ctrneg(2) - rh*r1(2)
              ctrneg(3) = ctrneg(3) - rh*r1(3)
            else
              sumpos = sumpos + rh
              ctrpos(1) = ctrpos(1) + rh*r1(1)
              ctrpos(2) = ctrpos(2) + rh*r1(2)
              ctrpos(3) = ctrpos(3) + rh*r1(3)
            endif
  350     continue
c
        else
c         Periodic: WS not same as unit cell
          ilmn = 1
          do 370 nr=1,nrbox
            rsq = (r(1,nr)+rl(1,ilmn))**2 + (r(2,nr)+rl(2,ilmn))**2
     $          + (r(3,nr)+rl(3,ilmn))**2
c
            if( rsq.gt.rsqws .and. nlnn.gt.1 )then
c             Check if some periodic lattice image of charge is closer
              do 360 il=1,nlnn
                rsql = (r(1,nr)+rl(1,il))**2 + (r(2,nr)+rl(2,il))**2
     $               + (r(3,nr)+rl(3,il))**2
                if( ABS( rsql - rsq ) .lt. vsmall )then
c                 Should deal with points on/near boundaries better.
c                 This hack at least should make slab and bulk results
c                 match by biasing boundaries toward earlier latvecs.
                elseif( rsql.lt.rsq )then
                  rsq = rsql
                  ilmn = il
                endif
  360         continue
            endif
c
            r1(1) = r(1,nr) + rl(1,ilmn)
            r1(2) = r(2,nr) + rl(2,ilmn)
            r1(3) = r(3,nr) + rl(3,ilmn)
c
            rh = rho(iv0+nr)
            sumrho = sumrho + rh
            dipole(1) = dipole(1) + rh*r1(1)
            dipole(2) = dipole(2) + rh*r1(2)
            dipole(3) = dipole(3) + rh*r1(3)
c
            sumdif = sumdif + ABS(rh)
            if( rh.lt.zero )then
              sumneg = sumneg - rh
              ctrneg(1) = ctrneg(1) - rh*r1(1)
              ctrneg(2) = ctrneg(2) - rh*r1(2)
              ctrneg(3) = ctrneg(3) - rh*r1(3)
            else
              sumpos = sumpos + rh
              ctrpos(1) = ctrpos(1) + rh*r1(1)
              ctrpos(2) = ctrpos(2) + rh*r1(2)
              ctrpos(3) = ctrpos(3) + rh*r1(3)
            endif
  370     continue
c
        endif
c
c    Close up box loop:
c
  380 continue
c
      if( sumpos.ne.zero )then
        do  j=1,3
          ctrpos(j) = ctrpos(j)/sumpos
        enddo
      endif
      if( sumneg.ne.zero )then
        do  j=1,3
          ctrneg(j) = ctrneg(j)/sumneg
        enddo
      endif
c
c  Install grid weighting factors to get final result:
      sumrho = weight*sumrho
      dipole(1) = weight*dipole(1)
      dipole(2) = weight*dipole(2)
      dipole(3) = weight*dipole(3)
c
      sumdif = weight*sumdif
      sumpos = weight*sumpos
      sumneg = weight*sumneg
c
c  Take defect charge out of dipole calculation:
      qdefct = sumrho + elchrg
      if( qdefct.gt.zero )then
        dipole(1) = dipole(1) - qdefct*ctrpos(1)
        dipole(2) = dipole(2) - qdefct*ctrpos(2)
        dipole(3) = dipole(3) - qdefct*ctrpos(3)
      else
        dipole(1) = dipole(1) - qdefct*ctrneg(1)
        dipole(2) = dipole(2) - qdefct*ctrneg(2)
        dipole(3) = dipole(3) - qdefct*ctrneg(3)
      endif
c
      dipsq = zero
      do  idim=1,3
        dipsq = dipsq + dipole(idim)**2
      enddo
c
c
c Complete calculation of dipole location, and offset by the
c origin "origws()" of the dipole computation
      do  j=1,3
        rdipole(j) = ( ctrpos(j) + ctrneg(j) ) / two
     $   + origws(j)
      enddo
c
 9020 format(1x,a,2f20.10)
 9030 format(1x,a,3f15.8)
      if( lstout.gt.1 .and. dipsq.gt. 1.d-4 )then
        write(IWR,9020) 'MOM: ########### dipole analysis ###########'
        write(IWR,9020) 'MOM:  grid charge=',sumrho,qdefct
        write(IWR,9030) 'MOM:  grid dipole=', dipole
        write(IWR,9030) 'MOM:  dipole origin=',rdipole
        if( lstout.gt.3 )then
          write(IWR,9020) 'MOM:  grid charge movement=',sumdif
          write(IWR,9020) 'MOM:  grid positive charge=',sumpos
          write(IWR,9030) 'MOM:  positive centroid=',ctrpos
          write(IWR,9020) 'MOM:  grid negative charge=',sumneg
          write(IWR,9030) 'MOM:  negative centroid=',ctrneg
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
