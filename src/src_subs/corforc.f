c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CORFORC
c
c
      subroutine CORFORC( IWR, ndim, frccor, strcor,
     $ natm,ntyp,nrd, nlat,  itypa,numshl, ratm,rlat,
     $ hh, weight, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ nrcor, radmsh,corden, y2,u,
     $ vxcgrd, r, wkb )
c---------------------------------------------------------------
c Purpose: compute parts of force due to partial core charge
c
c Written: Peter A. Schultz, 29-August-1995
c
c Revision history:
c  24May00-PAS/2.43h: bugfix of ndim=1,2 code
c  21Dec99-PAS/2.42: cosmetic clean-up
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  17Dec98-PAS/2.30: clean force/stress arrays
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  30Jan98-PAS/2.21: adapted for boxed grid
c   1May97-PAS/2.17: naming conventions, new grid coord rtn
c  15Dec95-PAS/2.14: shell-by-shell ranges
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  frccor(3,natm), strcor(3,3)
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp)
      DIMENSION  ratm(3,natm),rlat(3,nlat), hh(3,3)
      DIMENSION  ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  nrcor(ntyp), radmsh(nrd,*),corden(nrd,*)
      DIMENSION  vxcgrd(*)
c scratch arrays:
      DIMENSION  y2(*),u(*)
      DIMENSION  r(3,*), wkb(*)
c local declarations:
      DIMENSION  str(3,3)
      DIMENSION  ratom(3)
      DIMENSION  r0box(3),boxctr(3)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      call MKZERO( 9, strcor )
      call MKZERO( 3*natm, frccor )
c
c >>>>>>>>>> Set up grid:
c
      irpot = 0
      do 30 ibox=1,nbox
c
c  Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
        iv0 = irpot
        irpot = irpot + nrbox
c
        nobox = noboxs(ibox)
        if( nobox.le.0 ) goto 30
c
c  Generate coordinates
c
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
        call MKRBOX( hh,r0box, idra,idrb,idrc, r(1,iv0+1) )
c
   30 continue
c
c >>>> Atom loop:
c
      inext = 1
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
c
        irec0 = irec00
        inbox0 = inbox00
        iatshl0 = iatshl00
        do  ishl=1,numshl(ityp)
          inbox00 = 2*inbox00
          if( inrec*(iatshl00/inrec).eq.iatshl00 )then
            irec00 = irec00 + 1
            inbox00 = 1
          endif
          iatshl00 = iatshl00 + 1
        enddo
c
        nrsc = nrcor(ityp) + 1
c check if this center has any core charge associated with it
        if( nrsc.eq.1 ) goto 1000
c
c Interpolate core density
        ddy0 = 1.d31
        ddy1 = 1.d31
        call SPLINE( nrsc, radmsh(1,ityp),corden(1,ityp),
     $   ddy0,ddy1,y2,u )
        rngcor = radmsh(nrsc,ityp)
c
        fx = zero
        fy = zero
        fz = zero
        call MKZERO( 9, str )
c
c  Cycle over boxes:
c
        irpot = 0
        do 300 ibox=1,nbox
c
c  Get basic box info
          call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $     ira0,irb0,irc0, idra,idrb,idrc )
          nrbox = idra*idrb*idrc
c
          iv0 = irpot
          irpot = irpot + nrbox
c
          nobox = noboxs(ibox)
          if( nobox.le.0 ) goto 300
c
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
c
c  Set up for scratch space below:
          ir0 = 0
          irmag = ir0
          idcor = irmag + nrbox
c
c >>>> Lattice loop:
c
          do 900 ilat=1,nlat
c
            rngbox = zero
            do  jr=1,3
              ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
              rngbox = rngbox + (ratom(jr)-boxctr(jr))**2
            enddo
            rngbox = SQRT(rngbox) - boxrad
            if( rngbox.gt.rngcor ) goto 900
c
c  Generate r**2 on mesh ...
            do  nr=1,nrbox
              wkb(irmag+nr) = SQRT( (r(1,iv0+nr)-ratom(1))**2
     $                            + (r(2,iv0+nr)-ratom(2))**2
     $                            + (r(3,iv0+nr)-ratom(3))**2 )
            enddo
c
c interpolate to get deriv. of core density from radial core density
            call SPLIND( nrsc, radmsh(1,ityp),corden(1,ityp),y2,
     $                   nrbox,wkb(irmag+1),wkb(idcor+1) )
c
c to prevent divide-by-zeroes
            do  nr=1,nrbox
              if( wkb(irmag+nr) .lt. 1.d-8 ) wkb(irmag+nr) = 1.d-8
            enddo
c
            if( ndim.eq.0 )then
              do  nr=1,nrbox
                dx = r(1,iv0+nr) - ratom(1)
                dy = r(2,iv0+nr) - ratom(2)
                dz = r(3,iv0+nr) - ratom(3)
                rfac = wkb(idcor+nr)*vxcgrd(iv0+nr)/wkb(irmag+nr)
                fx = fx + rfac*dx
                fy = fy + rfac*dy
                fz = fz + rfac*dz
              enddo
            elseif( ndim.eq.1 )then
              do  nr=1,nrbox
                dx = r(1,iv0+nr) - ratom(1)
                dy = r(2,iv0+nr) - ratom(2)
                dz = r(3,iv0+nr) - ratom(3)
                rfac = wkb(idcor+nr)*vxcgrd(iv0+nr)/wkb(irmag+nr)
                fx = fx + rfac*dx
                fy = fy + rfac*dy
                fz = fz + rfac*dz
                str(1,1) = str(1,1) + rfac*dx*dx
              enddo
            elseif( ndim.eq.2 )then
              do  nr=1,nrbox
                dx = r(1,iv0+nr) - ratom(1)
                dy = r(2,iv0+nr) - ratom(2)
                dz = r(3,iv0+nr) - ratom(3)
                rfac = wkb(idcor+nr)*vxcgrd(iv0+nr)/wkb(irmag+nr)
                fx = fx + rfac*dx
                fy = fy + rfac*dy
                fz = fz + rfac*dz
                str(1,1) = str(1,1) + rfac*dx*dx
                str(2,1) = str(2,1) + rfac*dy*dx
                str(1,2) = str(1,2) + rfac*dx*dy
                str(2,2) = str(2,2) + rfac*dy*dy
              enddo
            elseif( ndim.eq.3 )then
              do  nr=1,nrbox
                dx = r(1,iv0+nr) - ratom(1)
                dy = r(2,iv0+nr) - ratom(2)
                dz = r(3,iv0+nr) - ratom(3)
                rfac = wkb(idcor+nr)*vxcgrd(iv0+nr)/wkb(irmag+nr)
                fx = fx + rfac*dx
                fy = fy + rfac*dy
                fz = fz + rfac*dz
                str(1,1) = str(1,1) + rfac*dx*dx
                str(2,1) = str(2,1) + rfac*dy*dx
                str(3,1) = str(3,1) + rfac*dz*dx
                str(1,2) = str(1,2) + rfac*dx*dy
                str(2,2) = str(2,2) + rfac*dy*dy
                str(3,2) = str(3,2) + rfac*dz*dy
                str(1,3) = str(1,3) + rfac*dx*dz
                str(2,3) = str(2,3) + rfac*dy*dz
                str(3,3) = str(3,3) + rfac*dz*dz
              enddo
            endif
c
c Close lattice vector loop:
c
  900     continue
c
c Close box loop:
c
  300   continue
c
        frccor(1,iatm) = frccor(1,iatm) + weight*fx
        frccor(2,iatm) = frccor(2,iatm) + weight*fy
        frccor(3,iatm) = frccor(3,iatm) + weight*fz
        do  jd=1,3
          do  id=1,3
             strcor(id,jd) = strcor(id,jd) + weight*str(id,jd)
          enddo
        enddo
c
c Close atom loop:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
