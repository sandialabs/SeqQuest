c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CNTRRHO
c
c
      subroutine CNTRRHO( icountr,
     $ natm,ntyp,nshld, n1r,n2r,n3r,nptr,nlat,maxwf,
     $ itypa, numshl,lshel, ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ r, rv, rho )
c---------------------------------------------------------------
c Purpose: compute charge densities of the required angular
c          symmetries used to cancel the lowest multipole
c          moments of the 1-center charges in the solution
c          of Poisson's equation for strong pseudopot atoms
c
c Written: Peter A. Schultz, based on PJFs original routine.
c
c Revision history:
c   8Mar02-PAS/2.52: clean unused variables/literal constants d0
c  21Jun01-PAS/2.48: replace STOPs
c  19Dec99-PAS/2.42: clean-up; read/write big grid records
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  24Apr97-PAS/2.17: naming conventions + on/off criterion
c  15Dec95-PAS/2.14: shell-by-shell range + cleanup
c  12Oct93-PAS: loosen a bit for rough consistency with rho calc
c  17Sep93-PAS: more stringent cutoffs installed
c   6Apr90-PAS: made more efficent and to need less disk
c---------------------------------------------------------------
c
c  On input, mesh coordinates  assumed carried in r()
c  Should install z-cutoffs of shorter range as well
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(ntyp), numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  ratm(3,natm),hh(3,3),rlat(3,nlat)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
c input/output:
      DIMENSION  r(3,nptr)
c scratch:
      DIMENSION  rv(3,nptr)
      DIMENSION  rho(*)
c local declarations:
      DIMENSION  number(3),ratom(3)
      DATA  number / 1,4,10 /
      DATA  zero,half  / 0.d0,0.5d0 /
c
c >>>> EXECUTABLE CODE:
c
      argmx = half*cutgaus
c
c    Loop over centers:
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
        if( lmx1ctr(ityp).lt.0 ) goto 1000
c
        al1c = almnv1c(ityp)
        twoal1c = al1c + al1c
c this cutoff ill-defined with current wvfcn cutoff implementation
        rsqcut = cutgrid / (0.75d0*al1c)
c       rsqcut = cutgrid / (0.5d0*(al1c+almin))
c
        aldpi = al1c/pi
        scoef = SQRT( aldpi )*aldpi
        pcoef = twoal1c*scoef
        dcoef = twoal1c*pcoef
        numlmx = number(lmx1ctr(ityp)+1)
c
        nr0 = n1r*n2r*(imin3-1)
        nrset = numlmx*nrat
        if( (nrset+5*nrat).gt.2*maxwf ) call STOPXERR( 'cntrrho1' )
c
        do  nr=1,nrset
          rho(nr) = zero
        enddo
        ng = nrset
        nrsq = ng + nrat
        nx = nrsq + nrat
        ny = nx + nrat
        nz = ny + nrat
c
        do 900 ilat=1,nlat
c
          do  jr=1,3
            ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
          enddo
c
          do  nr=1,nrat
            rho(ng+nr) = zero
            rho(nx+nr) = r(1,nr0+nr) - ratom(1)
            rho(ny+nr) = r(2,nr0+nr) - ratom(2)
            rho(nz+nr) = r(3,nr0+nr) - ratom(3)
            rho(nrsq+nr) = rho(nx+nr)**2 + rho(ny+nr)**2 + rho(nz+nr)**2
          enddo
          rsqmin = rho(nrsq+1)
          rsqmax = rho(nrsq+1)
          do  nr=1,nrat
            if( rho(nrsq+nr).lt.rsqmin ) rsqmin = rho(nrsq+nr)
            if( rho(nrsq+nr).gt.rsqmax ) rsqmax = rho(nrsq+nr)
          enddo
          if( rsqmin.gt.rsqcut ) goto 900
c
          if( al1c*rsqmax.lt.argmx )then
            do  nr=1,nrat
              rho(ng+nr) = EXP( -al1c*rho(nrsq+nr) )
            enddo
          else
            do  nr=1,nrat
              if( rho(nrsq+nr).lt.rsqcut )
     $        rho(ng+nr) = EXP( -al1c*rho(nrsq+nr) )
            enddo
          endif
c
          if( numlmx.eq.1 )then
c
            do  nr=1,nrat
              rho(nr) = rho(nr) + scoef*rho(ng+nr)
            enddo
c
          elseif( numlmx.eq.4 )then
c
            do  nr=1,nrat
              gauss = pcoef*rho(ng+nr)
              rho(       nr) = rho(       nr) + gauss*(scoef/pcoef)
              rho(  nrat+nr) = rho(  nrat+nr) + gauss*rho(nx+nr)
              rho(2*nrat+nr) = rho(2*nrat+nr) + gauss*rho(ny+nr)
              rho(3*nrat+nr) = rho(3*nrat+nr) + gauss*rho(nz+nr)
            enddo
c
          elseif( numlmx.eq.10 )then
c
            do  nr=1,nrat
              x = rho(nx+nr)
              y = rho(ny+nr)
              z = rho(nz+nr)
              gauss = rho(ng+nr)
              rho(       nr) = rho(       nr) + scoef*gauss
              rho(  nrat+nr) = rho(  nrat+nr) + pcoef*gauss* x
              rho(2*nrat+nr) = rho(2*nrat+nr) + pcoef*gauss* y
              rho(3*nrat+nr) = rho(3*nrat+nr) + pcoef*gauss* z
              rho(4*nrat+nr) = rho(4*nrat+nr) + dcoef*gauss* x*y
              rho(5*nrat+nr) = rho(5*nrat+nr) + dcoef*gauss* y*z
              rho(6*nrat+nr) = rho(6*nrat+nr) + dcoef*gauss* z*x
              rho(7*nrat+nr) = rho(7*nrat+nr) + dcoef*gauss* x*x
              rho(8*nrat+nr) = rho(8*nrat+nr) + dcoef*gauss* y*y
              rho(9*nrat+nr) = rho(9*nrat+nr) + dcoef*gauss* z*z
            enddo
c
          endif
c
c    Close lattice vector loop:
c
  900   continue
c
c    Write out counter-charges:
c
        nrset0 = 0
        do  lmx=1,numlmx
C         write(icountr)  (rho(nr),nr=nrset0+1,nrset0+nrat)
          call WRITBIG( icountr, nrat, rho(nrset0+1) )
          nrset0 = nrset0 + nrat
        enddo
c
c    Close center loop:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
