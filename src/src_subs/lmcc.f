c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LMCC
c
c
      subroutine LMCC( IWR, igrdref,  lvlout, lmcc_on, ion_opt,
     $ esion,vacplus, eslmref, elchrg, ndim,nsym,
     $ dolvatm, rprim,rscale, orig,rsym,origws,rchrg, rmatsym,
     $ n1r,n2r,n3r,nptr, hh,weight, ndbox,nbox,
     $ wk1,wk2,wk3,wk4,wk5 )
c---------------------------------------------------------------
c Purpose: generate grid densities and potentials for the
c          local moment countercharge (LMCC) treatment of
c          supercell electrostatics.
c
c Written: Peter A. Schultz, 17-May-2001 for v2.47
c
c Revision history:
c  13May11-PAS/2.62: prepping inserting of 2D slab charges
c  15Feb11-PAS/2.62: extract slab dipole vacuum potential
c  10Sep05-PAS/2.59: chemical potential perspective diagnostics
c   4Oct03-PAS/2.57: set/require reference be non-spin; cosmetics
c   5May03-PAS/2.55: fix LMCC + spin bug (up+dn delrho)
c---------------------------------------------------------------
c
c Notes:
c  Implements the method for treating supercell moments described
c  in:  P.A. Schultz,  PRB 60, 1551 (1999);  PRL 84, 1942 (2000).
c
c on input
c   wk1 contains the electron density (with moments) on the grid
c       and is destroyed
c on output
c   wk2 contains the lmcc density on grid
c   wk3 contains the lmcc potential on grid
c Grid fields are presumed sorted into boxes, and returned boxed
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, half=0.5d0 )
c Input arrays:
      LOGICAL    lmcc_on, dolvatm
      DIMENSION  rprim(3,3),rscale(*)
      DIMENSION  orig(3),rsym(3),origws(3),rchrg(*)
      DIMENSION  rmatsym(*)
      DIMENSION  hh(3,3), ndbox(2,2,3)
c
c big arrays:
      DIMENSION  wk1(nptr),wk2(nptr),wk3(nptr),wk4(nptr),wk5(nptr)
c
c Local declarations:
      DIMENSION  dipole(3), rdipole(3)
      DIMENSION  dipref(3), rdipref(3)
      DIMENSION  vect(3), vec0(3)
      DIMENSION  zchrglm(3), zchrgcc(3)
cxxx
      LOGICAL    debugdl
c
      DATA  vec0 / 0.d0,0.d0,0.d0 /
      DATA  toev / 13.605805d0 /
c
c >>>> EXECUTABLE CODE:
c
cxxx:
      debugdl = .false.
      if( debugdl )then
        write(IWR,*) 'DBG/LMCC: esion, dble-layer self-energy - clean'
      endif
      lmcc_on = .false.
      esion = zero
      vacplus = zero
      eslmref = zero
c  Format for diagnostics:
 9020 format(1x,a,3f20.10)
cxxx:
 9013 format(1x,a,6f15.6)
c
c Start by taking moment of total density
c
      call MOMENT( IWR, elchrg, dipole,rdipole, ndim,
     $ origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $ wk1,    wk4 )
cxxx:start
      if( debugdl )then
        qgrid = zero
        do  i=1,nptr
          qgrid = qgrid + wk1(i)
        enddo
        qgrid = weight*qgrid
        write(IWR,9020) 'DBG/LMCC: delta-rho net charge=',qgrid
      endif
cxxx:end
c -->  rho-io  r-3bs
c
      call SYMVEC( 1, dipole, vec0, nsym, rmatsym )
      call SYMVEC( 1,rdipole, rsym, nsym, rmatsym )
cxxx:
      if( debugdl )then
        write(IWR,9013) 'DBG/LMCC: delrho dipole=',dipole
        write(IWR,9013) 'DBG/LMCC: delrho rdipole=',rdipole
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  If have reference density, retrieve and compute dipoles.
c  I will want two dipoles when I come out of the next section.
c  The first is a "defect" dipole for which a local spherical gaussian
c  countercharge array will be used, and the second is a "vacuum"
c  dipole with boundaries in, of course, vacuum.  Note that molecules
c  can be classed in either because they have boundaries in vacuum,
c  but they also use spherical gaussian countercharges.  Here, I class
c  them with the defect terms.  The defect dipole gets put into
c  the vector dipole(), while the vacuum dipole (now the only one
c  we consider is the slab dipole) gets put into dipref().
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( igrdref .gt. 0 )then
c
        write(IWR,*)
        if( lvlout.gt.1 )then
c         Report dipole moment from full delta-density of this system:
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     vect, rdipole )
          write(IWR,9020) 'full dipole moment(z*bohr)=',dipole
          write(IWR,9020) 'full dipole center: rdipole=',vect
        endif
c
c       Retrieve reference delta-density, and put into boxes:
        REWIND( unit=igrdref )
        call READBIG( igrdref, nptr, wk4 )
        call INTOBOX( n1r,n2r,n3r, wk4, ndbox, wk5 )
c
c       Compute moments from reference system density
c
        call MOMENT( IWR, elchrg, dipref,rdipref, ndim,
     $   origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $   wk5,    wk4 )
c -->    rho-io  r-3bs
c
        call SYMVEC( 1, dipref, vec0, nsym, rmatsym )
        call SYMVEC( 1,rdipref, rsym, nsym, rmatsym )
c
        if( lvlout.gt.1 )then
c         Report dipole from reference system density
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     vect, rdipref )
          write(IWR,9020) 'reference dipole(z*bohr)=',dipref
          write(IWR,9020) 'reference rdipole=',vect
        endif
c
c       Take difference of current density from reference:
        do  nr=1,nptr
          wk3(nr) = wk1(nr) - wk5(nr)
        enddo
c
c       Take moment of difference of density from reference:
c
        call MOMENT( IWR, elchrg, dipole,rdipole, ndim,
     $   origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $   wk3,    wk4 )
c -->    rho-io  r-3bs
c
        call SYMVEC( 1, dipole, vec0, nsym, rmatsym )
        call SYMVEC( 1,rdipole, rsym, nsym, rmatsym )
c
        if( lvlout.gt.1 )then
c         Report lmcc defect dipole:
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     vect, rdipole )
          write(IWR,9020) 'net dipole(z*bohr)=',dipole
          write(IWR,9020) 'net mom: rdipole=',vect
        endif
c
      else
c
        if( lvlout.gt.4 )then
c         Report full dipole from this system:
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     vect, rdipole )
          write(IWR,9020) 'full dipole moment(z*bohr)=',dipole
          write(IWR,9020) 'full dipole center: rdipole=',vect
        endif
c
c       Clear out undefined parts of dipole
        do  idim=1,ndim
          dipole(idim) = zero
          rdipole(idim) = zero
        enddo
c
c       In case have slab calculation, save moment in reference variables
        do  idim=1,3
          dipref(idim) = dipole(idim)
          rdipref(idim) = rdipole(idim)
        enddo
c
        if( ndim.eq.0 )then
c
c         Symmetrize net dipole
          call SYMVEC( 1, dipole, vec0, nsym, rmatsym )
          call SYMVEC( 1,rdipole, rsym, nsym, rmatsym )
          if( nsym.gt.1 .and. lvlout.gt.1 )then
            call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $       vect, rdipole )
            write(IWR,*)
            write(IWR,9020) 'net dipole(z*bohr)=',dipole
            write(IWR,9020) 'net mom: rdipole=',vect
          endif
c
        else
c
c         Clear out defect lmcc moment for non-molecules
c         (This compels skip of elchrg=0, point dipole LMCC)
          do  idim=1,3
            dipole(idim) = zero
            rdipole(idim) = zero
          enddo
        endif
c
      endif
c
c  Compute magnitude of lmcc "defect" dipole moment:
      dipsq = zero
      do  idim=1,3
        dipsq = dipsq + dipole(idim)**2
      enddo
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      We have all the dipole moments, now go and generate the
c      LMCC charge density and LMCC potential for defect LMCC
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Clear out local moment counter-charge density and potential
      call MKZERO( nptr, wk2 )
      call MKZERO( nptr, wk3 )
c
      if( elchrg .ne. zero )then
c       Charged system (point charge): compute LMCC density and potential
c
        if( ion_opt .eq. 2 )then
c         Initialize charge position to WS (dipole) origin
          do  i=1,3
            rchrg(i) = origws(i)
          enddo
        endif
c
        if( IABS( ion_opt ) .ne. 3 )then
c         (Shunt slab 2D planar charge, ion_opt=3, to later in routine)
          call IONPOT( ion_opt,elchrg,dipole,ndim,esion,rchrg,
     $     rprim, n1r,n2r,n3r,  hh,weight, ndbox,nbox,
     $     wk2,    wk3,    wk1 )
c -->      rhion-o vion-o  r-3bs
c
          if( lvlout.gt.0 .and. ion_opt.gt.0 )then
            call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $       vect  ,rchrg )
            write(IWR,*)
            write(IWR,9020) 'Location of LMCC charge=',vect
          endif
c
          lmcc_on = .true.
        endif
c
      elseif( dipsq .gt. 4.d-4 .and. ion_opt.gt.0 )then
c       Neutral system, with dipole: compute counter-moment charge density:
c       NB: Do not know how to specify dipole for point-dipole in bulk.
c
        if( ndim.eq.0 .or. ndim.eq.3 )then
c         (Shunt slab planar dipole treatment to later in routine)
          call DIPPOT( ion_opt, dipole,rdipole, ndim, esion,
     $     rprim, n1r,n2r,n3r,  hh,weight, ndbox,nbox,
     $     wk2,    wk3,    wk1 )
c -->      rhion-o vion-o  r-3bs
c
          if( lvlout.gt.0 )then
            write(IWR,*)
            write(IWR,9020) 'LMCC local dipole, |dip|=', SQRT( dipsq )
          endif
          lmcc_on = .true.
        endif
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Compute excluded volume reference term via fft/local b.c. exchange:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( igrdref.gt.0 .and. lmcc_on )then
c
c       Retrieve reference delta electrostatic potl; put into boxes:
        call READBIG( igrdref, nptr, wk1 )
        call INTOBOX( n1r,n2r,n3r, wk1, ndbox, wk4 )
c
c       First, a simple diagnostic check, v(ref):rho(ref)
        esnsref = half*weight*DDOT( nptr, wk4,1, wk5,1 )
        write(IWR,9020) 'Reference system esns(Ry)=',esnsref
c
c       Take G=0 net charge error out of reference delta density:
        rhave = zero
        do  nr=1,nptr
          rhave = rhave + wk5(nr)
        enddo
        rhave = rhave / DBLE( nptr )
        do  nr=1,nptr
          wk5(nr) = wk5(nr) - rhave
        enddo
c
c       Take pbc reference LMCC integral over all space: rho(ion):v(ref)
        epionvref = weight*DDOT( nptr, wk2,1, wk4,1 )
c
c       Take local reference LMCC integral over WS cell: v(ion)*rho(ref)
        evionpref = weight*DDOT( nptr, wk3,1, wk5,1 )
        eslmref = epionvref - evionpref
        eslmxev = toev*eslmref
        if( lvlout.gt.1 )then
          write(IWR,9020) 'LMCC reference int(Ry): pivx(pbc)=',
     $     epionvref
          write(IWR,9020) 'LMCC reference int(Ry): vipx(loc)=',
     $     evionpref
          write(IWR,9020) 'LMCC total reference int(Ry,eV): pivx-vipx=',
     $     eslmref,eslmxev
        endif
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Treat slab 2D planar dipole or charge via vacuum LMCC
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      esionsl = zero
      if( ndim.eq.2 .and. ion_opt .ne. 0 )then
c       Slab with 2D electrostatic corrections (i.e., planar, not point)
c
c       Compute magnitude of slab (normal) dipole:
        diprefsq = zero
        do  idim=ndim+1,3
          diprefsq = diprefsq + dipref(idim)**2
        enddo
c
        if( elchrg.eq.zero .and. diprefsq .gt. 4.d-4 )then
c         Neutral slab, with a planar dipole
c
          call DIPSLAB( ion_opt, ndim,
     $     dipref,rdipref,  esionsl,vacplus,
     $     rprim, n1r,n2r,n3r,  hh,weight, ndbox,nbox,
     $     wk2,     wk3,     wk1 )
c -->      rhlmcc-o vlmcc-o  r-3b
c
c         Let user know slab correction is on:
          if( ion_opt.eq.2 )then
            write(IWR,*) 'Slab dipole removed using vacuum LMCC'
            write(IWR,*) 'Vacuum level[Ry]@z(+),z(-)=',vacplus,-vacplus
c           **** NB ****
c           Combined ion/dipole defect and slab vacuum LMCC suspect
            esion = esionsl
            lmcc_on = .true.
          endif
c
c         End neutral slab dipole branch
c
        elseif( elchrg.ne.zero .and. IABS( ion_opt ) .eq. 3 )then
c
c         Charged slab, with implicit planar Gaussian countercharge in vacuum
c
c         Put counter-charge plane at specific location (stored in second half of rchrg)
c         This will ultimately come from an external source (fluids-DFT, MD, etc)
          do  i=1,3
            zchrglm(i) = rchrg(i)
            zchrgcc(i) = rchrg(3+i)
          enddo
          call SLABCHRG( ion_opt,ndim,
     $     elchrg, origws,zchrglm,zchrgcc,
     $     dipref,rdipref,  esionsl,vacplus,
     $     rprim, n1r,n2r,n3r,  hh,weight, ndbox,nbox,
     $     wk2,     wk3,     wk1 )
c -->      rhlmcc-o vlmcc-o  r-3b
c
           esion = esionsl
           write(IWR,9020) 'Planar charge dble-layer energy=',esionsl
cxxx:`
      write(IWR,*) 'DBG/LMCC Charged plane SLABCHRG not fully vetted'
c
cxxx:start - compute moment of countercharge
      call MOMENT( IWR, elchrg, dipole,rdipole, ndim,
     $ origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $ wk2,    wk4 )
c
      qgrid = zero
      do  i=1,nptr
        qgrid = qgrid + wk2(i)
      enddo
      qgrid = weight*qgrid
      write(IWR,9013) 'DBG/LMCC: rho_lmcc net charge=',qgrid
      call SYMVEC( 1, dipole, vec0, nsym, rmatsym )
      call SYMVEC( 1,rdipole, rsym, nsym, rmatsym )
cxxx:
      write(IWR,9013) 'DBG/LMCC: rho_lmcc dipole=',dipole
      write(IWR,9013) 'DBG/LMCC: rho_lmcc rdipole=',rdipole
cxxx:end
c 
c         Record location of slab planar charge into rchrg()
          do  i=1,3
            rchrg(i) = zchrglm(i)
          enddo
cxxx
          write(IWR,9013)'DBG/LMCC: z_lm(internal)=',zchrglm
          write(IWR,9013)'DBG/LMCC: z_cc(internal)=',zchrgcc
c
          lmcc_on = .true.
c
c         Close double-layer LMCC branch ...
c
        endif
c
c       Close slab (ndim=2) LMCC treatments
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
