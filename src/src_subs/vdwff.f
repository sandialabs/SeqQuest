c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWFF
c
c
      subroutine VDWFF( engyvdw, frcvdw,strvdw, c6vdw,r0vdw,
     $ c6ijvdw,r0ijvdw,
     $ natm,ntyp, nlat  ,nlat2c, itypa, znuc, alamin,
     $ ratm, rlat )
c---------------------------------------------------------------
c Purpose: compute van der Waals correction using pair potentials
c
c Written: Peter A. Schultz, xx-February-2013, for 2.65
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Two forms of pair-potential vdW corrections are currently implemented:
c  (1) DFT-D2, per Stefan Grimme, J. Comput. Chem. 27, 1787 (2006).
c  (2) ULG, per H. Kim, J.-M. Choi, and WAG III, JPhChL 3, 360 (2012).
c
c  Input vdW parameters to this routine are assumed to be Ry and bohr.
c  Output units are also in Ry atomic units
c
c  Note that meanings of d_vdw and r0vdw are different in these two methods.
c   d_vdw is used in different functional forms for damp short-range interations.
c   r0vdw in DFT-D@ is an atomic vdW radius, in ULG it is a dimer bond minimum.
c   rc_vdw os a global long-range cutoff for evaluation of these terms
c  Note that code will use overlap cutoffs to limit range of interacions,
c  as the lattice vectors are set up for this (and <10^-6 terms will not matter)
c
c  This version is non-parallel, all processors do whole thing
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  frcvdw(3,*), strvdw(3,3)
c
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  znuc(ntyp), alamin(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat)
      DIMENSION  c6vdw(ntyp), r0vdw(ntyp)
      DIMENSION  c6ijvdw(ntyp,ntyp), r0ijvdw(ntyp,ntyp)
c
      LOGICAL     DO_VDW
c
c Local declarations:
c
      DIMENSION   rij(3), fij(3)
      LOGICAL     dbgvdw
      DATA        dbgvdw / .false. /
c lstout controls diagnostic output 0-none,1-minimal,n>1-more
      DATA  lstout  / 0 /
      CHARACTER*6  fcnal
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
c Start with a blank slate:
      engyvdw = zero
      call MKZERO( 3*natm, frcvdw )
      call MKZERO( 3*3   , strvdw )
c
c Check if we have a vdW correction to evaluate
      if( .not. DO_VDW() ) RETURN
      call FLGETIWR( IWR )
      call VDWGETNAME( fcnal )
      write(IWR,*) 'VDWFF: van der Waals C6 corrections=',fcnal
c
      call VDWGET_S( sfac )
      s_vdw = sfac
      call VDWGET_RC( rcfac )
      rc_vdw = rcfac
c
      if( dbgvdw )then
        call VDWGETNAME( fcnal )
        call VDWGET_D( dfac )
        ivdw = IVDWTYPE()
 9020   format(1x,a,i3,1x,a6, 3f12.6)
        write(IWR,9020) 'VDW-FF: i/vdw.s.d.rc=',
     $   ivdw       , fcnal, sfac,dfac,rcfac
 9021   format(2x,a,2i6)
 9022   format(2x,a,(9f18.6))
        write(IWR,9022) 'VDWFF: c6vdw()=',(c6vdw(ityp),ityp=1,ntyp)
        write(IWR,9022) 'VDWFF: r0dww()=',(r0vdw(ityp),ityp=1,ntyp)
      endif
c
c Set up the local cutoff criterion
      argmx = half*cutgaus
      emax = cutslow
c
c Loop over atom pairs, collect terms
c
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        if(dbgvdw) write(IWR,9021) 'jatm,jtyp=',jatm,jtyp
  
        if( znuc(jtyp).eq.zero ) goto 2000
c
        c6j = c6vdw(jtyp)
        r0j = r0vdw(jtyp)
        if(dbgvdw) write(IWR,9022) 'cvdw(atom): c6j,r0j=',c6j,r0j
        if( c6j .le. zero ) goto 2000
c
        alaminj = alamin(jtyp)
c
        do 1000 iatm=jatm,natm
          ityp = itypa(iatm)
          if(dbgvdw) write(IWR,9021) '  iatm,ityp=',iatm,ityp
          if( znuc(ityp).eq.zero ) goto 1000
c
          c6i = c6vdw(ityp)
          r0i = r0vdw(ityp)
          if( c6i .le. zero ) goto 1000
c
c         Construct cross-potentials, as needed
C          call VDWIJMIX( c6i,c6j, c6ij, r0i,r0j,r0ij )
          c6ij = c6ijvdw(ityp,jtyp)
          r0ij = r0ijvdw(ityp,jtyp)
          if(dbgvdw) write(IWR,9022)'  cvdw(pair): c6ij,r0ij=',c6ij,r0ij
c
c         Check that this pairwise interaction is allowed
          if( c6ij .le. zero ) goto 1000
c         Determine maximum range of interactions, using orbital overlap cutoffs
          alamini = alamin(ityp)
          rijmax = SQRT( emax*(alamini+alaminj) / (alamini*alaminj) )
          rcutoff = MIN( rijmax, rc_vdw )
c
          do 900 ilat=1,nlat
c
            rij2 = zero
            do  ir=1,3
              rij(ir) = ratm(ir,jatm) - (ratm(ir,iatm)+rlat(ir,ilat))
              rij2 = rij2 + rij(ir)**2
            enddo
            rijmag = SQRT( rij2 )
c           This ensures atom does not see itself, and cuts off long-range vdW
            if( rijmag .lt. 0.1d0 .or. rijmag .gt. rcutoff ) goto 900
c
c           Compute vdw term
            call VDWPAIR( e_vdw,de_vdw, rijmag,c6ij,r0ij )
C           write(IWR,914) '    pair vdw: ilat,rij,e_vdw=',ilat,rijmag,e_vdw
C  914      format(1x,a,i5,f9.4,3x,1pd15.8)
c
            if( iatm .eq. jatm )then
              e_vdw = e_vdw / two
              de_vdw = de_vdw / two
            endif
c
c           Add contribution to energy force and stress
            engyvdw = engyvdw + e_vdw
c
            do  ir=1,3
              fij(ir) = de_vdw*(rij(ir)/rijmag)
            enddo
            do  ir=1,3
              frcvdw(ir,jatm) = frcvdw(ir,jatm) - fij(ir)
              frcvdw(ir,iatm) = frcvdw(ir,iatm) + fij(ir)
            enddo
c
            do  jr=1,3
              do  ir=1,3
                strvdw(ir,jr) = strvdw(ir,jr) + fij(ir)*rij(jr)
              enddo
            enddo
c
c           Close Atom i + latvec loops:
  900     continue
c
 1000   continue
c
c       Close Atom j loop:
 2000 continue
c
c Scale the computed vdw results by the (functional-dependent) global scale factor
c
      if( dbgvdw ) write(IWR,'(a,f20.8)')
     $ 'End vdwff: Unscaled engyvdw, sfactor=', engyvdw,s_vdw
c
      engyvdw = engyvdw*s_vdw
      call DSCAL( 3*natm, s_vdw, frcvdw,1 )
      call DSCAL( 3*3   , s_vdw, strvdw,1 )
c
c    That's all Folks!
c
      RETURN
      END
