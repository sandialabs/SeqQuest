c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GGASTRESS
c
c
      subroutine GGASTRESS( strgga, do_gga, do_spin,
     $ idrhfl0,idrhofl,idrhsfl,igrhofl,
     $ ndim, nptr,weight,xcrhocut,
     $ rhoa, rhob, g_dxc, g2_dxc2, gg_d2dn,
     $ drhoa, drhob )
c---------------------------------------------------------------
c Purpose: compute the extra GGA-specific terms needed in stress
c
c Written: Peter A. Schultz, 5-April-2002, for v2.53
c
c Revision history:
c  17Spt08-RJM/2.62: added opposite spin passing to GRIDXSTR needed for am05sp 
c  27Jun07-PAS/2.60: task-parallel - master-only
c   9Oct03-AEM+PAS/2.57: spin-dependent stress added for GGA
c  11Jun03-PAS/2.55: add Laplacian stress term (BLYP)
c   5Apr02-PAS/2.53: modeled off of "GRIDDFT"
c---------------------------------------------------------------
c
c Notes:
c   I need to acquire the extra term in the stress formula that
c   appears because of the gradient-rho that appears in the GGA
c   fcnal [See: A. Dal Corso and R. Resta, PRB 50, 4327 (1994)].
c
c   Let drh[i] = the i'th component of the density gradient:
c      drh[i] = d[rho]/dx[i]
c   Then need the integral over the unit cell:
c     stress(i,j) = INT( ( d[rho*exc]/d[drh[i]] ) * drh[j] )
c   I compute the potential in such a way that I get
c     dxc(r) = d[rho*exc]/d[drh[i]] / drh[i]
c   so that I need only to complete the integral as
c     stress(i,j) = INT( dxc(r) * drh[i](r) * drh[j](r) )
c   Hence, I need to compute dxc(r) in the functional, and sum
c   over the grid the product with the i'th and j'th component
c   of the density gradient.
c
c   For the Laplacian term need something analogous:
c     strlap(i,j) = INT( dxc2(r) * d2_rh(r)/d_x_i*d_x_j )
c
c   This calculation is moot in LDA.
c
c   Spin polarized GGA adds two spin-resolved forms of 1st derivative
c   integrals, and splits the Laplacian second derivative term
c   into two spin resolved pieces (thx to Ann Mattsson).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, half=0.5d0, one=1.d0 )
c
c Output:
      DIMENSION  strgga(3,3)
c Input::
      LOGICAL    do_gga, do_spin
c Input/scratch: total xc-density (if spin, spin-up and spin-dn)
      DIMENSION  rhoa(nptr), rhob(nptr)
c Scratch:
      DIMENSION  g_dxc(nptr),g2_dxc2(nptr),gg_d2dn(nptr)
      DIMENSION  drhoa(nptr,3), drhob(nptr,3)
c
c Local:
      DIMENSION  strdxc(3,3), strlap(3,3)
      DIMENSION  strsav(3,3),strdel(3,3)
      LOGICAL    do_laplac
      DATA       lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( lstout.gt.0 ) call FLGETIWR( IWR )
      if( IWR.le.0 ) lstout = 0
 9013 format(1x,a,3f16.8)
c
      call MKZERO( 9, strgga )
      call MKZERO( 9, strdxc )
      call MKZERO( 9, strlap )
c
c Get MP info
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
c Only the master will do this work
      if( iproc .ne. master ) goto 990
c
      if( .not. do_gga .or. ndim .eq. 0 ) RETURN
c
c     Note if need concern ourselves with Laplacian in the energy fcnal
      call DFTGETLAP( do_laplac )
c
      if( do_gga )then
c       Retrieve grho, g2rho, and gpggp (gga gradients) from igrhofl
        REWIND( igrhofl )
        call READBIG( igrhofl, nptr, g_dxc )
        call READBIG( igrhofl, nptr, g2_dxc2 )
        call READBIG( igrhofl, nptr, gg_d2dn )
      endif
c
      if( do_spin )then
c       Spin-polarized GGA stress * * * * * * * * * * * * * * * * * * *
c
c       Retrieve grho-up, grho-*dn*, g2rho-up from igrhofl
        call READBIG( igrhofl, nptr, drhoa(1,1) )
        call READBIG( igrhofl, nptr, drhoa(1,3) )
        call READBIG( igrhofl, nptr, drhoa(1,2) )
c
c       Compute correlation stress fcnal on grid
        call GRIDCSTR( xcrhocut, nptr, weight,
     $   rhoa,  rhob,  g_dxc,   g2_dxc2,   gg_d2dn,
c -->    rhoa-i rhob-i g_dxc-io g2_d2up-io gg_d2dn-io
     $   drhoa(1,1),drhoa(1,3),drhoa(1,2), drhob(1,1),drhob(1,2) )
c -->    grhoa-i    grhob-i    g2rhoa-i    dxcup-o    dxcdn-o
c
c       Retrieve gpggp-up from igrhofl; grho-up,g2rho-up survived above
        call READBIG( igrhofl, nptr, drhoa(1,3) )
c
c       Compute spin-up exchange stress fcnal on grid
        call GRIDXSTR( xcrhocut, nptr, weight,
     $   rhoa,  rhob,drhoa(1,1),drhoa(1,2),drhoa(1,3),
c -->    rho-i  grho-i     g2rho-i    ggrho-i
     $   drhob(1,1) )
c -->    dxc-io[up]
c
c       Retrieve grho-dn, g2rho-dn, and gpggp-dn from igrhofl
        call READBIG( igrhofl, nptr, drhoa(1,1) )
        call READBIG( igrhofl, nptr, drhoa(1,2) )
        call READBIG( igrhofl, nptr, drhoa(1,3) )
c
c       Compute spin-dn exchange stress fcnal on grid
        call GRIDXSTR( xcrhocut, nptr, weight,
     $   rhob,  rhoa,drhoa(1,1),drhoa(1,2),drhoa(1,3),
c -->    rho-i  grho-i     g2rho-i    ggrho-i
     $   drhob(1,2) )
c -->    dxc-io[dn]
c
c       Have stress functionals, now the fun part of taking integrals
c
c       Move some stress fcnals (dxc-up/dn) to more convenient locations
        call DCOPY( nptr, drhob(1,1),1, rhoa,1 )
        call DCOPY( nptr, drhob(1,2),1, rhob,1 )
c
c       Position the idrh grid density derivatives files
        REWIND( idrhfl0 )
        REWIND( idrhofl )
        REWIND( idrhsfl )
c
c       Construct 1st derivatives of total (up+dn) xc-density
c
c       Retrieve 1st derivatives of reference+core density from idrhfl0
        do  iu=1,ndim
          call READBIG( idrhfl0, nptr, drhoa(1,iu) )
        enddo
c       Retrieve 1st derivatives of total scf delta-density from idrhofl
        do  iu=1,ndim
          call READBIG( idrhofl, nptr, drhob(1,iu) )
        enddo
c       Add d1drho to d1rho0 to get full (up+dn) xc-density derivatives
        call DAXPY( ndim*nptr, one, drhob,1, drhoa,1 )
c
        call DCOPY( 9, strdxc,1, strsav,1 )
c
c       Integrate stress contribution involving total 1st derivatives
c
        call DFTSTR( ndim, strdxc, nptr,weight,
     $   g_dxc, drhoa )
c -->    dxc-i  drho-3i
c
        if( lstout.gt.2 )then
          call DCOPY( 9, strdxc,1, strdel,1 )
          call DAXPY( 9, -1.0d0, strsav,1, strdel,1 )
          write(IWR,9013) 'GGASTR: dftstr-1rh =',(strdel(id,id),id=1,3)
          call DCOPY( 9, strdxc,1, strsav,1 )
        endif
c
c       Construct 1st derivatives of total up/dn-spin xc-densities
c
c       Remove del-density derivs to recover reference xc-density derivs
        call DAXPY( ndim*nptr, -one, drhob,1, drhoa,1 )
c       Now some tight-packed manipulation to get up/dn at same time
        do  iu=1,ndim
c         Make copy of d1rho0[iu]
          call DCOPY( nptr, drhoa(1,iu),1, g_dxc,1 )
c         Retrieve d1drho-up[iu] from idrhsfl
          call READBIG( idrhsfl, nptr, drhoa(1,iu) )
c         Make d1drho-dn[iu] = d1drho[iu] - d1drho-up[iu]
          call DAXPY( nptr, -one, drhoa(1,iu),1, drhob(1,iu),1 )
c         Add half of d1rho0 to d1drho-up/dn to get total up/dn-derivs
          call DAXPY( nptr, half, g_dxc,1, drhoa(1,iu),1 )
          call DAXPY( nptr, half, g_dxc,1, drhob(1,iu),1 )
        enddo
c
c       Integrate stress contribution involving up-spin 1st derivatives
c
        call DFTSTR( ndim, strdxc, nptr,weight,
     $   rhoa, drhoa )
c -->    dxc-i drho-3i
c
        if( lstout.gt.2 )then
          call DCOPY( 9, strdxc,1, strdel,1 )
          call DAXPY( 9, -1.0d0, strsav,1, strdel,1 )
          write(IWR,9013) 'GGASTR: dftstr-1up =',(strdel(id,id),id=1,3)
          call DCOPY( 9, strdxc,1, strsav,1 )
        endif
c
c       Integrate stress contribution involving dn-spin 1st derivatives
c
        call DFTSTR( ndim, strdxc, nptr,weight,
     $   rhob, drhob )
c -->    dxc-i drho-3i
c
        if( lstout.gt.2 )then
          call DCOPY( 9, strdxc,1, strdel,1 )
          call DAXPY( 9, -1.0d0, strsav,1, strdel,1 )
          write(IWR,9013) 'GGASTR: dftstr-1dn =',(strdel(id,id),id=1,3)
          call DCOPY( 9, strdxc,1, strsav,1 )
        endif
c
        if( do_laplac )then
c         Integrate stress contribution involving up/dn 2nd derivatives
c
c         Position idrhfl0,idrhofl,idrhsfl to beginning of 2nd derivs
          if( ndim .lt. 3 )then
c           Step through remaining 1st derivatives in idrh files
            do  id=ndim+1,3
              call READSKP( idrhfl0, nptr )
              call READSKP( idrhofl, nptr )
              call READSKP( idrhsfl, nptr )
            enddo
          endif
c
c         Following assumes idrh files positioned to read 2nd derivs
          call DFTSTR2SP( idrhfl0,idrhofl,idrhsfl,
     $     ndim, strlap, nptr,weight,
     $     g2_dxc2, gg_d2dn,  drhoa(1,1),drhoa(1,2),drhoa(1,3) )
c -->      dxc2up-i dxc2dn-i  ddrho0-s   ddrhoa-s   ddrhob-s
c
          if( lstout.gt.1 )
     $    write(IWR,9013) 'GGASTR: dftstr-2tot=',(strlap(id,id),id=1,3)
        endif
c
        if( lstout.gt.1 )then
          write(IWR,9013) 'GGASTR: dftstr-1tot=',(strdxc(id,id),id=1,3)
        endif
c
c       End spin-GGA stress * * * * * * * * * * * * * * * * * * * * * *
      else
c       Non-spin GGA stress * * * * * * * * * * * * * * * * * * * * * *
c
c       Compute stress functionals dxc(r), dxc2(r) on grid
        call GRIDSTR( xcrhocut, nptr, weight,
     $   rhoa,   g_dxc,    g2_dxc2,   gg_d2dn )
c-->     rho-io  g_dxc-io  g2_dxc2-io gpggp-i
c
c       Position density derivative files
        REWIND( unit=idrhfl0 )
        REWIND( unit=idrhofl )
c
c       Retrieve 1st derivatives of reference+core density from idrhfl0
        do  iu=1,ndim
          call READBIG( idrhfl0, nptr, drhoa(1,iu) )
        enddo
c
c       Retrieve 1st derivatives of scf delta-density from idrhofl
        do  iu=1,ndim
          call READBIG( idrhofl, nptr, drhob(1,1) )
c          ... and combine to get total density derivatives
          do  nr=1,nptr
            drhoa(nr,iu) = drhoa(nr,iu) + drhob(nr,1)
          enddo
        enddo
c
c       Integrate stress contribution involving 1st derivatives
        call DFTSTR( ndim, strdxc, nptr,weight,
     $   g_dxc, drhoa )
c -->    dxc-i  drho-3i
c
        if( do_laplac  )then
c         Integrate stress contribution involving 2nd derivatives
c
c         Position idrhfl0,idrhofl to beginning of 2nd derivs
          if( ndim .lt. 3 )then
c           Step through remaining 1st derivatives in idrh files
            do  id=ndim+1,3
              call READSKP( idrhfl0, nptr )
              call READSKP( idrhofl, nptr )
            enddo
          endif
c
c         Following assumes idrh files positioned to read 2nd derivs
          call DFTSTR2( idrhfl0,idrhofl,idrhsfl,
     $     ndim, strlap, nptr,weight,
     $     g2_dxc2, drhoa(1,1),drhoa(1,2) )
c -->      dxc2-i   ddrho0-s   ddrho-s
c
        endif
c
        if( lstout.gt.1 )then
          write(IWR,9013) 'GGASTR: dftstr-2tot=',(strlap(id,id),id=1,3)
          write(IWR,9013) 'GGASTR: dftstr-1tot=',(strdxc(id,id),id=1,3)
        endif
c
c       End non-spin GGA stress * * * * * * * * * * * * * * * * * * * *
      endif
c
c Finally, load grid integral results into the output stress:
c (note minus sign, given by Eq.24 of Dal Corso and Resta)
c
      do  jd=1,ndim
        do  id=1,ndim
          strgga(id,jd) = strgga(id,jd) - strdxc(id,jd)
     $                  + strlap(id,jd)
        enddo
      enddo
c
c If non-master, drop to here:
  990 continue
c
      if( nprocs .gt. 1 )then
        lenmsg = 9
        call MPBCAST8( master, lenmsg, strgga, icomm )
      endif
      if( lstout.gt.0 )then
        write(IWR,9013) 'GGASTR: ggastress=',(strgga(id,id),id=1,3)
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
