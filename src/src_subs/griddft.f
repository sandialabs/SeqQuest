c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDDFT
c
c
      subroutine GRIDDFT( IWR, igrd0fl,igrhofl,
     $ i1cgrfl, anyxc1c,
     $ anycore, do_gga,do_spin, iter,lvlout,
     $ nptr, weight,xcrhocut, exctot,excslo,corslow,eharris,
     $ rho0,   rhoc,    rhoa,  vxca,  exca,  exc0,
     $                  rhob,  vxcb,  excb,  wksp1, wksp2 )
c---------------------------------------------------------------
c Purpose: evaluate xc-potentials on grid for scf
c
c Written: Peter A. Schultz, 19-July-2001, for 2.49 (spin)
c
c Revision history:
c  29Apr14-PAS/2.66: master-only igrd0fl
c  17Spt08-RJM/2.62:: modified variable passing both spins to AM05SP
c  28Jan08-PAS/2.60e: patch igrhofl rewind to master-only
c  26Jun03-PAS/2.56: pass g2rho[up] into spin-gga corr (blypsp)
c   7Apr03-PAS/2.55: fix rhotot(up,dn) output; extract core
c  19Jul01-PAS/2.49: pushed down from main program
c---------------------------------------------------------------
c
c Notes:
c  The big grid arrays are used multiply as scratch.  The names
c  they are given reflect a compromise between (1) their contents
c  on input to, or output from, this routine, and (2) how they
c  are used in the trickier parts of this routine.  The intent
c  is to be as compact as possible, and still be understandable.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, half = 0.5d0, one = 1.d0 )
      PARAMETER  ( two = 2.d0 )
c
      LOGICAL  anyxc1c
      LOGICAL  anycore, do_gga,do_spin
c input arrays ...
c  ... reference density
      DIMENSION  rho0(nptr)
c input/output ...
c  ... delta density, up/dn
      DIMENSION  rhoa(nptr), rhob(*)
c output arrays ...
c  ... xc-potentials for up and dn
      DIMENSION  vxca(nptr),exca(nptr),exc0(nptr)
      DIMENSION  vxcb(*),   excb(*)
c scratch space ...
c  ... core density
      DIMENSION  rhoc(*)
      DIMENSION  wksp1(*), wksp2(*)
c
c >>>> EXECUTABLE CODE
c
 9020 format(1x,a,3f20.10)
c
c How much of core/ref density to add?
      spinfac = one
      if( do_spin ) spinfac = half
c
c Add rho-zero reference density into delta (up/dn) density:
      call DAXPY( nptr, spinfac, rho0,1, rhoa,1 )
      if( do_spin ) call DAXPY( nptr, spinfac, rho0,1, rhob,1 )
c
c We are assumeing that rho0 has been read off igrd0fl before here ...
      cortot = zero
      if( anycore )then
c       Recover core density from igrd0fl ...
        call MPBCREADBIG( igrd0fl, nptr, rhoc )
        cortot = weight*DDOT( nptr, one,0, rhoc,1 )
c        ... and add into (up/dn) density for xc-calc
        call DAXPY( nptr, spinfac, rhoc,1, rhoa,1 )
        if( do_spin ) call DAXPY( nptr, spinfac, rhoc,1, rhob,1 )
      endif
c
      if( do_gga )then
c       Retrieve grho, g2rho, and gpggp (gga gradients) from igrhofl:
        call MPBCREWIND( igrhofl )
        call MPBCREADBIG( igrhofl, nptr, vxca )
        call MPBCREADBIG( igrhofl, nptr, exca )
        call MPBCREADBIG( igrhofl, nptr, exc0 )
        if( do_spin )then
c         And retrieve grho[up] and grho[dn] from igrhofl:
          call MPBCREADBIG( igrhofl, nptr, wksp1 )
          call MPBCREADBIG( igrhofl, nptr, excb )
        endif
      endif
c
c Compute xc-functionals on grid:
c
      if( do_spin )then
c
c       Spin-polarized
c
        if( do_gga )then
c
c         Retrieve g2rho-up from igrhofl ...
          call MPBCREADBIG( igrhofl, nptr, wksp2 )
c          ... and put it into vxcb for GRIDCGGA:
          call DCOPY( nptr, wksp2,1, vxcb,1 )
c
c         Compute correlation part of GGA functional ...
          call GRIDCGGA( xcrhocut, weight,nptr,
     $     rhuptot, vcuptot,ecuptot,  rhdntot, vcdntot,ecdntot,
     $     rhoa, vxca,  exca,  exc0,
c -->(i)   rh-up grho   g2rho  gpggp
c <--(o)   rh-up vxc-up exc-up
     $     rhob, vxcb,   excb,   wksp1 )
c -->(i)   rh-dn g2-up   grho-dn grho-up
c <--(o)   rh-dn vxc-dn  exc-dn  grho-up
c
c          ... retrieve gpggp-up from igrhofl:
c             (note grho-up survived from gridcgga, in wksp1,
c              and g2rho-up continued in wksp2)
          call MPBCREADBIG( igrhofl, nptr, exc0 )
c
c          ... compute up-spin GGA exchange:
          call GRIDXGGA( xcrhocut, weight,nptr,
     $     rhuxtot, vxuptot, exuptot,
     $     rhoa,   rhob, wksp1,  wksp2,   exc0,     vxca,    exca )
c -->      rhup-i  grhup-i g2rhup-i gpggpup-i vxcup-io excup-io
c
c          ... retrieve grho-dn, g2rho-dn, and gpggp-dn from igrhofl:
          call MPBCREADBIG( igrhofl, nptr, wksp1 )
          call MPBCREADBIG( igrhofl, nptr, wksp2 )
          call MPBCREADBIG( igrhofl, nptr, exc0 )
c
c          ... compute dn-spin GGA exchange:
          call GRIDXGGA( xcrhocut, weight,nptr,
     $     rhdxtot, vxdntot, exdntot,
     $     rhob,   rhoa, wksp1,  wksp2,   exc0,     vxcb,    excb )
c -->      rhdn-i  grhdn-i g2rhdn-i gpggpdn-i vxcdn-io excdn-io
c
          vctot = vcuptot + vcdntot
          ectot = ecuptot + ecdntot
          vxtot = vxuptot + vxdntot
          extot = exuptot + exdntot
c
          if( lvlout.gt.3 )then
            write(IWR,9020) 'gga up-spin exch energy, vx:rho=',
     $       exuptot,vxuptot
            write(IWR,9020) 'gga dn-spin exch energy, vx:rho=',
     $       exdntot,vxdntot
            write(IWR,9020) 'gga up-spin corr energy, vc:rho=',
     $       ecuptot,vcuptot
            write(IWR,9020) 'gga dn-spin corr energy, vc:rho=',
     $       ecdntot,vcdntot
          endif
c
          vxcuptot = vxuptot + vcuptot
          excuptot = exuptot + ecuptot
          vxcdntot = vxdntot + vcdntot
          excdntot = exdntot + ecdntot
c
        else
c
          call GRIDSLDA( xcrhocut, weight,nptr,
     $     rhuptot, vxcuptot,excuptot,  rhdntot, vxcdntot,excdntot,
     $     rhoa,   vxca,    exca,     rhob,   vxcb,  excb )
c -->      rh-up/i vxc-up/o exc-up/o  rh-dn/i vxc-dn/o exc-dn/o
c
        endif
c
        rhotot = rhuptot + rhdntot
        vxctot = vxcuptot + vxcdntot
        exctot = excuptot + excdntot
c
      else
c
c       Spin-restricted
c
        if( do_gga )then
c
          call GRIDGGA( xcrhocut, weight,nptr,
     $     rhotot, vxctot,exctot,
     $     rhoa, vxca,  exca,   exc0 )
c -->(i)   rho-i grho-i g2rho-i gpggp-i
c <--(o)         vxc-o  exc-o
c
        else
c
          call GRIDLDA( xcrhocut, weight,nptr,
     $     rhotot, vxctot,exctot,
     $     rhoa, vxca,  exca )
c -->      rho-i vxc-o  exc-o
c
        endif
c
        rhuptot = rhotot / two
        rhdntot = rhuptot
c
      endif
c
      rhotot = rhotot - cortot
      rhuptot = rhuptot - cortot/two
      rhdntot = rhdntot - cortot/two
      eharris = exctot - vxctot
c
c  Retrieve Vxc0 from igrd0fl, and remove from Vxc to get del-Vxc ...
      call MPBCREADBIG( igrd0fl, nptr, exc0 )
      call DAXPY( nptr, -one, exc0,1, vxca,1 )
      if( do_spin ) call DAXPY( nptr, -one, exc0,1, vxcb,1 )
c      ... and need core:Vxc0 integral for Harris functional
      if( anycore ) vcorfst = weight*DDOT( nptr, rhoc,1, exc0,1 )
c
c  Retrieve Exc0 from igrd0fl
      call MPBCREADBIG( igrd0fl, nptr, exc0 )
      call DAXPY( nptr, -one, exc0,1, exca,1 )
      exc0rho = weight*DDOT( nptr, rhoa,1, exc0,1 )
      if( do_spin )then
        call DAXPY( nptr, -one, exc0,1, excb,1 )
        exc0rhb = weight*DDOT( nptr, rhob,1, exc0,1 )
        exc0rho = exc0rho + exc0rhb
      endif
      excslo = exctot - exc0rho
c
c  Reconstitute slow-density by taking rho0(+rhoc) out of total:
      if( anycore ) call DAXPY( nptr, one, rhoc,1, rho0,1 )
      call DAXPY( nptr, -spinfac, rho0,1, rhoa,1 )
      if( do_spin ) call DAXPY( nptr, -spinfac, rho0,1, rhob,1 )
c  rho0() is considered destroyed at this point
c
      if( anyxc1c )then
c       Non-standard/non-working code; use rho0() as scratch
c       Retrieve dynamic fast xc grid fields from i1cgrfl ...
        REWIND( i1cgrfl )
c        ... get -vxc1c
        call READBIG( i1cgrfl, nptr, rho0 )
c        ... and remove from grid-vxc to get slow-vxc
        call DAXPY( nptr, one, rho0,1, vxca,1 )
c        ... get -exc1c
        call READBIG( i1cgrfl, nptr, rho0 )
c        ... and remove from grid-exc to get slow-exc
        call DAXPY( nptr, one, rho0,1, exca,1 )
c        ... and add to Exc0 to get "fast" Exc:
        call DAXPY( nptr, -one, rho0,1, exc0,1 )
      endif
c
      if( lvlout.gt.0 )then
        write(IWR,*) ' '
        if( do_spin )then
          write(IWR,9020) 'grid-integrated up-density total =',rhuptot
          write(IWR,9020) 'grid-integrated dn-density total =',rhdntot
        endif
        write(IWR,9020) 'grid-integrated valence-density total =',rhotot
      endif
      if( lvlout.gt.2 )then
        if( do_spin )then
          write(IWR,9020) 'grid-integrated up-spin xc energy=',excuptot
          write(IWR,9020) 'grid-integrated dn-spin xc energy=',excdntot
        endif
        write(IWR,9020) 'grid-integrated total xc-energy =     ', exctot
        write(IWR,9020) 'grid-integrated total slow xc-energy =', excslo
        if( do_spin )then
          write(IWR,9020) 'grid-integrated up-spin vxc energy=',vxcuptot
          write(IWR,9020) 'grid-integrated dn-spin vxc energy=',vxcdntot
        endif
        write(IWR,9020) 'grid-integrated total vxc:rho energy=',vxctot
      endif
      if( anycore )then
        corfast = weight*DDOT( nptr, rhoc,1, exc0,1 )
c
        corslow = spinfac*weight*DDOT( nptr, rhoc,1, exca,1 )
c       Need to pull core:Vxc intergral out of Harris term
        vxccor = vcorfst + spinfac*weight*DDOT( nptr, rhoc,1, vxca,1 )
        if( do_spin )then
          corslob = spinfac*weight*DDOT( nptr, rhoc,1, excb,1 )
          vxccorb = spinfac*weight*DDOT( nptr, rhoc,1, vxcb,1 )
          corslow = corslow + corslob
          vxccor = vxccor + vxccorb
        endif
        coregrd = corfast + corslow
        if( lvlout.gt.1 )then
          write(IWR,9020)     'grid-int core charge  =', cortot
          if( lvlout.gt.2 )then
            write(IWR,9020)   'grid-int core exc total =', coregrd
            if( lvlout.gt.3 )then
              write(IWR,9020) 'grid-int core fast exc  =', corfast
              write(IWR,9020) 'grid-int core slow exc  =', corslow
              write(IWR,9020) 'grid-int core vxc total =', vxccor
            endif
          endif
        endif
        eharris = eharris + vxccor
      endif
c
      if( iter.eq.0 .and.lvlout.gt.1 )
     $ write(IWR,9020) 'Harris Fcnal XC piece=',eharris
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End
