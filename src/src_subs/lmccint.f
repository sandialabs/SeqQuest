c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LMCCINT
c
c
      subroutine LMCCINT( IWR, igrdref,  lvlout, lmcc_on,
     $ esnsgr, esion,eslmref,  ndim, nptr, weight,
     $ rhoslo, vesslo, rholmcc, veslmcc, delrho )
c---------------------------------------------------------------
c Purpose: evaluate integrals needed to compute contributions
c          to LMCC energy.
c
c Written: Peter A. Schultz, 30-May-2001 for v2.47
c
c Revision history:
c  22May11-PAS/2.62: cosmetics, so far
c  10Sep05-PAS/2.59: chemical potential perspective diagnostics
c  30May01-PAS/2.47: pushed down from main program
c---------------------------------------------------------------
c
c Notes:
c  Implements the method for treating supercell moments described
c  in:  P.A. Schultz,  PRB 60, 1551 (1999);  PRL 84, 1942 (2000).
c
c on input
c   rhoslo()  = fft density
c   vesslo()  = fft es potential
c   rholmcc() = lmcc density
c   veslmcc() = lmcc es potential
c   delrho()  = delta-rho [rho - rho0(ref atoms)]
c on output
c   rhoslo()  = full slow density (fft+lmcc)
c   vesslo()  = full slow es potential (fft+lmcc)
c This routine does not care whether grid fields come boxed or
c not, as long as all the fields are the same way.
c
c Want to evaluate convolution INT[delrho:slorho] where
c   delrho = del = fast + fft + ion
c   slorho = slo = fft + ion
c and want to evaluate all ion:X terms using X density and
c ion potentials, and use a special evaluation for the ion:ion
c term in case it is provided.  The way I do this is:
c (density left, potential right)
c
c    del:slo = (fast+fft+ion):(fft+ion) "esdelslo"
c        - ion:fft + fft:ion - ion:ion(grid) + ion:ion(spec)
c        "epionvslo" "evionpslo"  "esiongr"     "esion"
c        + ion:fft[ref] [reference system lmcc integrals ell]
c     Note: esion was computed previously in "ionpot" or "dippot"
c I also would like to compute the "chemical potential" shift,
c which is the LMCC grounding of the energy to crystal reference
c electron reservoir:
c    esulmcc = - rho_lm:ves_def + v_lm:rho_def
c              + rho_lm:ves_ref - v_lm:rho_ref
c The reference terms were done earlier, and sum is passed
c in "eslmref", while the defect terms are evaluated here.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, half=0.5d0, one=1.d0 )
c
c Input arrays:
      LOGICAL    lmcc_on
c
c Input/Output arrays:
      DIMENSION  rhoslo(nptr),vesslo(nptr),  delrho(nptr)
c Input arrays:
      DIMENSION  rholmcc(nptr),veslmcc(nptr)
c
c Local declarations:
      DATA  toev / 13.605805d0 /
cxxx      DATA  lvldbg / 0 /
      DATA  lvldbg / 0 /
c
c >>>> EXECUTABLE CODE:
c
 9020 format(1x,a,3f20.10)
      lvllst = MAX( lvlout, lvldbg )
c
      epionvslo = half*weight*DDOT( nptr, rholmcc,1, vesslo ,1 )
      evionpslo = half*weight*DDOT( nptr, rhoslo ,1, veslmcc,1 )
      eslmdef = -epionvslo + evionpslo
      esulmcc = eslmdef + eslmref
      esiongr  = half*weight*DDOT( nptr, rholmcc,1, veslmcc,1 )
c
c  Combine lmcc into slow densities, and potentials:
      call DAXPY( nptr, one, rholmcc,1, rhoslo,1 )
      call DAXPY( nptr, one, veslmcc,1, vesslo,1 )
c
c  Take the last integral I need ...
      esdelslo = half*weight*DDOT( nptr, delrho ,1, vesslo ,1 )
c
c  And put it all together ...
      esnsgr = esdelslo + esulmcc - esiongr + esion
c
      if( lvllst .gt. 3 )then
        write(IWR,'(1x,a)') 'LMCCINT energy breakdown:'
        write(IWR,9020) '  epionvslo(pivd)  = ',-epionvslo
        write(IWR,9020) '  evionpslo(vipd)  = ', evionpslo
        write(IWR,9020) '  esdelslo(pdvd)   = ', esdelslo
        write(IWR,9020) '  grid ion self-energy esiongr   =',-esiongr
        write(IWR,9020) '  grid ion self-energy(analytic)`=', esion
        write(IWR,9020) '  total grid es energy=',esnsgr
      endif
      eslmdev = toev*eslmdef
      if( lvlout.gt.2 )then
        write(IWR,9020) 'LMCC defect int(Ry): pivd(pbc)=',
     $   epionvslo
        write(IWR,9020) 'LMCC defect int(Ry): vipd(loc)=',
     $   evionpslo
        write(IWR,9020) 'LMCC total defect int(Ry,eV): -pivd+vipd=',
     $   eslmdef,eslmdev
      endif
      if( ndim.gt.0 )then
        if( igrdref.gt.0 )then
          if( lvlout.gt.1 )then
            write(IWR,9020) 'LMCC total reference int(Ry)=', eslmref
          endif
        else
          write(IWR,*) '***** Warning: No reference term in LMCC'
        endif
      endif
      esulmev = toev*esulmcc
      if( lvlout.gt.0 )then
        write(IWR,9020) 'LMCC chemical potential shift(Ry,eV)=',
     $   esulmcc,esulmev
      endif
c
c    That's all Folks!
c
      RETURN
      END
