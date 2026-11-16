c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CELLOPT
c
c
      subroutine CELLOPT( iucstat, ihistfl,
     $ do_cell, ndim, iuctype, icellscheme,strsbroy,
     $ iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $ strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $ strlast, engytotl,
     $ rprim, orig,rscale, strtot, nsym,rmatsym,
     $ nk,veck, dolvatm, natm,ratm,ratm0 )
c---------------------------------------------------------------
c Purpose: initialize user-adjustable parameters for cell opt
c
c Written: Rick P. Muller and Peter A. Schultz, 1-May-2002, v2.53
c
c Revision history:
c   1Oct13-PAS/2.65a: patch dscal args, (b)5Feb14-ACP patch another dscal arg
c   5Apr13-PAS/2.65: axial constraints; independent cell factors
c  28Jun12-PAS?2.63: fix steepest (had damped when hole accelerates)
c  23Sep08-PAS/2.62: output cell volume
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c   8Aug03-PAS/2.56: modify k-vector handling
c  20May03-PAS/2.55: history file added
c  10Jan03-PAS/2.55: limited-history Broyden
c  23Jun02-PAS/2.53e: oops - forgot to symmetrize stress tensor!
c  19May02-PAS/2.53b: 1D and 2D cell opt code made to work
c   1May02-PAS: Rick Muller code pushed down from main program
c---------------------------------------------------------------
c
c The input scf stress "strtot" is assumed to be in Rydberg.
c
c The user-specified convergence criterion "cellconv" and the
c external stress "strxtrnl" - and the various quantities
c determining the blend - use units of GPa (for bulk, and the
c appropriate units generated in "str2gpa" and shown in the
c string "strlbl" for other dimensions).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, half=0.5d0, one=1.d0, two=2.d0 )
      PARAMETER  ( three=3.d0 )
c
c Input:
      LOGICAL    do_cell
c      ... scf stress (in Ry), and external stress (in GPa)
      DIMENSION  strtot(3,3), strxtrnl(3,3)
      DIMENSION  strlast(3,3), ucstepfac(3,3)
      DIMENSION  rmatsym(*)
c Output:
      DIMENSION  rprim(3,3)
      DIMENSION  orig(3),rscale(3)
c Input/Output coordinate arrays modified with lattice vectors:
      DIMENSION  veck(3,*)
      LOGICAL    dolvatm
      DIMENSION  ratm(3,*),ratm0(3,*)
c
c Local:
      DIMENSION  strgpa(3,3), oldprim(3,3), epsilon(3,3)
      DIMENSION  vtmp(3)
      DIMENSION  ipiv(3)
      DIMENSION  delprim(9),optprim(9),wt3(9),wt4(3,3)
      CHARACTER*8  strlbl
      LOGICAL    debug_cell
      DIMENSION  ioptaxis(3)
c
c >>>> EXECUTABLE CODE:
c
      iucstat = 0
      if( ndim.eq.0 .or. .not. do_cell ) RETURN
c
c Set flag for debug output messages:
      debug_cell = .false.
 9020 format(1x,a,3f20.10)
 9030 format(1x,5x,3f12.6)
c
      call FLGETIWR( IWR )
c
      write(IWR,9010)
 9010 format(/1x,79('='),
     $ /5x,'Cell optimization driver v1.0 ',
     $ /5x,'Written by R.P. Muller - April 2002',
     $ /5x,'Revised by P.A. Schultz - May 2002',
     $ /5x,'Last update:  6Apr13-PAS (independent components)',
     $ /1x,49('-') )
 9090 format(1x,'End cell optimization step' / 79('=') / )
c
c Convert stress into the appropriate units for dimension ...
      call MKZERO( 9, strgpa )
      call STR2GPA( ndim, rprim, strtot, strgpa, strlbl )
c
c And symmetrize matrix:
      do  jd=2,ndim
        do  id=1,jd-1
          strij = ( strgpa(id,jd) + strgpa(jd,id) ) / two
          strgpa(id,jd) = strij
          strgpa(jd,id) = strij
        enddo
      enddo
c
      if( debug_cell )then
c      ... and output the stress, both raw, and converted:
        write(IWR,*) 'Total stress tensor (scf) in Ry and ',strlbl
        do  jd=1,ndim
          if( ndim.eq.1 )then
            write(IWR,'(2(2x, f13.6))') (strtot(id,jd),id=1,ndim),
     $                                  (strgpa(id,jd),id=1,ndim)
          elseif( ndim.eq.2 )then
            write(IWR,'(2(2x,2f13.6))') (strtot(id,jd),id=1,ndim),
     $                                  (strgpa(id,jd),id=1,ndim)
          elseif( ndim.eq.3 )then
            write(IWR,'(2(2x,3f13.6))') (strtot(id,jd),id=1,ndim),
     $                                  (strgpa(id,jd),id=1,ndim)
          endif
        enddo
      endif
      call QHIST( ihistfl, 'STRE' , ndim , rdata, strgpa )
c
c Rework stress to apply constraints to cell optimization:
c
      if( iuctype .eq. 1 )then
        write(IWR,*) 'Isotropic volume expansion only'
        press = ( strgpa(1,1) + strgpa(2,2) + strgpa(3,3) ) / three
        call MKZERO( 9, strgpa )
        strgpa(1,1) = press
        strgpa(2,2) = press
        strgpa(3,3) = press
c
      elseif( iuctype .eq. 2 )then
        write(IWR,*) 'Fixed-volume cell optimization'
c       Here, we remove the net pressure from the stress tensor
        press = ( strgpa(1,1) + strgpa(2,2) + strgpa(3,3) ) / three
        do  idim=1,3
          strgpa(idim,idim) = strgpa(idim,idim) - press
        enddo
c
      elseif( iuctype .ge. 3 .and. iuctype .le. 9 )then
        write(IWR,*) 'Enforce axial constraints'
c       Off-diagonal terms are zeroed ...
        do  idim=1,ndim-1
          do  jdim=idim+1,ndim
            strgpa(idim,jdim) = zero
            strgpa(jdim,idim) = zero
          enddo
        enddo
c        ... and axial constraints were binary-encoded into the option
        iucopt = iuctype - 2
        ioptaxis(1) = 1
        if( AND(iucopt,1) .eq. 0 )then
c         Freeze out x-optimization
          ioptaxis(1) = 0
          strgpa(1,1) = zero
        endif
        iucopt = iucopt/2
        ioptaxis(2) = 1
        if( AND(iucopt,1) .eq. 0 )then
c         Freeze out y-optimization
          ioptaxis(2) = 0
          strgpa(2,2) = zero
        endif
        iucopt = iucopt/2
        ioptaxis(3) = 1
        if( AND(iucopt,1) .eq. 0 )then
c         Freeze out z-optimization
          ioptaxis(3) = 0
          strgpa(3,3) = zero
        endif
        call STOPXERR('uc axial constraints not ready')
      endif
c
c Add external desired stress to scf stress to get stress to
c be eliminated in cell optimization
c
      strxsum = zero
      do  jd=1,ndim
        do  id=1,ndim
          strgpa(id,jd) = strgpa(id,jd) + strxtrnl(id,jd)
          strxsum = strxsum + ABS( strxtrnl(id,jd) )
        enddo
      enddo
      if( strxsum.gt. 1.d-6 )then
c       Write out stress (different from scf) being eliminated:
        write(IWR,*) 'Net stress tensor, external stress, in ',strlbl
        do  jd=1,ndim
          if( ndim.eq.1 )then
            write(IWR,'(2(2x, f13.6))') (strgpa(id,jd),id=1,ndim),
     $                                  (strxtrnl(id,jd),id=1,ndim)
          elseif( ndim.eq.2 )then
            write(IWR,'(2(2x,2f13.6))') (strgpa(id,jd),id=1,ndim),
     $                                  (strxtrnl(id,jd),id=1,ndim)
          elseif( ndim.eq.3 )then
            write(IWR,'(2(2x,3f13.6))') (strgpa(id,jd),id=1,ndim),
     $                                  (strxtrnl(id,jd),id=1,ndim)
          endif
        enddo
      endif
c
c Test for cell stress convergence:
c
      ijd = 0
      strs_this = zero
      do  jd=1,ndim
        do  id=1,ndim
          ijd = ijd + 1
          wt3(ijd) = strgpa(id,jd)
          if( ABS( strgpa(id,jd) ) .gt. ABS( strs_this ) )
     $     strs_this = strgpa(id,jd)
        enddo
      enddo
      nrlxdim = ijd
c
      call CONVCHK( IWR, 'stress', iucstat,iucstep,maxucstep,cellconv,
     $ strsmax,strsrms,  1, nrlxdim, wt3 )
c
c     We have the absolute value of largest stress in strsmax
c     and its signed value in strs_this, which we will use for different purposes
c
      call QHIST( ihistfl, 'MAXS' , ndata, strsmax, wt3 )
      call QHIST( ihistfl, 'RMSS' , ndata, strsrms, wt3 )
c
c     iucstat = 0(cont); 1(steps); 2(conv)
c
      if( debug_cell ) write(IWR,'(a,3i6)') 'iucstat,step,iter=',
     $                                    iucstat,iucstep,iuciter
c
      if( iucstat .ne. 2 )then
c
c       Not converged - generate new cell
c
c       Save the old lattice vectors in oldprim
        do  jd=1,3
          do  id=1,3
            oldprim(id,jd) = rprim(id,jd)
          enddo
        enddo
c
        if( iucstep.eq.1 )then
          ucstepf = ucstepfac(1,1)
c         On initial step, limit strain to less than strnmax
          epsmax = ABS( ucstepf*strsmax )
          if( epsmax .gt. strnmax )then
c           Shrink step size
            ascal = (strnmax/epsmax)
            call DSCAL( 9, ascal, ucstepfac,1 )
C            ucstepfac = ucstepf*( strnmax / epsmax )
            write(IWR,9020) 'ucstepfac reduced to (1/GPa):',
     $                       ucstepf*ascal
          endif
        endif
c
c       icellscheme==0: Start with steepest descent, and switch
c                       to broyden when strsmax < strsbroy
c       icellscheme==1: Do steepest descent
c       icellscheme==2: Do broyden
c
        if( (icellscheme.eq.0) .and. (strsmax.lt.strsbroy) )
     $   icellscheme = 2
c
        if( (icellscheme.eq.0) .or. (icellscheme.eq.1) )then
c
c         ***** STEEPEST DESCENT *****
c
          if( iucstep.gt.1 )then
c           Past first step: adjust the step size (compare stress)
c           (sorry, the energy is the wrong criterion here)
c           (but, oops, so is the absolute value of the stress!)
            do  id=1,ndim
              do  jd=id,ndim
                strs_this = strgpa(id,jd)
                strs_last = strlast(id,jd)
                ucstepf = ucstepfac(id,jd)
                if( strs_this*strs_last .lt. zero )then
c                 Step overshot (stress reversed):
                  if( ABS(strs_this) .gt. 0.400*ABS(strs_last)
     $             .and. ABS(strs_this) .gt. 0.400*cellconv )then
c                   Step overshot by more than ~50%, decrease step fac
                    ucstepf = ucstepf*cellstepdecr
                  endif
                elseif( ABS(strs_this) .gt. 0.400*ABS(strs_last)
     $            .and. ABS(strs_this) .gt. 0.400*cellconv )then
c                 Step undershot, and by a lot, increase step increment
                  ucstepf = ucstepf*cellstepincr
                endif
                ucstepfac(id,jd) = ucstepf
                if( id.ne.jd) ucstepfac(jd,id) = ucstepf
              enddo
            enddo
          endif
c         Record the current stress for use in next step update
          call DCOPY( 9, strgpa,1, strlast,1 )
c
          write(IWR,'(a,6f10.6)') 'steepest descent ucstepf (1/GPa)=',
     $                  ((ucstepfac(id,jd),jd=id,ndim),id=1,ndim)
c
c         Form the strain matrix epsilon = (I+h\sigma) and
c         then multiply by the lattice vectors
c
          call MKZERO( 9, epsilon )
          do  jd=1,ndim
            epsilon(jd,jd) = one
            do  id=1,ndim
              epsilon(id,jd) = epsilon(id,jd)
     $                       + ucstepfac(id,jd)*strgpa(id,jd)
            enddo
          enddo
          if( debug_cell )then
            write(IWR,*) 'Strain matrix (steep d.) before symmetrize'
            do  id=1,ndim
              write(IWR,'(3f11.6)') (epsilon(id,jd),jd=1,ndim)
            enddo
          endif
c         Generate candidate new lattice vectors L' in rprim
          do  id=1,ndim
            do  jd=1,ndim
              rprim(id,jd) = zero
              do  kd=1,ndim
                rprim(id,jd) = rprim(id,jd)
     $                       + epsilon(id,kd)*oldprim(kd,jd)
               enddo
            enddo
          enddo
c
c         Reset iterative method step counter:
          iuciter = 1
c
c         End steepest descent update
        else
c
c         ***** BROYDEN *****
c
          write(IWR,*) 'cell opt: broyden iteration#',iuciter
c
c         Generate candidate strain matrix
          call MKZERO( 9, epsilon )
          do  jd=1,ndim
            epsilon(jd,jd) = one
            do  id=1,ndim
              epsilon(id,jd) = epsilon(id,jd)
     $                       + ucstepfac(id,jd)*strgpa(id,jd)
            enddo
          enddo
          if( debug_cell )then
            write(IWR,*) 'Strain matrix (broyden) before symmetrize'
            do  id=1,ndim
              write(IWR,'(3f11.6)') (epsilon(id,jd),jd=1,ndim)
            enddo
          endif
c         Generate candidate new lattice vectors L' in rprim
          do  id=1,ndim
            do  jd=1,ndim
              rprim(id,jd) = zero
              do  kd=1,ndim
                rprim(id,jd) = rprim(id,jd)
     $                       + epsilon(id,kd)*oldprim(kd,jd)
              enddo
            enddo
          enddo
c
c         Pack rprim, and cell change matrix (L'-L) for blend
          chgmax = zero
          ijd = 0
          do  jd=1,ndim
            do  id=1,ndim
              ijd = ijd + 1
              optprim(ijd) = oldprim(id,jd)
              delprim(ijd) = rprim(id,jd) - oldprim(id,jd)
              chg = ABS( delprim(ijd) )
              if( chg.gt.chgmax ) chgmax = chg
            enddo
          enddo
          nblend = ijd
c
          call BROYDDJ( IWR, 'cell', iuciter,nhistuc,
     $     ucblend,chgmax,nblend,
     $     delprim, optprim,  wt3,  wt4 )
c -->      vv1-i    vv2-io    vv3-s vv4-s
c
          iuciter = iuciter + 1
c
c         Unpack new primitive vectors:
          ijd = 0
          do  jd=1,ndim
            do  id=1,ndim
              ijd = ijd + 1
              rprim(id,jd) = optprim(ijd)
            enddo
          enddo
c
c         End broyden update
        endif
c
        if( iuctype .eq. 2 )then
          write(IWR,*) 'Fixed volume update enforced'
c         This is the iffy one, as update will not be clean.
          call CROSS( rprim(1,2),rprim(1,3), vtmp )
          volold = ABS( oldprim(1,1)*vtmp(1) + oldprim(2,1)*vtmp(2)
     $                + oldprim(3,1)*vtmp(3) )
          call CROSS( rprim(1,2),rprim(1,3), vtmp )
          volnew = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $                + rprim(3,1)*vtmp(3) )
          volchange = volold/volnew
          write(IWR,*) 'Volume change factor=',volchange
          afactor = volchange**(0.3333333333333333d0)
          call DSCAL( 9, afactor, rprim,1 )
        endif
c
c       Purify the move based on any symmetry; reuse the epsilon matrix
c       for the strain matrix
c
        if( nsym.ge.1 )then
c         Copy matrices into temp space and take transpose
c         We want to solve e*L = L'. If we take the transpose of both
c         sides, we get (e*L)T = (L)T*(e)T = (L')T. Thus, we can solve
c         (L)T*e=(L')T for e, since e is symmetric.
c
          if( debug_cell )then
            write(IWR,*) 'Old and (raw) new lattice vectors:'
            do  jd=1,ndim
              write(IWR,'(6f13.6)') (oldprim(id,jd),id=1,ndim),
     $                                (rprim(id,jd),id=1,ndim)
            enddo
          endif
c
c         Load full 3x3 matrix for use in dgesv:
          do  id=1,3
            do  jd=1,3
              epsilon(id,jd) = rprim(jd,id)
              wt4(id,jd) = oldprim(jd,id)
            enddo
          enddo
c
c          Solve L*epsilon = L' for the strain matrix "epsilon"
           info = 0
           call DGESV( 3,3, wt4,3,ipiv,epsilon,3,info )
c
c          Purify the strain matrix (in "epsilon")
           if( debug_cell .and. nsym.gt.1 )then
             write(IWR,*) 'Strain matrix before symmetrize (again)'
             do  id=1,ndim
               write(IWR,'(3f11.6)') (epsilon(id,jd),jd=1,ndim)
             enddo
           endif
           call SYMSTR( ndim, epsilon, nsym, rmatsym )
c
c          Apply constraints to get squeeze what we want into this:
           if( iuctype .eq. 1 )then
             write(IWR,*) 'Isotropic cell expansion enforced'
             ascale = ( epsilon(1,1) + epsilon(2,2) + epsilon(3,3) )
     $                                / three
             dif = ABS( ascale-epsilon(1,1) )
     $           + ABS( ascale-epsilon(2,2) )
     $           + ABS( ascale-epsilon(3,3) )
             if( dif .gt. 1.d-4 )then
c              Something went wrong, tell everyone:
               write(IWR,*) 'Strain matrix for isotropic update:'
               do  id=1,ndim
                 write(IWR,'(3f11.6)') (epsilon(id,jd),jd=1,ndim)
               enddo
               call STOPXERR( 'isotropic cell update not isotropic' )
             endif
             call MKZERO( 9, epsilon )
             epsilon(1,1) = ascale
             epsilon(2,2) = ascale
             epsilon(3,3) = ascale
           elseif( iuctype .eq. 2 )then
             write(IWR,*) 'Fixed volume update not checked'
c            This is the iffy one, as update will not be clean.
C             call CROSS( rprim(1,2),rprim(1,3), vtmp )
C             volold = ABS( oldprim(1,1)*vtmp(1) + oldprim(2,1)*vtmp(2)
C     $                   + oldprim(3,1)*vtmp(3) )
C             call CROSS( rprim(1,2),rprim(1,3), vtmp )
C             volnew = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
C     $                   + rprim(3,1)*vtmp(3) )
C             volchange = volold/volnew
C             write(IWR,*) 'Volume change factor=',volchange
C             afactor = volchange**(0.3333333333333333d0)
C             call DSCAL( 9, afactor, epsilon,1 )
           elseif( iuctype .ge. 3 .and. iuctype .le. 9 )then
             write(IWR,*) 'Axial constraints enforced'
             do  id=1,nd
               if( ioptaxis(id).eq.0 ) epsilon(id,id) = zero
               do  jd=1,nd
                 if( id.ne.jd ) epsilon(id,jd) = zero
               enddo
             enddo
           endif
c
           if( debug_cell )then
             write(IWR,*) 'Strain matrix after symmetrization'
             do  id=1,ndim
               write(IWR,'(3f11.6)') (epsilon(id,jd),jd=1,ndim)
             enddo
           endif
c
c          Reapply epsilon to L (rprim): L' = L*epsilon
           do  id=1,ndim
             do  jd=1,ndim
               rprim(id,jd) = zero
               do  kd=1,ndim
                 rprim(id,jd) = rprim(id,jd) +
     $            epsilon(id,kd)*oldprim(kd,jd)
              enddo
            enddo
          enddo
c
c         End symmetry purification
        endif
c
        if( debug_cell )then
          write(IWR,*) 'Old and (purified) new lattice vectors'
          do  jd=1,ndim
            write(IWR,'(6f13.6)') (oldprim(id,jd),id=1,ndim),
     $                              (rprim(id,jd),id=1,ndim)
          enddo
c         Leave a record of the cell volume for new cell:
          call CROSS( rprim(1,2),rprim(1,3), vtmp )
          volum = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $               + rprim(3,1)*vtmp(3) )
          write(IWR,'(1x,a,f16.6)') 'New Cell volume(bohr^3)=',volum
        endif
c
c       Update the atomic coords, k-points, etc.:
        call CELLRUPD( do_cell, iucstat, dolvatm,
     $   ndim, rprim, oldprim, orig,rscale,
     $   natm,ratm,ratm0 )
c
c       Print output and iterate
c
        iucstep = iucstep + 1
c
        if( debug_cell ) write(IWR,'(a8,i3,6f13.6)')
     $   'CellOpt',iucstep,engytotl,strsmax,
     $    strs_this,strs_last,oldprim(1,1),rprim(1,1)
c
c       Save this max stress for next step descent criterion
        strs_last = strs_this
c
      else
c       Converged!
c
        if( iucstat .eq. 2 )then
          write(IWR,*) 'CellOpt: Converged'
        elseif( iucstat .eq. 1 )then
          write(IWR,*) 'CellOpt: Not converged - out of iterations'
        endif
c
        if( debug_cell ) write(IWR,'(a8,i3,5f13.6/)')
     $   'CellOpt', iucstep,engytotl,zero,
     $           strsmax,zero,rprim(1,1)
c
        if( debug_cell )then
          write(IWR,*)
          write(IWR,*) 'Final Cell Parameters'
          do  jd=1,ndim
            write(IWR,'(5x,3f20.10)') (rprim(id,jd),id=1,ndim)
          enddo
c         Leave a record of the cell volume for new cell:
          call CROSS( rprim(1,2),rprim(1,3), vtmp )
          volum = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $               + rprim(3,1)*vtmp(3) )
          write(IWR,'(1x,a,f16.6)') 'Final Cell volume(bohr^3)=',volum
          write(IWR,*)
          write(IWR,*) 'Final Stress Tensor'
          do  jd=1,ndim
            write(IWR,'(5x,3f20.10)') (strgpa(id,jd),id=1,ndim)
          enddo
        endif
c
      endif
c
      write(IWR,9090)
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> UCSTAMP
c
c
      subroutine UCSTAMP( IWR, iucstat, iucstep,iuciter, rprim )
c---------------------------------------------------------------
c Purpose: write the new cell vector to the output file
c
c Written: Peter A. Schultz, 23-Septemeber-2002, v2.62
c
c Revision history:
c  23Sep08-PAS/2.62: extracted from main routine; 23Jan11-add LAT out
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  rprim(3,3)
c Local:
      DIMENSION  vtmp(3)
c
c >>>> EXECUTABLE CODE:
c
      if( IWR.lt. 0 ) goto 999
c
      write(IWR,9430) 
 9430 format(1x,79('#') )
c
      if( iucstat.lt.2 )then
        write(IWR,9531) 'NEW',iucstep,iuciter
      else
        write(IWR,9531) 'FINAL',iucstep, iuciter
      endif
 9531 format(1x,a,' UNIT CELL for step',i4,', (blend step',i4,')')
c
      write(IWR,*) 'primitive lattice vectors'
      write(IWR,9532)  ((rprim(id,jd),id=1,3),jd=1,3)
 9532 format(5x,3f20.10)
c
      write(IWR,9430)
c
c     Leave a record of the cell volume for new cell:
      call CROSS( rprim(1,2),rprim(1,3), vtmp )
      volum = ABS( rprim(1,1)*vtmp(1) + rprim(2,1)*vtmp(2)
     $           + rprim(3,1)*vtmp(3) )
      if( iucstat .lt. 2 )then
c       Not converged
        write(IWR,'(1x,a,f16.6)') 'New Cell volume(bohr^3)=',volum
      else
c       Final converged volume
        write(IWR,'(1x,a,f16.6)') 'FINAL Cell volume(bohr^3)=',volum
      endif
c
c     Put out lattice vectors in condensed form:
      write(IWR,'(1x,a,        3(1x,3f9.4)                )')
     $              'NEWCELL', ((rprim(i,j),i=1,3),j=1,3)
c
  999 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
