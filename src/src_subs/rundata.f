
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUNDATA
c
c
      subroutine RUNDATA( IDAT,IWR,  doflag,
     $ itstart,itstop,meth_bl,nhiste,igesopt,
     $ ioccdft,iflcdft, elecno, eigpop,
     $ convgr,convs,convsl,convii, conv1c, alfast,
     $ icloseocc,
     $ edegen,etemp, scfblnd,scfbl2,scfconv,
     $ nspin,nk,nstate, norb,nocc0k,
     $ opt_spin, meth_sp,conv_sp,nstep_sp,blend_sp,
     $ dorelax, natm,nafrc,iforce,ifdefct,
     $          atmnm, natmnm,
     $          nslics,islics,islicats,vslics,
     $          iatframe,
     $          nhistg,gblend,gconv,tstep,igstart,igstop,
     $ do_cell, do_pvcalc, ndim, iuctype, icellscheme,strsbroy,
     $          iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $          strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $          strlast,
     $ idosolv, scut,
     $ nbranchd,nbranch,bvl,bsymb )
c---------------------------------------------------------------
c Purpose: read in (optional) data for Quest run
c
c Revision history:
c   7Dev22-PAS/2.69: Koopmans/occdft capability
c   3Aug17-PAS/2.68: unified slics data passing (68c) 7Jan20 meth{,bl}_asd
c   5Jun15-PAS/2.67: option to select atoms by indices
c   8Apr15-PAS/2.66: e-blend options; lnwrit lines
c  24Sep13-PAS/2.65a: pass flag for enthalpy (E+PV) calc
c  03May12-ACP+AHE/2.64: Band structure input (PAS-cosmetics)
c  16Nov10-PAS/2.62: tweak diagnostic output
c  10Dec08-PAS/2.61j: fix error in read of spsteps
c  23Sep08-PAS/2.62: option to freeze atoms along vector; isotropic cell
c  21Jan08-PAS/2.61: igesopt moved to rundata
c  10Jan08-PAS/2.61: silence SLIC unless active; spin optimization
c  11Jul07-PAS/2.60: disable/replace "last iter"/"percent" keywords
c   4Oct05-RPM/2.59: idosolv and parameters set from input file
c  14Feb05-PAS/2.59: defect k-sampling
c   7May04-RPM/2.58: add MD input stuff
c  25Apr04-PAS/2.58: purge "igtype" variable
c  25Nov03-PAS/2.58: remove scf blend reset feature
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  14Jan03-PAS/2.55: limited-history into elec, geom, and cell
c  28Jan02-PAS/2.54: pass fixed-frame atom data
c  16May02-PAS/2.53: cutoff param input moved out of "setdata"
c   3May02-PAS/2.53: install cell optimization option input
c   2Oct01-PAS/2.50: damped dynamics params
c  21Jun01-PAS/2.48: replace STOPs
c   5May01-PAS/2.47: clean echo, transfer out defaults/initialize
c   3May00-PAS/2.44: scf blend "scfbl2" for follow-on geometries
c                    occup temperature moved here; blend keyword
c   3Mar00-PAS/2.43: full branching options list
c  31Mar99-PAS/2.33: clean out green fcn stuff
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  26Oct98-PAS/2.27: restart geometry run, and end iter restart.
c  14Jul98-PAS/2.23: blend chol option; optionalize iter nums.
c  31Mar98-PAS/2.22: fix bug - passing natm into geomdat, ges
c  17Sep97-PAS/2.21: blend restart iteration, input formats
c  14Apr94-PAS/2.10
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c currently, the doflag option is unassigned
      LOGICAL    doflag
c geometry parameter stuff:
      LOGICAL    dorelax
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
      DIMENSION  iforce(*)
cRMV slic projected coordinate stuff
      DIMENSION  islics(*),islicats(*),vslics(*)
      DIMENSION  iatframe(3)
c cell parameter stuff
      LOGICAL    do_cell, do_pvcalc
      DIMENSION  strxtrnl(*), strlast(3,3),ucstepfac(3,3)
c spin optimization stuff
      LOGICAL    opt_spin
c occdft stuff, elecno, and room for eigpop(nstate,nk,nspin)
      DIMENSION  elecno(*)
      DIMENSION  eigpop(*)
c cutoff info:
      DIMENSION  convs(8)
c band structure data
      CHARACTER*3  bvl
      CHARACTER*2  bsymb(nbranchd)
c
c Local declarations:
      PARAMETER  ( nline = 128 )
      CHARACTER    label*(128), readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
      DIMENSION    vtmp(4)
      DATA         zero / 0.d0 /
c Defines blending type
      DATA    methbl_broy,methbl_asd / 1,2 /
c
c >>>> EXECUTABLE CODE:
c
c Input/output data formats:
 8000 format(a72)
 9000 format(1x,a)
 9010 format(1x,3f12.8)
 9030 format(1x,1p4d16.8)
 9130 format(1x,'***** Format change:',
     $ /5x,'Keyword >',a,'< disabled, replaced by',
     $ /5x,'keyword >',a,'<',
     $ /5x,'Modify/update input file and submit again')
 9210 format(1x,2i4)
c
      wantlbl = 'run ph'
      read(IDAT,8000)  label
      call LNWRIT( IWR, label, nline )
      if( readlbl.ne.wantlbl ) goto 1300
c
  100 continue
c
      read(IDAT,8000,err=999,end=999)  label
      call LNWRIT( IWR, label, nline )
c
      if    ( readlbl.eq. 'end ru' .or.
     $        readlbl.eq. 'end_ru' .or.
     $        readlbl.eq. 'end of' )then
c       Finished with run phase data:
        goto 999
c
      elseif( readlbl.eq. 'do fla' )then
c       Set free flag ON
        doflag = .true.
      elseif( readlbl.eq. 'no fla' )then
c       Set free flag OFF
        doflag = .false.
c
      elseif( readlbl.eq. 'first ' )then
c       Set initial iteration number
        read(IDAT,*   )  itstart
        write(IWR,9210)  itstart
c
      elseif( readlbl.eq. 'iterat' )then
c       Set final iteration number
        read(IDAT,*   )  itstop
        write(IWR,9210)  itstop
c
      elseif( readlbl.eq. 'last i' )then
c       Disabled - catch and notify
        write(IWR,9130) 'last iteration','iterations'
        call STOPXERR( 'Outdated keyword, use >iterations<' )
c
      elseif( readlbl.eq. 'blmeth' )then
c       Select a Method for Blending
        read(IDAT,8000,err=999,end=999)  label
        call LNWRIT( IWR, label, nline )
        call STRPARS( label, 1,60, il1,il2 )
c       We assume here that method is a 6-char string
        wantlbl = label(il1:il1+5)
        if( wantlbl.eq.'BROYDE' .or. wantlbl.eq.'Broyde' .or.
     $      wantlbl.eq.'broyde' )then
          meth_bl = methbl_broy
        elseif( wantlbl.eq.'ASD' .or. wantlbl.eq.'asd' )then
          meth_bl = methbl_asd
          nhiste = 0
        endif
c
      elseif( readlbl.eq. 'histor' )then
c       Number of steps to keep in scf blending history
        read(IDAT,*   )  nhiste
        if( meth_bl .eq. methbl_asd )then
          write(IWR,*) '# For ASD blend, nhiste set to 0'
          nhiste = 0
        endif
        write(IWR,9210)  nhiste
c
      elseif( readlbl.eq. 'restar' )then
c       Restart blend scheme at specified iteration number
        write(IWR,*) '***** restart blend feature removed'
        write(IWR,*) '***** use history feature instead'
        call STOPXERR( 'restart blend feature removed' )
c
      elseif( readlbl.eq. 'states' )then
c       Number of states for eigensolve/density computation
        read(IDAT,*   )  nstate
        write(IWR,9210)  nstate
        if( ioccdft.gt.0 )then
c         Oops, this must come before occ-dft setup
          write(IWR,*) '***** ERROR: must have states before occdft'
          call STOPXERR( 'states must be specified before occdft' )
        endif
c
      elseif( readlbl.eq. 'temper' )then
c       Temperature for level occupations
        read(IDAT,*   )  etemp
        write(IWR,9010)  etemp
c
      elseif( readlbl.eq. 'kclose' )then
c       Set almost closed occs - fixed e-count within k-pt
        icloseocc = 1
      elseif( readlbl.eq. 'closed' )then
c       Set closed occs  - 0K, fixed e-count within k-pt
        icloseocc = 2
      elseif( readlbl.eq. 'occdft')then
c       Read occ-dft setup input, either from full file or k=1 setup
c       The existing file should have the full k-points
c       If not, we build a nk=1 version here, expand later when know nk
        nk1 = 1
        call SETOCCDFT( IDAT,IWR, ioccdft,iflcdft,
     $   nspin, nk1, nstate,  elecno,
     $   eigpop )
c
      elseif( readlbl.eq. 'edegen' )then
c       Energy difference criterion for degeneracy of states
        read(IDAT,*   )  edegen
        write(IWR,9010)  edegen
c
      elseif( readlbl.eq. 'no ges' )then
c       Simple overlapping atom guess for new geom scf.
        igesopt = 0
      elseif( readlbl.eq. 'do ges' )then
c       Use full dHam guess for scf ges for new geom scf.
        igesopt = 2
      elseif( readlbl.eq. 'guess ' )then
c       Specify type of scf guess for new geometry
        read(IDAT,*)  igesopt
        write(IWR,'(1x,i4,5x,a)')  igesopt,' = option for scf ges'

      elseif( readlbl.eq. 'blend ' )then
c       Blend ratio for initial scf step
        read(IDAT,*   )  scfblnd
        write(IWR,9010)  scfblnd
      elseif( readlbl.eq. 'percen' )then
c       Disable - catch and notify
        write(IWR,9130) 'percent blend','blend ratio'
        call STOPXERR( 'Outdated keyword, use >blend ratio<' )
c
      elseif( readlbl.eq. 'scfbl2' )then
c       Blend ratio for initial scf step, for subsequent geoms
        read(IDAT,*   )  scfbl2
        write(IWR,9010)  scfbl2
c
      elseif( readlbl.eq. 'conver' )then
c       Convergence criterion for scf (del-Ham)
        read(IDAT,*   )  scfconv
        write(IWR,9010)  scfconv
c
      elseif( readlbl.eq. 'alfast' )then
c       Boundary between slow/fast gaussians
        read(IDAT,*   )  alfast
        write(IWR,9010)  alfast
      elseif( readlbl.eq. 'idosol' )then
c       idosolv controls which solver is used in the eigensolver
        read(IDAT,* ) idosolv
        write(IWR,9210) idosolv
      elseif( readlbl(1:4).eq. 'scut')then
c       scut is the value below which eigenvalues of S are discarded
c        when idosolv=3
        read(IDAT,* ) scut
        write(IWR,9010) scut
c
      elseif( readlbl(1:3).eq.'cut' )then
c       Various cutoffs for integrals
c
        if(     readlbl .eq. 'cutii ' )then
          read(IDAT,*   )  convii
          write(IWR,9030)  convii
        elseif( readlbl .eq. 'cutslo' )then
          read(IDAT,*   )  convsl
          write(IWR,9030)  convsl
        elseif( readlbl .eq. 'cut2s ' )then
          read(IDAT,*   )  conv2s
          write(IWR,9030)  conv2s
        elseif( readlbl .eq. 'cutset' )then
          read(IDAT,*   )  cutfac
          write(IWR,9030)  cutfac
          do  ic=1,4
            convs(ic) = cutfac*convs(ic)
          enddo
        elseif( readlbl .eq. 'cutfrc' )then
          read(IDAT,*   )  cutfac
          write(IWR,9030)  cutfac
          do  ic=5,8
            convs(ic) = cutfac*convs(ic)
          enddo
        elseif( readlbl .eq. 'cutfac' )then
          read(IDAT,*   )  cutfac
          write(IWR,9030)  cutfac
          do  ic=1,8
            convs(ic) = cutfac*convs(ic)
          enddo
        elseif( readlbl .eq. 'cutall' )then
C          read(IDAT,*   )  (convs(ic),ic=1,8)
          read(IDAT,*   )  (vtmp(i),i=1,4)
          write(IWR,9030)  (vtmp(i),i=1,4)
          do  i=1,4
            if( vtmp(i).gt.zero ) convs(i) = vtmp(i)
          enddo
          read(IDAT,*    )  (vtmp(i),i=1,4)
          write(IWR,9030)  (vtmp(i),i=1,4)
          do  i=1,4
            if( vtmp(i).gt.zero ) convs(i+4) = vtmp(i)
          enddo
        elseif( readlbl .eq. 'cutgrd' )then
          read(IDAT,*   )  convgr
          write(IWR,9030)  convgr
        elseif( readlbl .eq. 'cut1c ' )then
          read(IDAT,*   )  conv1c
          write(IWR,9030)  conv1c
        else
c         Do not recognize cut option
          goto 1399
        endif
c
c ***************************
c  Spin optimization options
c ***************************
c
      elseif( readlbl.eq. 'spconv' )then
c       Convergence criterion (polorization, in e) for spin opt
        read(IDAT,*   ) rdat
        write(IWR,9010) rdat
        if( rdat .gt. zero ) conv_sp = rdat
c
      elseif( readlbl.eq. 'spblen' )then
c       Number of steps allowed for spin optimization
        read(IDAT,*   ) rdat
        write(IWR,9010) rdat
        if( rdat .gt. zero .and. rdat .lt. 10.d0 ) blend_sp = rdat
c
      elseif( readlbl.eq. 'spstep' )then
c       Step limit for spin optimization
        read(IDAT,*   ) itmp
        write(IWR,9210) itmp
        if( itmp .gt. 0 ) nstep_sp = itmp
c
      elseif( readlbl.eq. 'spmeth' )then
c       Select a Method for spin optimization
        read(IDAT,8000,err=999,end=999)  label
        call LNWRIT( IWR, label, nline )
        call STRPARS( label, 1,60, il1,il2 )
c       We assume here that method is a 6-char string
        wantlbl = label(il1:il1+5)
        if( wantlbl.eq.'LINEAR' .or. wantlbl.eq.'linear' .or.
     $   wantlbl.eq.'Linear' )then
          meth_sp = 1
        else
c         Safest is to default to "SIMPLE" method
          meth_sp = 0
        endif
c
c ***************************
c  Geometry update options
c ***************************
c
      elseif( readlbl.eq. 'geomet' )then
c       Geometry relaxation parameter specifications:
        call GEOMDAT( IWR,IDAT, dorelax,
     $   natm,nafrc,iforce,ifdefct,
     $   atmnm, natmnm,
     $   nhistg,gblend,gconv,tstep,igstart,igstop,
     $   nslics,islics,islicats,vslics,   iatframe )
c
      elseif( readlbl .eq. 'md dat' .or.
     $        readlbl .eq. 'dynami' )then
c       Molecular dynamics input (either "md data" or "dynamics":
        call QMDDAT( IWR,IDAT )
c
c ***************************
c  Cell optimization options
c ***************************
c
      elseif( readlbl(1:5).eq. 'cell ' )then
c       Cell optimization parameter specifications:
        call CELLDAT( IWR,IDAT,
     $   do_cell, do_pvcalc, ndim, iuctype, icellscheme,strsbroy,
     $   iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $   strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $   strlast )
c
      elseif( readlbl(1:5).eq.'bands') then
c       Band structure data
        call BANDDATA( IDAT,IWR, nbranchd,nbranch,bvl,bsymb )
c
c ***************************
c  And the alternative ...
c ***************************
c
      else
c       Do not know what is happening:
        call STOPXERR( 'rundata - unknown input command' )
c
      endif
c
      goto 100
c
  999 continue
c
c ***************************
c  Finished, do checks
c ***************************
c
c Prevent stupidities, 21jan93-PAS:
      if( nstate.gt.norb ) nstate = MIN( nocc0k+10, norb )
      if( nstate.lt.nocc0k )then
        write(IWR,*) '>>>>> electron pairs, states=',nocc0k,nstate
        call STOPXERR( 'nstate-m - too few occupied states' )
      endif
c
c Record the cutoffs to the output file:
      write(IWR,'(a)') 'cutXXX for code (grd,slo,ii,1c)'
      write(IWR,9030)   convgr,convsl,convii,conv1c
      write(IWR,'(a)') 'cutall for code (s,t,nl,loc;frc)'
      write(IWR,9030)  (convs(ic),ic=1,8)
      if( icloseocc.eq.2 )
     $write(IWR,'(a,f12.6)') 'etemperature(Ry)=',etemp
c
      write(IWR,*)
      write(IWR,*)  '##### output #####'
      write(IWR,*)
c
      RETURN
c
c >>>>>> ERROR HANDLING:
c
 1300 continue
      write(IWR,*) 'Want label >',wantlbl,'<, found >',readlbl,'<'
 1399 continue
      call STOPXERR( 'rundata - bad input' )
c
c    That's all Folks!
c
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
