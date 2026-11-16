c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GEOMDAT
c
c
      subroutine GEOMDAT( IWR,IDAT, dorelax,
     $ natm,nafrc,iforce,ifdefct,
     $ atmnm, natmnm,
     $ nhistg,gblend,gconv,tstep,igstart,igstop,
     $ nslics,islics,islicats,vslics,    iatframe )
c---------------------------------------------------------------
c Purpose: read definition of geometry relaxation
c
c Written: Peter A. Schultz, 26-September-1997, v2.21
c
c Revision history:
c  31Jul17-PAS/2.68: layers and vgfgix into slics
c   5Jun15-PAS/2.67: select atom sequences by atom labals
c   9Apr15-PAS/2.66: lnwrite lines
c   5May11-PAS/2.62: enforce vgfix after grelax/gfixed requirement
c  20Feb09-PAS/2.62: refined options for constraining atoms along vector
c   9Sep08-PAS/2.62: option to freeze atoms along vector
c  21Jan08-PAS/2.61: mv igesopt to rundata
c  25Feb06-PAS/2.60: "do ges" added (to partner with "no ges")
c  31Aug05-RMV/2.60: add options for projected coordinates - SLIC
c  15Mar05-PAS/2.59: skip blank line, rather than treat as error
c   7May04-PAS/2.58: offer MD section here, too.
c  25Apr04-PAS/2.58: purge "igtype"
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  10Jan03-PAS/2.55: limited-history option
c  28Jan02-PAS/2.54: fixed frame atom data
c   6May02-PAS/2.53: reconcile with cell optimization
c   2Oct01-PAS/2.50: damped dynamics options, and time step
c  21Jun01-PAS/2.48: replace STOPs
c   5May01-PAS/2.47: clean echo
c  29Jul00-PAS/2.45: choose atoms to be relaxed
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  26Oct98-PAS/2.27: restart geometry
c   8Apr98-PAS/2.22: add igesopt for carrying ges option
c  31Mar98-PAS/2.21: bug fix for last/first relax atoms
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input declarations
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
c Output declarations
      LOGICAL    dorelax
      DIMENSION  iforce(natm)
cRMV slic data
      DIMENSION  islics(4,*),islicats(*),vslics(3,*)
c
      DIMENSION  iatframe(3)
c
c Local declarations:
      PARAMETER  ( nline = 128 )
      CHARACTER    label*(nline), readlbl*6, glabel*8
      EQUIVALENCE  (label,readlbl)
      DIMENSION    vtmp(3)
      DIMENSION    itmp(3)
c
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a128)
 9000 format(1x,a)
c
c  Set flag for whether force array cleared:
      itypfix = 0
      islicv = 0
      islicat = 0
c
  100 continue
      read(IDAT,8000,err=999,end=999)  label
      call LNWRIT( IWR, label, nline )
c
      if    ( readlbl.eq. 'stop' )then
        write(IWR,*) 'stop  relaxation'
        goto 999
c
      elseif( readlbl.eq. 'end ge' .or.
     $        readlbl.eq. 'end_ge' )then
c       End of geometry parameter input
        goto 990
c
      elseif( readlbl.eq. 'relax ' )then
c       Set relaxation to be on
        dorelax = .true.
        goto 990
c
      elseif( readlbl.eq. 'grelax' )then
c       From all-fixed, choose sequence of atoms to relax:
cpas    Note: need to fix grelax/gfixed/slics to allow general mixing
        if( itypfix .eq.  1 )then
          call STOPXERR( 'grelax-1/do not mix grelax and gfixed' )
        elseif( nslics .ne. 0 )then
          call STOPXERR( 'SLICS must be after all grelax' )
        elseif( itypfix .eq. 0 )then
c         Preset to fix everybody
          nafrc = 0
          do  iatm=1,natm
            iforce(iatm) = 0
          enddo
        endif
        itypfix = -1
c        read(IDAT,*)  irlx1,irlxn
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( 2, label,nline, natm,atmnm,natmnm, itmp )
        irlx1 = itmp(1)
        irlxn = itmp(2)
        write(IWR,'(2i6,a)')  irlx1,irlxn,' => relax atom sequence'
C        call LNWRIT( IWR, label, nline )
        if( irlx1.lt.1 .or. irlxn.gt.natm .or. irlxn.lt.irlx1 )then
          call STOPXERR( 'grelax-2/bad atom sequence' )
        endif
        do  iatm=irlx1,irlxn
          if( iforce(iatm) .eq. 0 )then
            iforce(iatm) = 1
            nafrc = nafrc + 1
          endif
        enddo
c
      elseif( readlbl.eq. 'gfixed' )then
c       From relax-all, fix atoms in place:
cpas    Note: need to fix grelax/gfixed/slics to allow general mixing
        if( itypfix .eq. -1 )then
          call STOPXERR( 'gfixed-1/do not mix gfixed and grelax' )
        elseif( nslics .ne. 0 )then
          call STOPXERR( 'SLICS must be after all gfixed' )
        elseif( itypfix .eq. 0 )then
c         Preset to relax everybody
          nafrc = natm
          do  iatm=1,natm
            iforce(iatm) = 1
          enddo
        endif
        itypfix = 1
C        read(IDAT,*)  irlx1,irlxn
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( 2, label,nline, natm,atmnm,natmnm, itmp )
        irlx1 = itmp(1)
        irlxn = itmp(2)
        write(IWR,'(2i5,a)')  irlx1,irlxn,' => fixed atom sequence'
        if( irlx1.lt.1 .or. irlxn.gt.natm .or. irlxn.lt.irlx1 )then
          call STOPXERR( 'gfixed-2/bad atom sequence' )
        endif
        do  iatm=irlx1,irlxn
          if( iforce(iatm) .eq. 1 )then
            iforce(iatm) = 0
            nafrc = nafrc - 1
          endif
        enddo
c
      elseif( readlbl(1:4) .eq. 'bond' )then
c       Constrain bond length
        islictyp = 4
        nslics = nslics + 1
        if( nslics.gt.natm ) call STOPXERR('Too many constraints')
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( 2     , label,nline, natm,atmnm,natmnm,
     $                itmp )
        write(IWR,'(1x,16i5)')  (itmp(i),i=1,2)
        islic = nslics
        islics(1,islic) = islictyp
        islics(2,islic) = itmp(1)
        islics(3,islic) = itmp(2)
c
        do  ia=1,2
          iatm = itmp(ia)
          if( iforce(iatm).eq.0 ) nafrc = nafrc + 1
          if( iforce(iatm).gt.1 ) call STOPXERR('constraint conflict')
          iforce(iatm) = islictyp
        enddo
c
      elseif( readlbl(1:5) .eq. 'vgfix' )then
c       Apply constraints connecting atoms along a vector
        if( readlbl(6:6) .eq. 'p' .or. readlbl(6:6) .eq. 'e' )then
c         Fix atoms in plane, no motion along vector, free in-plane
          islictyp = 5
        elseif( readlbl(6:6) .eq. 'd' )then
c         Fix distance between atom planes along vector, free in-plane
          islictyp = 6
        elseif( readlbl(6:6) .eq. 'a' )then
c         Limit atoms along axis, with free distance between planes
          islictyp = 7
        elseif( readlbl(6:6) .eq. 'l' )then
c         Limit atoms along axis, with fixed distance between planes
          islictyp = 8
        else
          call STOPXERR( 'Bad option for vgfixX (X=plane,dist,axis)' )
        endif
c
        read(IDAT,*)  nvgfix, vtmp
        write(IWR,'(i5,3x,3f12.8,a)')  nvgfix, vtmp,
     $   ' = # of atoms constrained with this vector, list is ...'
c
        vmag = vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2
        if( vmag .lt. 0.01 ) call STOPXERR( 'vgfix - zero vector' )
        vmag = DSQRT( vmag )
        do  ir=1,3
          vtmp(ir) = vtmp(ir) / vmag
        enddo
c
        nslics = nslics + 1
        if( nslics.gt.natm ) call STOPXERR('Too many constraints')
c
        isl1 = islicat + 1
C        read(IDAT,*)  (islicats(iatm),iatm=isl1,isln)
        read(IDAT,8000,end=913,err=913)  label
        if(     nvgfix.eq.-1 )then
c         Pick out a (single) sequence of atoms
          call INDATOM( 2     , label,nline, natm,atmnm,natmnm,
     $                  itmp )
          write(IWR,'(1x,16i5)')  itmp(1), itmp(2)
          if( itmp(1).lt.1 .or. itmp(2).lt.itmp(1) .or. itmp(2).gt.natm)
     $     call STOPXERR( 'vgfix - atoms out of boumds' )
          nvgatms = itmp(2) - itmp(1) + 1
          isln = islicat + nvgatms
          if( isln.gt.natm )call STOPXERR( 'Too many constrained atoms')
          iatm = itmp(1)
          do  isl=isl1,isln
            islicats(isl) = iatm
            iatm = iatm + 1
          enddo
        elseif( nvgfix.eq.0 )then
c         Select all the atoms
          nvgatms = natms
          isln = islicat + nvgatms
          if( isln.gt.natm )call STOPXERR( 'Too many constrained atoms')
          iatm = 1
          do  isl=isl1,isln
            islicats(isl) = iatm
            iatm = iatm + 1
          enddo
        elseif( nvgfix.ge.1 .and. nvgfix.le.natm )then
          nvgatms = nvgfix
          isln = islicat + nvgatms
          if( isln.gt.natm )call STOPXERR( 'Too many constrained atoms')
          call INDATOM( nvgfix, label,nline, natm,atmnm,natmnm,
     $                  islicats(isl1) )
          write(IWR,'(1x,16i5)')  (islicats(iatm),iatm=isl1,isln)
        else
          call STOPXERR( 'vgfix - bad number input' )
        endif
c
        isln = islicat + nvgatms
        if( isln.gt.natm )call STOPXERR('Too many constrained atoms')
        islicat = isln
        islic = nslics
        islicv = islicv + 1
c
        islics(1,islic) = islictyp
        islics(2,islic) = nvgatms
        islics(3,islic) = isl1  
        islics(4,islic) = islicv
        vslics(1,islicv) = vtmp(1)
        vslics(2,islicv) = vtmp(2)
        vslics(3,islicv) = vtmp(3)
c
c
        do  ivgatm=isl1,isln
          iatm = islicats(ivgatm)
          if( iatm.lt.1 .or. iatm.gt.natm )then
            call STOPXERR( 'vgfixed atom <1 or >natm')
          endif
c         Record this atom as a plane-fixed atom:
cpas:  vgfix-need to reconcile order of grelax/gfixed/vgfixed
          if( iforce(iatm).eq.0 ) nafrc = nafrc + 1
          if( iforce(iatm).gt.1 ) call STOPXERR('constraint conflict')
          iforce(iatm) = islictyp
        enddo
c
      elseif( readlbl(1:6) .eq. 'slinky' .or.
     $        readlbl(1:6) .eq. 'layers' )then
c       Apply constraints connecting atoms along a vector
        islictyp = 9
c
        read(IDAT,*)  nslink, vtmp
        write(IWR,'(i5,3x,3f12.8,a)')  nslink, vtmp,
     $   ' = # of atoms in slinky layers, list is ...'
c
        vmag = vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2
        if( nslink.lt.1 .or. nslink.gt.natm .or. vmag .lt. 0.01 )then
          call STOPXERR( 'slinky/layers - bad input' )
        endif
        vmag = DSQRT( vmag )
        do  ir=1,3
          vtmp(ir) = vtmp(ir) / vmag
        enddo
c
        nslics = nslics + 1
        isl1 = islicat + 1
        isln = islicat + nslink
        if( nslics.gt.natm ) call STOPXERR('Too many constraints')
        if( isln.gt.natm ) call STOPXERR('Too many total vgfix atoms')
        islicat = isln
        islic = nslics
        islicv = islicv + 1
c
        islics(1,islic) = islictyp
        islics(2,islic) = nslink
        islics(3,islic) = isl1  
        islics(4,islic) = islicv
        vslics(1,islicv) = vtmp(1)
        vslics(2,islicv) = vtmp(2)
        vslics(3,islicv) = vtmp(3)
c
C        read(IDAT,*)  (islicats(iatm),iatm=isl1,isln)
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( nslink, label,nline, natm,atmnm,natmnm,
     $                islicats(isl1) )
        write(IWR,'(1x,16i5)')  (islicats(iatm),iatm=isl1,isln)
c
        do  ivgatm=isl1,isln
          iatm = islicats(ivgatm)
          if( iatm.lt.1 .or. iatm.gt.natm )then
            call STOPXERR( 'slinky atom <1 or >natm')
          endif
c         Record this atom as a plane-fixed atom:
cpas:  vgfix-need to reconcile order of grelax/gfixed/vgfixed
          if( iforce(iatm).eq.0 ) nafrc = nafrc + 1
          if( iforce(iatm).gt.1 ) call STOPXERR('constraint conflict')
          iforce(iatm) = islictyp
        enddo
c
      elseif( readlbl.eq. 'only r' )then
c       Set only first n, or last n, atoms to be relaxed.
        read(IDAT,*)  irelax
        irlx1 = 1
        irlxn = natm
        if( irelax.lt.0 .and. irelax .ge. -natm )then
          write(IWR,'(i4,a)')  irelax,' => last n atoms relaxed'
          irlx1 = natm + irelax + 1
        elseif( irelax.gt.0 .and. irelax.le.natm )then
          write(IWR,'(i4,a)')  irelax,' => first n atoms relaxed'
          irlxn = irelax
        else
          write(IWR,'(i4,a)')  irelax,' => all atoms relaxed'
        endif
c       Preset to relax nobody:
        do  iatm=1,natm
          iforce(iatm) = 0
        enddo
c       Now set to relax specified set of atoms:
        nafrc = 0
        do  iatm=irlx1,irlxn
          iforce(iatm) = 1
          nafrc = nafrc + 1
        enddo
c
      elseif( readlbl.eq. 'frame ' )then
c       Atoms which define a fixed frame
C        read(IDAT,*) iatmfix,iatmvec,iatmpln
        islictyp = 3
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( 3     , label,nline, natm,atmnm,natmnm,
     $                itmp )
        iatmfix = itmp(1)
        iatmvec = itmp(2)
        iatmpln = itmp(3)
        if( iatmfix.ge.1 .and. iatmfix.le.natm .and.
     $      iatmvec.ge.1 .and. iatmvec.le.natm .and.
     $      iatmpln.ge.1 .and. iatmpln.le.natm .and.
     $      iatmfix.ne.iatmvec .and. iatmvec.ne.iatmpln .and.
     $      iatmpln.ne.iatmfix )then
          iatframe(1) = iatmfix
          iatframe(2) = iatmvec
          iatframe(3) = iatmpln
          write(IWR,'(1x,3i5,a)')  iatframe,' =fixed,vec,plane atoms'
        else
          call STOPXERR( 'GEOMDAT: ill-defined frame/bad atoms' )
        endif
        do  isl=1,3
          iatm = itmp(isl)
          if( iforce(iatm).gt.1 ) call STOPXERR('constraint conflict')
          if( iforce(iatm).eq.0 ) nafrc = nafrc + 1
          iforce(iatm) = islictyp
        enddo
c
      elseif( readlbl.eq. 'defect' )then
c       Atom to which defect force is to be applied
C        read(IDAT,*) ifdefct
        read(IDAT,8000,end=913,err=913)  label
        call INDATOM( 1, label,nline, natm,atmnm,natmnm, itmp )
        if( ifdefct.lt.1 .or. ifdefct.gt.natm ) ifdefct = 0
        write(IWR,'(1x,i4)')  ifdefct
c
      elseif( readlbl.eq. 'gblend' )then
c       "Force constant" for initial blend geometry
        read(IDAT,*)  gblend
        write(IWR,'(1x,f15.8)')  gblend
c
      elseif( readlbl.eq. 'gconv ' )then
c       Convergence criterion for geom: largest force magnitude.
        read(IDAT,*)  gconv
        write(IWR,'(1x,f15.8)')  gconv
c
      elseif( readlbl.eq. 'gstart' )then
c       Starting step for geometry relaxation
        read(IDAT,*)  igstart
        write(IWR,'(1x,i4)')  igstart
c
      elseif( readlbl.eq. 'gsteps' )then
c       Total number of steps in geometry relaxation
        read(IDAT,*)  igstop
        write(IWR,'(1x,i4)')  igstop
c
      elseif( readlbl.eq. 'ghisto' )then
c       Total number of steps to use in step updates
        read(IDAT,*)  nhistg
        if( nhistg.lt.0 ) nhistg = 0
        write(IWR,'(1x,i4)')  nhistg
c
      elseif( readlbl.eq. 'no ges' )then
c       Simple overlapping atom guess for new geom scf.
        call STOPXERR( 'no ges - moved to run phase data' )
        igesopt = 0
      elseif( readlbl.eq. 'do ges' )then
c       Use full dHam guess for scf ges for new geom scf.
        call STOPXERR( 'do ges - moved to run phase data' )
        igesopt = 2
      elseif( readlbl.eq. 'guess ' )then
c       Specify type of scf guess for new geometry
        call STOPXERR( 'guess - moved to run phase data' )
        read(IDAT,*)  igesopt
        write(IWR,'(1x,i4,5x,a)')  igesopt,' = option for scf ges'
c
      elseif( readlbl.eq. 'gmetho' )then
c       Pick a method for geometry stepping
        read(IDAT,8000,err=913,end=913)  label
        call LNWRIT( IWR, label, nline )
        call STRPARS( label, 1,60, ig1,ig2 )
c       We assume here that glabel is an 8-char string
        glabel = label(ig1:ig1+7)
        call GTYPSET( glabel )
c
      elseif( readlbl .eq. 'md dat' .or.
     $        readlbl .eq. 'dynami' )then
c        Molecular dynamics input.
         call QMDDAT( IWR,IDAT )
c
      elseif( readlbl.eq. 'timest' )then
c       Get time step for dynamics
        read(IDAT,*)  tstep
        if( tstep .le. zero ) tstep = one
        write(IWR,'(f12.4)')  tstep
c
cRMV:start
cDEV-Use these to reset and communicate slic params
      elseif( readlbl(1:6) .eq. 'slpar1' )then
        read(IDAT,*)  slicpar 
        write(IWR,'(f12.4,a)')  slicpar, ' reset sphere/scalephi'
        call SLICPSET( 'PAR1  ', slicpar,ipslic )
      elseif( readlbl(1:6) .eq. 'slpar2' )then
        read(IDAT,*)  slicpar 
        write(IWR,'(f12.4,a)')  slicpar, ' reset sphere/scaletheta'
        call SLICPSET( 'PAR2  ', slicpar,ipslic )
      elseif( readlbl(1:6) .eq. 'slpar3' )then
        read(IDAT,*)  slicpar 
        write(IWR,'(f12.4,a)')  slicpar, ' is unused'
        call SLICPSET( 'PAR3  ', slicpar,ipslic )
      elseif( readlbl(1:6) .eq. 'slpar4' )then
        read(IDAT,*)  slicpar 
        write(IWR,'(f12.4,a)')  slicpar, ' is unused'
        call SLICPSET( 'PAR4  ', slicpar,ipslic )
c         ... and we have PAR1-PAR9 available in the module
      elseif( readlbl(1:6) .eq. 'fank0 ' )then
        write(IWR,*) 'fank0   - use raw force on sphere anchor'
        call SLICPSET( 'FANK0 ', slicpar,ipslic )
      elseif( readlbl(1:6) .eq. 'fankf ' )then
        write(IWR,*) 'fankf   - add full 1-fold force onto anchor'
        call SLICPSET( 'FANKF ', slicpar,ipslic )
      elseif( readlbl(1:6) .eq. 'fankb ' )then
        write(IWR,*) 'fankb   - add bond-project 1-fold frc onto anchor'
        call SLICPSET( 'FANKB ', slicpar,ipslic )
c
      elseif( readlbl(1:4) .eq. 'slic' )then
c       Projection options for SLIC local cooordinate transforms
c
        write(IWR,*) '# Selected Local Internal Coordinates (SLIC)'
c       Get id of projected atom, type of projection, anchor atoms
C        read(IDAT,*) iatmproj,ityproj,iank1,iank2
        read(readlbl(5:5),'(i1)',err=713)  ityproj
        if( ityproj .eq. 1 .or. ityproj.eq.2 )then
c         1=Cylinder and 2=Sphere are only valid options
          read(IDAT,8000,end=913,err=913)  label
          call INDATOM( 3     , label,nline, natm,atmnm,natmnm,
     $                  itmp )
          iatmproj = itmp(1)
          iank1  = itmp(2)
          iank2  = itmp(3)
          goto 720
        endif
  713   continue
        write(IWR,*) 'ERROR reading SLIC type'
        write(IWR,*) ' Keyword = slic1 .or. slic2  ie. not just slic'
        write(IWR,*) '  projected-atom  anchor1-atom  anchor2-atom '
        call STOPXERR( 'bad slic definition' )
  720   continue
c
c       Set projection flag for atom iatmproj
        write(IWR,'(1x,3i6,2x,a)') iatmproj, iank1, iank2,
     $                          '= atom#,  ank1#, ank2#'
        if( ( iatmproj .le. 0 .or. iatmproj .gt. natm ) .or.
     $      ( iank1    .le. 0 .or. iank1    .gt. natm ) .or.
     $      ( iank2    .le. 0 .or. iank2    .gt. natm ) )then
          call STOPXERR( 'unknown atom(s) in slic' )
        endif
        if( iatmproj.eq.iank1 .or. iatmproj.eq.iank2 .or.
     $      iank1.eq.iank2 )then
          call STOPXERR( 'atoms in slic definition not distinct' )
        endif
c
c:      Load into (not yet integrated  unified-SLICS data
        islictyp = ityproj
        nslics = nslics + 1
        islic = nslics
        islics(1,islic) = islictyp
        islics(2,islic) = iatmproj
        islics(3,islic) = iank1
        islics(4,islic) = iank2
        if( iforce(iatmproj).gt.1 )call STOPXERR('constraints conflict')
        if( iforce(iatmproj).eq.0 ) nafrc = nafrc + 1
c       Set islictyp=2 for iforce(), to distinguish it from other relax
        iforce(iatmproj) = 2
c
cRMV:end
      elseif( readlbl.ne. '      ' )then
c
        call STOPXERR( 'geomdata/unknown input command' )
c
      endif
c
      goto 100
c
  913 continue
      write(IWR,*)  'ERROR reading/parsing geometry section'
      call STOPXERR( 'error in relaxation input' )
c
  990 continue
      write(IWR,9991) igstop,gblend,gconv, nafrc
 9991 format(
     $ /1x,'>>>> Relax params (gstop,gblend,gconv,atoms-relaxd):'
     $ /1x,'>>>> ',i6,1x,f10.7,1x,f15.12,1x,i6 / )
c
       RETURN
c
  999 continue
      dorelax = .false.
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INDATOM
c
c
      subroutine INDATOM( mlist, line,nline, natm,atmnm,natmnm,
     $                    iatmlist )
c---------------------------------------------------------------
c Purpose: extract atom indices from line with atom specs
c
c Written: Peter A. Schultz, 5-June-2015, for 2.67
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c On input:
c     line*(nline) = character string, presumably with atom list
c     mlist = number of presumptive atoms on line
c           = -1 -> open-ended list, then mlist returns total atoms found
c On output:
c     iatmlist(iatm) = indices of atoms found in list
c                    = -1 upon atom-not-found or other error
c     mlist = total number of valid atoms found
c
      IMPLICIT NONE
c
c Input declarations
      INTEGER  mlist,natm,natmnm
      INTEGER  nline
      CHARACTER*(*)  line
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
c Output declaration
      INTEGER  iatmlist(*)
c Local declarations
      INTEGER  maxlist, nlist, iatm
      INTEGER  nl1,nl2,n1,n2,nmlen, nm1,nm2
c
c >>>> EXECUTABLE CODE
c
      if( mlist .eq. 0 )then
c       Nothing to look for
        RETURN
      elseif( mlist .lt. 0 )then
c       Open-ended, check against total number of atoms
        maxlist = natm
      else
        maxlist = mlist
        do  iatm=1,mlist
          iatmlist(iatm) = -1
        enddo
      endif
      nlist = 0
c n1 is first character in string, nlreadtyp self-indexes
      n1 = 1
c
  100 continue
      if( nlist .eq. maxlist ) goto 1000
      call LNREADTYP( line,n1,nline,  iatm,  natm,atmnm,natmnm )
      if( iatm .gt. 0 .and. iatm .le. natm )then
        nlist = nlist + 1
        iatmlist(nlist) = iatm
      else
c       We did not obtain a valid atom label
c       Either end of list, or error reading label in list
        goto 1000
      endif
      goto 100
c
 1000 continue
      if( mlist .eq. -1 ) mlist = nlist
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
