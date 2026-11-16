c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CELLDAT
c
c
      subroutine CELLDAT( IWR,IDAT,
     $ do_cell, do_pvcalc, ndim, iuctype, icellscheme,strsbroy,
     $ iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $ strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $ strlast )
c---------------------------------------------------------------
c Purpose: input options affecting cell relaxation
c
c Written: Peter A. Schultz, 3-May-2002, for v2.53 (cell opt)
c
c Revision history:
c   9Apr15-PAS/2.66: lnwrite lines
c  24Sep13-PAS/2.65a: install flag for enthalpy(E+PV) calc
c   6Apr130PAS/2.65: axial-based constrained cell options
c  23Sep08-PAS/2.62: isotropic cell opt
c  11Jul07-PAS/2.60: cosmetic
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  10Jan03-PAS/2.55: limited-history option
c  15May02-PAS/2.53: patch new changes
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL    do_cell, do_pvcalc
      DIMENSION  strxtrnl(3,3)
      DIMENSION  strlast(3*3), ucstepfac(3*3)
c
c Local:
      PARAMETER  ( nline = 128 )
      CHARACTER    label*(nline), readlbl*6, wantlbl*6, glabel*8
      EQUIVALENCE  (label,readlbl)
      DIMENSION    v(3)
c
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a128)
 9000 format(1x,a)
c
  100 continue
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
c      ***** Standard options *****
c
      if    ( readlbl.eq. 'end ce' .or.
     $        readlbl.eq. 'end_ce' )then
c       End of cell parameter input
        goto 990
c
      elseif( readlbl.eq. 'ucstep' )then
c       Maximum number of steps in cell relaxation
        read(IDAT,*)  maxucstep
        write(IWR,'(1x,i4)')  maxucstep
        if( maxucstep.le.0 ) call STOPXERR( 'uc steps max bad' )
c
      elseif( readlbl.eq. 'uc_con' )then
c       Converged stress (GPa) for cell relaxation
        read(IDAT,*)  fac
        write(IWR,'(1x, f15.8, a, f15.8, a )')
     $   fac, ' GPa (default=', cellconv, ' GPa)'
        if( fac.le.zero ) call STOPXERR( 'uc_conv < 0' )
        cellconv = fac
c
      elseif( readlbl.eq. 'isotro' )then
c       Just do an isotropic (fixed-shape) cell volume optimization
        iuctype = 1
c
      elseif( readlbl.eq. 'constr' )then
c       Choose some limiting constraints for cell shape optimization
        read(IDAT,8000)  label
        call LNWRIT( IWR, label, nline )
        call STRPARS( label, 1,60, il1,il2 )
        wantlbl = label(il1:il1+5)
        if(     wantlbl.eq.'isotro' )then
          iuctype = 1
        elseif( wantlbl.eq.'vfixed' )then
          iuctype = 2
c
c       NB: Urelax/Ufixed iuctype parsed via bitmap in cellopt
        elseif( wantlbl.eq.'xrelax' )then
          iuctype = 3
        elseif( wantlbl.eq.'xfixed' )then
          iuctype = 8
        elseif( wantlbl.eq.'yrelax' )then
          iuctype = 4
        elseif( wantlbl.eq.'yfixed' )then
          iuctype = 7
        elseif( wantlbl.eq.'zrelax' )then
          iuctype = 6
        elseif( wantlbl.eq.'zfixed' )then
          iuctype = 5
        elseif( wantlbl.eq.'arelax' )then
          iuctype = 9
        endif
c
      elseif( readlbl.eq. 'pressu' )then
c       Input external pressure (in GPa)
        call MKZERO( 9, strxtrnl )
        read(IDAT,*)  pressure
        write(IWR,'(1x,f15.8,a)')  pressure, ' GPa'
c       Hydrostatic pressure, on diagonal of stress tensor
        do  id=1,ndim
          strxtrnl(id,id) = pressure
        enddo
c       If doing pressure, will want enthalpy
        do_pvcalc = .true.
c
      elseif( readlbl.eq. 'uniaxi' )then
c       Input external uniaxial pressures (in GPa)
        call MKZERO( 9, strxtrnl )
        read(IDAT,*)  (v(id),id=1,ndim)
        write(IWR,'(1x,3f15.8,a)')  (v(id),id=1,ndim), ' GPa'
        do  id=1,ndim
          strxtrnl(id,id) = v(id)
        enddo
c
      elseif( readlbl.eq. 'stress' )then
c       Input external stress (in GPa)
        call MKZERO( 9, strxtrnl )
        read(IDAT,*)  pressure
        write(IWR,'(1x,f15.8,a)')  pressure, ' GPa'
        do  jd=1,ndim
          read(IDAT,*,err=1301)    (strxtrnl(id,jd),id=1,ndim)
          write(IWR,'(1x,3f15.8)') (strxtrnl(id,jd),id=1,ndim)
        enddo
c       Check that this stress tensor is symmetric:
        do  jd=1,ndim
          do id=jd+1,ndim
            if( strxtrnl(id,jd) .ne. strxtrnl(jd,id) )
     $       call STOPXERR( 'tsk,tsk - external stress not symmetric!' )
          enddo
        enddo
c
c      ***** Slightly non-standard options *****
c
      elseif( readlbl.eq. 'uchist' )then
c       Specify history to be used in blend opt methods
        read(IDAT,*)  nhistuc
c
      elseif( readlbl.eq. 'ucmeth' )then
c       Specify cell optimization method
        read(IDAT,8000,err=1300,end=1300)  label
        call LNWRIT( IWR, label, nline )
        call STRPARS( label, 1,60, ig1,ig2 )
c       We assume here that glabel is an 8-char string
        glabel = label(ig1:ig1+7)
        if( glabel.eq. 'BEST    ' )then
c         Choose steepest-to-broyden switch
          icellscheme = 0
        elseif( glabel.eq. 'STEEPEST' )then
c         Force steepest descent method
          icellscheme = 1
        elseif( glabel.eq. 'BROYDEN ' )then
c         Force broyden
          icellscheme = 2
        else
c         Unknown - stick with "best" default
          write(IWR,*) '***** Unknown cell opt method'
          write(IWR,*) 'Choices are: BEST,STEEPEST,BROYDEN'
          write(IWR,*) 'Default to BEST=steep-to-broy scheme'
          icellscheme = 0
        endif
c
      elseif( readlbl.eq. 'str_br' )then
c       Stress criterion to convert to broyden from steepest
        read(IDAT,*) fac
        write(IWR,'(1x,f15.8,a)')  fac, ' GPa'
        strsbroy = fac
c
c      ***** Non-standard options *****
c
      elseif( readlbl.eq. 'max_st' )then
c       Size limit on strain in initial cell step (unitless)
        read(IDAT,*)  fac
        write(IWR,'(1x,f15.8)')  fac
        if( fac.le.zero .or. fac.gt.0.5 )
     $   call STOPXERR( 'cell strain input is out-of-bounds' )
        strnmax = fac
c
      elseif( readlbl.eq. 'cell_f' )then
c       Stress (GPa) to strain factor for initial cell update (1/GPa)
        read(IDAT,*)  fac
        write(IWR,'(1x, f15.8, a, f15.8, a )')
     $   fac, ' GPa^-1 (default=', ucstepfac, ')'
        if( fac .ge. zero )
     $   call STOPXERR( 'cell_factor must be strictly negative')
        do  ij=1,3*3
          ucstepfac(ij) = fac
        enddo
c
      elseif( readlbl.eq. 'bulk_m' )then
c       Bulk_modulus guess for setting step factor (GPa)
        read(IDAT,*)  bmod
        write(IWR,'(1x,f15.8,a)')  bmod, ' GPa'
        if( bmod.le.zero ) call STOPXERR( 'bulk_modulus < 0' )
        fac = - ( one / (3.d0*bmod) )
        do  ij=1,3*3
          ucstepfac(ij) = fac
        enddo
c
      elseif( readlbl.eq. 'uc_ble' )then
c       Blend factor for initial broyden step (should be 1.0)
        read(IDAT,*)  fac
        write(IWR,'(1x,f15.8)')  fac
        if( fac.le.zero ) call STOPXERR( 'uc_blend < 0' )
        ucblend = fac
c
      else
c
        call STOPXERR( 'celldat - unknown keyword' )
c
      endif
c
      goto 100
c
  990 continue
      write(IWR,9991) cellconv
 9991 format(1x,'>>>> uc_conv (cell convergence)=',f15.8,' GPa' )
c
       RETURN
c
 1300 continue
      write(IWR,*) '***** ERROR: reading cell relax params'
      write(IWR,*) '>>>>> >end cell opt< to end cell input'
      call STOPXERR( 'celldat - error processing cell data' )
 1301 continue
      call STOPXERR( 'cellstr - error reading xtrnl stress' )
      do_cell = .false.
c
      RETURN
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
c End
