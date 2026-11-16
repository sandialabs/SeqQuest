c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUNOPT
c
c
      subroutine RUNOPT( IDAT,IWR, lvlout,mkrhopt, kparopt,
     $ dotests, dosetup,doiters,doforce,dorelax,do_cell,do_pvcalc,
     $ opt_spin, do_md, do_neb, do_post,
     $ do_bands, do_optic, do_nonscf,
     $ dosplit,doblas3, redusym, dolocxc,
     $ dospars,doeigen,doprjct,dorspsy,dodefct,
     $ do_kppsolve, do_psolv,
     $ ikoopman
     $ )
c---------------------------------------------------------------
c Purpose: options read in for setup phase, bulk code
c
c Revision history:
c  22Dec22-PAS/2.69: koopmans flag
c  24Sep13-PAS/2.65a: offer flag for enthalpy (E+PV) calculation
c  01May12-ACP/2.63: kppsolve (parallel-parallel) flag (do_kppsolve)
c  01May12-ACP+AHE/2.63: non-SCF calc flags (do_bands, do_optic)
c  22Nov08-PAS/2.63: k-parallelism flag; 05Dec12-PAS: psolv flag
c  10Jan08-PAS/2.61: spin optimization
c  11Jul07-PAS/2.60: disable/replace "input data" keyword
c   7May04-RPM/2.58: add MD flag
c  23Jul03-PAS/2.56: add option to set record length
c   4Jan03-PAS/2.55: add "setup data" = "input data"
c   5Jan02-PAS/2.52: do_cell flag input
c  21JUn01-PAS/2.48: replace STOPs
c   2May01-PAS/2.47: echo output level control, add do_neb
c   4May00-PAS/2.44: added Raczkowski's flags for real/recip
c                    space sym and linear scaling defect calc.
c  10Jun99-PAS/2.36: flags for blas3,lin scaling
c  14May99-PAS/2.35: option force dm/ef grid density build
c  31Mar99-PAS/2.33: purge green fcn, pass locxc and madel
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  29Jun98-PAS/2.23: flags for split basis and eigen/linscale
c  23Jun98-PAS/2.23: output level control
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL   dotests
      LOGICAL   dosetup,doiters,doforce,dorelax,do_cell,do_pvcalc
      LOGICAL   opt_spin, do_md, do_neb, do_post
      LOGICAL   do_bands, do_optic, do_nonscf
      LOGICAL   dosplit,doeigen,doblas3, redusym, dolocxc
      LOGICAL   dospars,doprjct,dorspsy,dodefct
      LOGICAL   do_kppsolve, do_psolv
c
c local declarations:
      CHARACTER  label*80, readlbl*8
      EQUIVALENCE  (label,readlbl)
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a)
 9000 format(1x,a)
 9130 format(1x,'***** Format change:'
     $ /5x,'The keyword >',a,'< has been replaced by'
     $ /5x,'the keyword >',a,'<'
     $ /5x,'Modify/update input file and submit again')
c
      write(IWR,'(a)') '##### Run options:'
c
  100 read(IDAT,8000)  label
      write(IWR,9000)  label
c
      if(     readlbl .eq. 'setup da' )then
        goto 990
c
      elseif( readlbl .eq. 'input da' )then
c       Disable - catch and notify
        write(IWR,9130) 'input data','setup data'
        call STOPXERR( 'Outdated keyword, use >setup data<' )
c
      elseif( readlbl .eq. 'set reco' )then
c       Change maximum record length in i/o
        read(IDAT,*)  maxrecl
        write(IWR,'(i10)')  maxrecl
        call IOSETRECL( maxrecl )
c
      elseif( readlbl .eq. 'output l' )then
        read(IDAT,*)  lvlout
        write(IWR,'(i5)')  lvlout
c
c  Compel type of grid density build?
      elseif( readlbl .eq. 'do dmrho' )then
        mkrhopt = 1
      elseif( readlbl .eq. 'no dmrho' )then
        mkrhopt = 2
      elseif( readlbl .eq. 'do efrho' )then
        mkrhopt = 2
      elseif( readlbl .eq. 'no efrho' )then
        mkrhopt = 1
c
c  Compel or forbid k-parallelism:
      elseif( readlbl .eq. 'do kpara' )then
        kparopt = -2
C        call KPFLAG_SET( kparopt )
      elseif( readlbl .eq. 'no kpara' )then
        kparopt = -1
C        call KPFLAG_SET( kparopt )
c
c  Choose to compel series of parallel solves?
      elseif( readlbl .eq. 'do psolv' )then
        do_psolv = .true.
      elseif( readlbl .eq. 'no psolv' )then
        do_psolv = .false.
c
      elseif( readlbl .eq. 'do tests' )then
        dotests = .true.
      elseif( readlbl .eq. 'no tests' )then
        dotests = .false.
c
      elseif( readlbl .eq. 'do setup' )then
        dosetup = .true.
      elseif( readlbl .eq. 'no setup' )then
        dosetup = .false.
c
      elseif( readlbl .eq. 'do iters' )then
        doiters = .true.
      elseif( readlbl .eq. 'no iters' )then
        doiters = .false.
c
      elseif( readlbl .eq. 'do force' )then
        doforce = .true.
      elseif( readlbl .eq. 'no force' )then
        doforce = .false.
c
      elseif( readlbl .eq. 'do relax' )then
        dorelax = .true.
      elseif( readlbl .eq. 'no relax' )then
        dorelax = .false.
c
      elseif( readlbl .eq. 'do cell ' )then
        do_cell = .true.
      elseif( readlbl .eq. 'no cell ' )then
        do_cell = .false.
c
      elseif( readlbl .eq. 'do entha' )then
        do_pvcalc = .true.
      elseif( readlbl .eq. 'no entha' )then
        do_pvcalc = .false.
c
      elseif( readlbl .eq. 'do neb  ' )then
        do_neb  = .true.
      elseif( readlbl .eq. 'no neb  ' )then
        do_neb  = .false.
c
      elseif( readlbl .eq. 'do spino' )then
        opt_spin = .true.
      elseif( readlbl .eq. 'no spino' )then
        opt_spin = .false.
c
cqmd: do md OR do dynamics?
      elseif( readlbl .eq. 'do md   ' .or.
     $        readlbl .eq. 'do dynam' )then
        do_md = .true.
        if (.not. doforce)
     $          write(IWR,*) "Turning on doforce because of dynamics"
        doforce = .true.
        doiters = .true. 
        dosetup = .true.
      elseif( readlbl .eq. 'no md   ' .or.
     $        readlbl .eq. 'no dynam' )then
        do_md = .false.
c
      elseif( readlbl .eq. 'do post ' )then
        do_post = .true.
      elseif( readlbl .eq. 'no post ' )then
        do_post = .false.
c
      elseif(readlbl .eq. 'do bands') then
           do_bands=.true.
           do_nonscf=.true.
      elseif(readlbl .eq. 'no bands') then
           do_bands=.false.
c
      elseif(readlbl.eq. 'do optic')then
         do_optic=.true.
         do_nonscf=.true.
      elseif(readlbl .eq. 'no optic') then
         do_optic=.false.
c
      elseif( readlbl .eq. 'do split' )then
        dosplit = .true.
      elseif( readlbl .eq. 'no split' )then
        dosplit = .false.
c
      elseif( readlbl .eq. 'do blas3' )then
        doblas3 = .true.
      elseif( readlbl .eq. 'no blas3' )then
        doblas3 = .false.
c
      elseif( readlbl .eq. 'do kppso' .or.
     $        readlbl .eq. 'do kpsol' )then
        do_kppsolve = .true.
      elseif( readlbl .eq. 'no kppso' .or.
     $        readlbl .eq. 'no kpsol' )then
        do_kppsolve = .false.
c
      elseif( readlbl .eq. 'redusym ' )then
        redusym = .true.
      elseif( readlbl .eq. 'keepsym ' )then
        redusym = .false.
c
      elseif( readlbl .eq. 'do locxc' )then
        dolocxc = .true.
      elseif( readlbl .eq. 'no locxc' )then
        dolocxc = .false.
c
c
      elseif( readlbl .eq. 'do spars' )then
        dospars = .true.
      elseif( readlbl .eq. 'no spars' )then
        dospars = .false.
c
      elseif( readlbl .eq. 'do eigen' )then
        doeigen = .true.
      elseif( readlbl .eq. 'no eigen' )then
        doeigen = .false.
c
      elseif( readlbl .eq. 'do prjct' )then
        doprjct = .true.
      elseif( readlbl .eq. 'no prjct' )then
        doprjct = .false.
c
      elseif( readlbl .eq. 'do rspsy' )then
        dorspsy = .true.
      elseif( readlbl .eq. 'no rspsy' )then
        dorspsy = .false.
c
      elseif( readlbl .eq. 'do defct' )then
        dodefct = .true.
      elseif( readlbl .eq. 'no defct' )then
        dodefct = .false.
c
      elseif( readlbl .eq. 'do koopm' )then
c       Koopmans - likely add a filename input here eventually
        ikoopman = 1
c       Compel a DM-density (not EF), to avoid gridrho/EF issue
        mkrhopt = 1
      elseif( readlbl .eq. 'no koopm' )then
        ikoopman = 0
c
      else
        write(IWR,*) '>>>>> Wanted: [do|no] [setup|iters|force|relax]'
        write(IWR,*) '>>>>> or [setup data], but got instead: ',readlbl
        call STOPXERR( 'optlabel/ bad runopt command option' )
      endif
c
c Go back and try to find more command options ...
      goto 100
c
  990 continue
c
      write(IWR,*) '>>>>> Options (set/scf/frc/relax/cell/neb/md)=',
     $  dosetup,doiters,doforce,dorelax,do_cell,do_neb,do_md
      write(IWR,*) '>>>>> use blas3=',doblas3
      write(IWR,*) '>>>>> split basis =',dosplit
      write(IWR,*) '>>>>> Output level set =',lvlout
c
c     For now, we will assure that if there are no band-structure
c     or optical properties requested, that do_nonscf is false
c
      if(.not.do_bands .and. .not. do_optic) do_nonscf=.false.
      write(IWR,*) '>>>> do_bands, do_optic, do_nonscf=',
     $  do_bands, do_optic, do_nonscf

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
c End:
