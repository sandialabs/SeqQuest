c*******************************************************************************
c GESSER  - routines to handle the scf-guessing machinery
c
C Routines:
c  GESSER   - construct and record guess info for next step
c  GESSET   - construct guess scf Ham for this step
c  GESOPEN  - open files needed for gesser
c  GESCLOSE - closed files used by gesser
c
c Revision history:
c  12Apr15-PAS/2.66: gesopen (use igesopt to determine if to open files)
c*******************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GESSER
c
c
      subroutine GESSER( iges0fl,igesfl, igesopt,
     $ do_geom, madeges, mat, nk, nspin,
     $ vij,       wk1 )
c      nmatskloc  mat ; nmatsloc= mat*nkloc*nspin
c---------------------------------------------------------------
c Purpose: construct guess for scf Hamiltonian in geom relax.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
c  12Nov09-ACP/2.63: k-parallel I/O, PAS: k-master only executes
c   2Feb09-PAS/2.63: HKsp+H memory for k-parallel
c  25Feb07-PAS/2.60: flag-based disable of guess scf in step
c   4Oct02-PAS/2.54: image-specific guess files for NEB
c   3May02-PAS/2.53: make work with cell optimization
c  23Jul01-PAS/2.49: push down from main code, add spin-pol
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER( one = 1.d0 )
c
c Input/Output arrays:
      LOGICAL    do_geom
c      DIMENSION  vij(mat,nk,nspin)
      DIMENSION  vij(mat,*)
c Scratch arrays:
      DIMENSION  wk1(mat)
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
      call KPMINE( nk, nkloc, nk0 )
c
c     Skip, if no k-points to work
      if( nkloc .eq. 0 ) goto 999
c     And only k-masters will be active ...
      if( iproc_k .ne. master_k ) goto 999
c
      if( do_geom .and. igesopt .eq. 2 )then
c
        REWIND( unit=iges0fl )
c
c       Skip if ges disabled/unavailable for this geom step:
        if( madeges .lt. 0 ) goto 90
c
c       Check that ges files in correct state for constructing ges
        read(iges0fl,err=90,end=90)  istatges
        if( istatges .eq. 0 )then
          call STOPXERR( 'GESSER: expected dHam, got dHam0' )
          REWIND( unit=iges0fl )
        endif
        goto 100
c
   90   continue
c       The ges file is empty or corrupted, or we ignore guess
        istatges = 0
c
  100   continue
        madeges = istatges
c
c       Store away current dHam0, stamped as 0-Ham-containing
        istatges = 0
        REWIND( unit=iges0fl )
        write(iges0fl)  istatges
        if( iproc_k .eq. master_k )then
c         Only kp-masters writes since within a k-group vij is the same.
c
          nkspin = 0
          do  ispin=1,nspin
            do  ikloc=1,nkloc
              call WRITBIG( iges0fl, mat, vij(1,ikloc+nkspin) )
            enddo
            nkspin = nkloc
          enddo
c
        endif
c
        if( madeges .gt. 0 )then
c         Construct guess Ham: augment current dHam0 with previous del-Ham
c
          REWIND( unit=igesfl )
c         Only the kp-masters do the big I/O and calculations.
          if( iproc_k .eq. master_k )then
c
            nkspin = 0
            do  ispin=1,nspin
              do  ikloc=1,nkloc
c               Get previous del-Ham  = ( dHam[scf] - dHam0 ) ...
                call READBIG( igesfl, mat, wk1 )
c     .          .. and add to current initial dHam0 to construct ges Ham
                call DAXPY( mat, one, wk1, 1, vij(1,ikloc+nkspin), 1 )
              enddo
              nkspin = nkloc
            enddo
c
          endif
        endif
c
      endif
c
  999 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GESSET
c
c
      subroutine GESSET( iges0fl,igesfl, igesopt,
     $ do_geom, madeges, mat, nk, nspin,
     $ vij, wk1 )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER( one = 1.d0 )
c
c Input/Output arrays:
      LOGICAL    do_geom
c      DIMENSION  vij(mat,nk,nspin)
      DIMENSION  vij(mat,*)
c Scratch arrays:
      DIMENSION  wk1(mat)
c
c >>>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
      call KPMINE( nk, nkloc, nk0 )
c
      if( nkloc .eq. 0 ) goto 1999
      if( iproc_k .ne. master_k ) goto 1999
c
      madeges = 0
      if( do_geom .and. igesopt .eq. 2 )then
c
        if( igesopt.gt.1 .and. do_geom )then
c         Construct del-Ham=(dHam(scf) - dHam0), and save
c
          REWIND( unit=iges0fl )
          read(iges0fl)  istatges
          if( istatges.ne.0 )
     $     call STOPXERR( 'GESSET: expected Ham0 in guess file' )
c
          REWIND( unit=igesfl )
c
          if( iproc_k .eq. master_k )then
c
            nkspin = 0
            do  ispin=1,nspin
              do  ikloc=1,nkloc
c               Retrieve dHam0 for this geometry ...
                call READBIG( iges0fl, mat, wk1 )
c               ... generate del-Ham = dHam[scf] - dHam0
                call DSCAL( mat, -one,  wk1, 1 )
                call DAXPY( mat,  one, vij(1,ikloc+nkspin), 1, wk1,1 )
c               ... and store away the del-Ham for the next geometry
                call WRITBIG( igesfl, mat, wk1 )
              enddo
              nkspin = nkloc
            enddo
          endif
c
c         Stamp dHam0 with record it is used up, del-Ham set ...
          madeges = 1
          istatges = 1
c         Everybody writes the flag.
          REWIND( unit=iges0fl )
          write(iges0fl)  istatges
c
        endif
c
      endif
c
c    That's all Folks!
c
 1999 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GESOPEN
c
c
      subroutine GESOPEN( image, iges0fl,igesfl,igesopt, nk )
c---------------------------------------------------------------
c Purpose: get appropriate ges file for this image
c
c Written: Peter A. Schultz, 3-October-2002, for 2.54 (neb)
c
c Revision history:
c  12Apr15-PAS/2.66: only open files if gesser is turned on
c  12Nov09-ACP/2.62: write separate files at k-local level
c   2Feb09-PAS/2.62: split into ges0+ges files
c---------------------------------------------------------------
c
c  Based on the image number, identify and open the appropriate
c  file for the scf guess for this image geometry
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*10  imgname
      CHARACTER*10  kptname
c
c Every MPI task opens it's own version 
c
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
      call KPMINE( nk, nkloc, nk0 )
c
c Set default unit=-1, to denote no file, so that gesclose knows:
      iges0fl = -1
      igesfl = -1
c     Skip, if no k-points to work
      if( nkloc .eq. 0 ) goto 999
c     Skip, if gessser is not on
      if( igesopt.ne.2 ) goto 999
c     And only k-masters will be active ...
      if( iproc_k .ne. master_k ) goto 999
c     
      if( nprocs_k.eq.0 ) then
        if( image.gt.0 )then
c         Image-specific ges file name
          call STRNUM( nstr, imgname, image,image )
          call FLOPENB( iges0fl, 'gesz'//imgname(1:nstr) )
          call FLOPENB( igesfl, 'ges'//imgname(1:nstr) )
        else
c         Blank name
          call FLOPENB( iges0fl, 'ges0' )
          call FLOPENB( igesfl, 'ges' )
        endif
      else
        call STRNUM( nstrk, kptname, iproc_k, iproc_k+1 )
        if( image.gt.0 )then
c         Image-specific ges file name
          call STRNUM( nstr, imgname, image,image )
          call FLOPENB( iges0fl,
     $         'gesz'//imgname(1:nstr)//'_'//kptname(1:nstrk) )
          call FLOPENB( igesfl,
     $         'ges'//imgname(1:nstr)//'_'//kptname(1:nstrk) )
        else
c         Blank name
          call FLOPENB( iges0fl, 'ges0'//'_'//kptname(1:nstrk) )
          call FLOPENB( igesfl, 'ges'//'_'//kptname(1:nstrk) )
        endif
      endif
c
  999 continue
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GESCLOSE
c
c
      subroutine GESCLOSE( iges0fl,igesfl )
c
c We do not need igesopt here. File units more definitive if files are open
c
      if( iges0fl .ge. 0 )then
c       igesfl exists iff iges0fl exists
        call FLCLOSE( iges0fl )
        call FLCLOSE( igesfl )
      endif
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
c End:
