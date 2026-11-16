c*******************************************************************************
c NEBSTAT - routines to try and handle image-related NEB failure gracefully
c C Routines: c  NEBSTAT  - the neb_master receiving status messages
c  NEBSTATOK- the image master reports to neb_master that all is ok
c  NEBSTATX - the image master reports failure to neb_master
c
c Written: Peter A. Schultz, 10-April-2015, for 2.66
c
c Revision history:
c  none
c*******************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBSTAT
c
c
      subroutine NEBSTAT( nodei, image, ierr, icomm )
c---------------------------------------------------------------
c Purpose: Neb_master checking status
c---------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER  nodei, image, ierr, icomm
c Local:
      INTEGER  lenmsg, msgneb(2)
c
      lenmsg = 2
      call MPRECVI( nodei, lenmsg, msgneb, icomm )
c
      image = msgneb(1)
      ierr  = msgneb(2)
c
c This might be the place to convert the error code to a useful message,
c  given that tne he nebstatx routine in here did the complement.
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBSTATOK
c
c
      subroutine NEBSTATOK( neb_master, image, icomm )
c---------------------------------------------------------------
c Purpose: Neb_master checking status
c---------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER  neb_master, image, ierr, icomm
c Local:
      INTEGER  lenmsg, msgneb(2)
c
      ierr = 0
      msgneb(1) = image
      msgneb(2) = ierr
c
      lenmsg = 2
      call MPSENDI( neb_master, lenmsg, msgneb, icomm )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBSTATX
c
c
      subroutine NEBSTATX( IWRNEB, statneb, image )
c---------------------------------------------------------------
c Purpose: report failure of image to neb_master
c
c Written: Peter A. Schultz, 10-April-2015, for 2.66
c
c Revision history:
c  none
c---------------------------------------------------------------
c Notes:
c  This routine is intended to have neb_master find out about errors
c   on remote images, and to have it, at least, record the failure
c   so that the job (and the user) is not in the dark about a hang.
c  The motivating issue is an scf failure on a remote image, which
c   causes the remote image to shut down, but the rest of the NEB,
c   including the neb_master, hung waiting to hear from this image,
c   wasting cpu and leaving the user in the dark.
c  The NEB in its full set of paths is ugly machinery, and it might
c   be rather challenging to shut the whole thing down cleanly.
C  Here we just want to guarantee that at least the master shuts down
c   cleanly, leaving a message for the user to know what happened.
c  Given that a failure might be reported from anywhere in the code
c   where there is a failure, it may not know enough to say
c   much other than it bombed (it might not know the image#).
c
c  This particular routine does not stop the run.
c  It assumes that maing code will stop after more diagnostics.
c  Purpose here is simply to leave user a useful error message.
c
      IMPLICIT NONE
      INTEGER  image, IWRNEB
      CHARACTER*(*)  statneb
c Local:
      INTEGER  neb_master, nodei, icomm, master, iproc, IWR, ierr
      INTEGER  lenmsg, msgneb(2)
c
c >>>> EXECUTABLE CODE
c
      call MPNODE0_G( neb_master )
      call MPNODE_G( nodei )
      call MPCOMM_G( icomm )
      call MPNODE0( master )
      call MPNODE( iproc )
c
      call FLGETIWR( IWR )
      if( nodei .eq. neb_master )then
c       THe image on the neb_master failed.
c       Record to the neb_out, and let main code handle this.
        write(IWR   ,*) 'Failure in image on master node'
        write(IWR   ,*) 'NEB therefore fails'
        write(IWR   ,*) 'image, failure code=',image,statneb
        write(IWRNEB,*) 'Failure in image on master node'
        write(IWRNEB,*) 'NEB therefore fails'
        write(IWRNEB,*) 'image, failure code=',image,statneb
      elseif( iproc .eq. master )then
c       Report failure to neb_master
        if( iproc .eq. master )then
          write(IWR,*) 'Failure on this image in NEB'
          write(IWR,*) 'Reporting failure to neb_master'
          write(IWR,*) 'image, failure code=',image,statneb
          msgneb(1) = image
          if( statneb .eq. 'scf' )then
            ierr = 1
          else
c           Default generic error
            ierr = -1
          endif
          msgneb(2) = ierr
          lenmsg = 2
          call MPSENDI( neb_master, lenmsg, msgneb, icomm )
        endif
      endif
c
c  Try to respond to last (ierror) sync in NEB
      call MPBCASTI( neb_master, 1, ierr, icomm )
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
