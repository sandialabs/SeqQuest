c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYF
c
c
c---------------------------------------------------------------
c Purpose: routines to handle the split-file Broyden file stuff
c          to enable truncating histories and 2GB file limits
c Written: Peter A. Schultz, 10-January-2003
c Revision history:
c  none
c---------------------------------------------------------------
c
      subroutine BROYF( IWR, itertot,itersav, itmaxbr, msize_it,
     $ nit_file,maxfl_broy, nfl_broy )
c---------------------------------------------------------------
c Purpose: determine number of files to divide Broyden into,
c          based on length of retained history, and 2GB limits
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c    2GB file limit = 2**31 = 2 147 483 648 bytes
c                   = 2**23 =   268 435 456 r*8 array spaces
      PARAMETER  ( maxflsiz =   268 435 456 )
c
      DIMENSION  ibrfls(*)
      CHARACTER  broyfl*80, broynum*3
      CHARACTER  sufx*(*)
c
c >>>> EXECUTABLE CODE:
c
c Determine limit on number of iterations per Broyden history file ...
c  ... number of iterations that fit into size limit of file
      nit_file = maxflsiz / msize_it
c  ... limit single file length to length of single complete history:
      if( nit_file.gt.itmaxbr ) nit_file = itmaxbr
c
c Determine number of files needed for broyden:
c  nfl_tot = number of files for all iterations:
      nfl_tot = (itertot - 1) / nit_file + 1
c  nfl_sav = # of files needed for select reduced number of iters:
      nfl_sav = (itersav - 2) / nit_file + 2
c  nfl_broy = # of files required
      nfl_broy = MIN( nfl_tot, nfl_sav )
      if( nfl_broy .ne. 1 ) write(IWR,*) ' BROY: nfiles=',nfl_broy
      if( nfl_broy.lt.1 .or. nfl_broy.gt.maxfl_broy )
     $ call STOPXERR( 'nfl_broy' )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYF0
c
c
      entry BROYF0( ibroyfl,
     $ ifl_broy,it_file,                   ibrfls )
c---------------------------------------------------------------
c Purpose: set broyden file pointer to very beginning
c---------------------------------------------------------------
c
      ifl_broy = 1
      ibroyfl = ibrfls(ifl_broy)
      it_file = 0
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYFSKIP
c
c
      entry BROYFSKIP( itprune,it_skip, ibroyfl,
     $ ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
c---------------------------------------------------------------
c Purpose: locate and set files to initial step in active blend
c          Skip through "itprune" steps to get to first active iter
c---------------------------------------------------------------
c
c Number of total files to skip through:
      nfl_skip = itprune / nit_file
c Number of file sets that are repeated
      nh_skip = nfl_skip / nfl_broy
c The file number in the current set of files:
      ifl_broy = nfl_skip - nh_skip*nfl_broy + 1
c Identify file unit number, and index in file of last skipped iter
      ibroyfl = ibrfls(ifl_broy)
      it_file = itprune - nfl_skip*nit_file
c
      it_skip = it_file
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYFIT
c
c
      entry BROYFIT( ibroyfl,
     $ ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
c---------------------------------------------------------------
c Purpose: step to next iteration in broyden history 
c---------------------------------------------------------------
c
      it_file = it_file + 1
      if( it_file.gt.nit_file )then
c       This file done, go to next file in sequence:
        ifl_broy = ifl_broy + 1
        if( ifl_broy.gt.nfl_broy ) ifl_broy = 1
        ibroyfl = ibrfls(ifl_broy)
        REWIND( ibroyfl )
        it_file = 1
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYFOPEN
c
c
      entry BROYFOPEN( sufx, nfl_broy, ibrfls )
c---------------------------------------------------------------
c Purpose: open files used to store Broyden history
c---------------------------------------------------------------
c
      do  ifl_br=1,nfl_broy
        call STRNUM( nstr, broynum, ifl_br, nfl_broy )
        broyfl  = 'broy' // broynum(1:nstr) // '_' // sufx
c
        call FLOPENB( ibrfl , broyfl  )
c
        ibrfls(ifl_br) = ibrfl
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYFCLOSE
c
c
      entry BROYFCLOSE( nfl_broy, ibrfls )
c---------------------------------------------------------------
c Purpose: close all files used to store Broyden history
c---------------------------------------------------------------
c
      do  ifl_br=1,nfl_broy
        ibrfl = ibrfls(ifl_br)
        call FLCLOSE( ibrfl )
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
