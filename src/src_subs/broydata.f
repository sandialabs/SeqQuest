c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROYDATA
c
c
c---------------------------------------------------------------
c Purpose: mini-module to manage storage and retrieval of
c          history data for broyden blend
c
c Written: Peter A. Schultz, 26-January-2007, for 2.60
c
c Revision history:
c  29Dec19-PAS/2.68a: blsav() array declarations in passthrough routines
c---------------------------------------------------------------
c
      subroutine BROYINIT( incorebr, iters,ithist,nhist,
     $ nvecbr,membroy )
c---------------------------------------------------------------
c Purpose: simple initialization of broyden data storage: do we
c          store on disk or in memory, and basic startup.
c
c Written: Peter A. Schultz, 26-January-2007, for 2.60
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
c
c Module data for file part of broyden data storage:
      PARAMETER  ( maxfl_broy=20 )
      COMMON /BROYFSTUFF/ nfl_broy,nit_file,ifl_broy,it_file,
     $                    ibrfls(maxfl_broy)
c      maxfl_broy = max num of files available for broyden data storage
c      nfl_broy   = number of files which are used to keep broyden data
c      nit_file   = number of iterations (entries) of data on each file
c      ifl_broy   = index (not file unit) of current file in sequence of files
c      it_broy    = local (to file) index of current data entry in current file
c      ibrfls()   = list of file units on which broyden data is stored
c
c Module data for memory part of broyden data storage:
      COMMON /BROYMSTUFF/ i0bl0, ibl0, i0bl2, ibl2, nbl2, itertotl
c      ibl0     = index to track save of vin/dv "broyfl0" data
c      ibl2     = index of first entry of u/dft data in memory
c      nbl2     = number of entries (iterations) of u/dft data in memory
c      itertotl = total (not truncated) index (iteration) of current entry
c
c Input:
c   incorebr < 0 - compel on-disk broyden, else attempt in-core
c   nvecbr  = length of r8 vector on local processor
c   nhist   = length of history (iterations) to be stored
c   membroy = length of storage space available for in-core histories
c
c >>>> EXECCUTABLE CODE:
c
      call FLGETIWR( IWR )
c Default to a disk-based (nfl_broy>0, i0bl0<0) scheme:
      nfl_broy = 1
      i0bl0 = -1
c
      if( incorebr .lt. 1 )then
c       We restrict to a file-based scheme
        if( IWR .ge. 0 ) write(IWR,*) 'EBROYDDJ: on-disk history'
        RETURN
      endif
c
c Determine if we have room to store full broyden history in memory:
      nbr0mem = nvecbr*2
      nbr2mem = nvecbr*2*nhist
      nbrmem = nbr0mem + nbr2mem
c
      if( nbrmem .le. membroy )then
c       We have enough free memory, use in-core broyden history
        nfl_broy = -1
c
c       Set up some defaults:
        i0bl0 = 0
        i0bl2 = i0bl0 + 2
        nbl2 = i0bl2 + 2*nhist
        itertotl = 0
c
c       Reset iters if have in-core history, and restart first iter
c       (i.e., history is empty!)
        if( ithist .eq. 0 )then
          iters = 0
          if( IWR .ge. 0 ) write(IWR,*) 'EBROYDDJ: in-core history'
        endif
      else
c       We keep broyden history on disk
c
        nfl_broy = 1
        i0bl0 = -1
        if( ithist .eq. 0 )then
          if( IWR .ge. 0 ) write(IWR,*) 'EBROYDDJ: on-disk history'
        endif
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY0OPEN
c
c
      subroutine BROY0OPEN( ibrfl0, broyfl )
c---------------------------------------------------------------
c Purpose: open storage location of vin/dv broyden data
c Written: Peter A. Schultz, 26-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
c
c Module data for file part of broyden data storage:
      PARAMETER  ( maxfl_broy=20 )
      COMMON /BROYFSTUFF/ nfl_broy,nit_file,ifl_broy,it_file,
     $                    ibrfls(maxfl_broy)
c Input: file name
      CHARACTER  broyfl*80
c
c >>>> EXECCUTABLE CODE:
c
      if( nfl_broy .gt. 0 )then
c       Open disk file
        call FLOPENB( ibrfl0, broyfl )
      else
c       Set file unit to -1 to signal no file
        ibrfl0 = -1
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY0CLOSE
c
c
      entry BROY0CLOSE( ibrfl0 )
c---------------------------------------------------------------
c Purpose: close file for v/dv storage
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibrfl0 .gt. 0 )then
c       Close disk file
        call FLCLOSE( ibrfl0 )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY0RWND
c
c
      subroutine BROY0RWND( ibrfl0 )
c---------------------------------------------------------------
c Purpose: set storage of vin/dv broyden data to beginning
c Written: Peter A. Schultz, 26-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
c
c Module data for file part of broyden data storage:
      PARAMETER  ( maxfl_broy=20 )
      COMMON /BROYFSTUFF/ nfl_broy,nit_file,ifl_broy,it_file,
     $                    ibrfls(maxfl_broy)
c
c Module data for memory part of broyden data storage:
      COMMON /BROYMSTUFF/ i0bl0, ibl0, i0bl2, ibl2, nbl2, itertotl
c
      DIMENSION  veclocl(*), wbroy(*)
c
c >>>> EXECCUTABLE CODE:
c
      if( ibrfl0 .ge. 0 )then
c       Rewind the disk ...
        REWIND( unit=ibrfl0 )
      else
c       Reset memory pointer to first entry - v()
        ibl0 = 1
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY0PUT
c
c
      entry BROY0PUT( ibrfl0, nveclocl, veclocl, wbroy )
c---------------------------------------------------------------
c Purpose: save v/dv to disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibrfl0 .gt. 0 )then
c       Write to disk ...
        call WRITBIG( ibrfl0, nveclocl, veclocl )
      else
c       Save to memory ...
        ibl = i0bl0 + ibl0
        call STASHVEC( ibl, nveclocl, veclocl, wbroy )
        ibl0 = ibl0 + 1
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY0GET
c
c
      entry BROY0GET( ibrfl0, nveclocl, veclocl, wbroy )
c---------------------------------------------------------------
c Purpose: retrieve v/dv from disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibrfl0 .gt. 0 )then
c       Read from disk ...
        call READBIG( ibrfl0, nveclocl, veclocl )
      else
c       Retrieve from  memory ...
        ibl = i0bl0 + ibl0
        call RETRVVEC( ibl, nveclocl, veclocl, wbroy )
        ibl0 = ibl0 + 1
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STASHVEC
c
c
      subroutine STASHVEC( ibl, nveclocl, veclocl, wbroy )
c---------------------------------------------------------------
c Purpose: save data to memory, retrieve data from memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  veclocl(*), wbroy(*)
c
c >>>> EXECUTABLE CODE:
c
      ipvec = 1 + nveclocl*(ibl-1)
      call DCOPY( nveclocl, veclocl,1, wbroy(ipvec),1 )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RETRVVEC
c
c
      entry      RETRVVEC( ibl, nveclocl, veclocl, wbroy )
c
      ipvec = 1 + nveclocl*(ibl-1)
      call DCOPY( nveclocl, wbroy(ipvec),1, veclocl,1 )
c
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routines that handle the broyden history data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c       
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2OPEN
c       
c     
      subroutine BROY2OPEN( itersm1,iterb,itmaxbr, sufx,
     $ ibroyfl, nveclocl )
c---------------------------------------------------------------
c Purpose: set up storage of broyden history: u() and dft()
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
c
c Module data for file part of broyden data storage:
      PARAMETER  ( maxfl_broy=20 )
      COMMON /BROYFSTUFF/ nfl_broy,nit_file,ifl_broy,it_file,
     $                    ibrfls(maxfl_broy)
c 
c Module data for memory part of broyden data storage:
      COMMON /BROYMSTUFF/ i0bl0, ibl0, i0bl2, ibl2, nbl2, itertotl
c
      CHARACTER  sufx*(*)
c
c >>>> EXECCUTABLE CODE:
c     
      if( nfl_broy .gt. 0 )then
c       Set up disk file storage:
c
c       Compute nfl_broy=number of files necessary to hold broyden data:
c       msize_it = how big array needed to save per iteration
c       itersav=itersm1 - saves *all* past iters in more and more files
c       itersav=iterb   - uses (and reuses) minimum possible files
c
        msize_it = 2*nveclocl
        itertot = itersm1
        itersav = iterb
c
        call FLGETIWR( IWR )
        call BROYF( IWR, itertot,itersav, itmaxbr, msize_it,
     $   nit_file,maxfl_broy, nfl_broy )
        call BROYFOPEN( sufx, nfl_broy, ibrfls )
        ibroyfl = ibrfls(1)
c
      else
c       Set number of files to zero
        nfl_broy = 0
        ibrfls(1) = -1
        ibroyfl = -1
c       Memory-based storage, starts in third entry, nbl2= history length
        ibl2 = 1
c
      endif
c
      RETURN
c
c       
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2CLOSE
c       
c     
      entry BROY2CLOSE( ibroyfl )
c---------------------------------------------------------------
c Purpose: close storage of broyden history: u() and dft()
c Written: Peter A. Schultz,  9-February-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibroyfl .gt. 0 )then
c       Close the disk files:
        call BROYFCLOSE( nfl_broy,ibrfls )
      else
c       Nothing to do for in-core Broyden, just "rewind"
        ibl2 = 1
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2RWND
c
c
      subroutine BROY2RWND( ibroyfl )
c---------------------------------------------------------------
c Purpose: set storage of vin/dv broyden data to beginning
c Written: Peter A. Schultz, 26-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
c
c Module data for file part of broyden data storage:
      PARAMETER  ( maxfl_broy=20 )
      COMMON /BROYFSTUFF/ nfl_broy,nit_file,ifl_broy,it_file,
     $                    ibrfls(maxfl_broy)
c
c Module data for memory part of broyden data storage:
      COMMON /BROYMSTUFF/ i0bl0, ibl0, i0bl2, ibl2, nbl2, itertotl
c
      DIMENSION  ulocl(*), dlocl(*), blsav(*)
c
c >>>> EXECCUTABLE CODE:
c
      if( ibroyfl .gt. 0 )then
c       Rewind the disks ...
        call BROYF0( ibroyfl, ifl_broy,it_file, ibrfls )
        REWIND( unit=ibroyfl )
      else
c       Reset memory pointer to first entry - u[it=1]
        ibl2 = 1
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2SKIP
c
c
      entry BROY2SKIP( itprune, ibroyfl, nveclocl )
c---------------------------------------------------------------
c Purpose: save u[it] and dft[it] to disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibroyfl .gt. 0 )then
c       Skip through files ...
        if( itprune .gt. 0 )then
c         Skip to beginning of active Broyden history:
          it_skip = itprune
          call BROYFSKIP( itprune,it_skip, ibroyfl,
     $     ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
c 
          REWIND( ibroyfl )
          do  it=1,it_skip
            call READSKP( ibroyfl, nveclocl )
            call READSKP( ibroyfl, nveclocl )
          enddo
        endif
c
      else
c       Skip through memory entries ...
        do  it=1,itprune
          ibl = i0bl2 + ibl2
          if( ibl .gt. nbl2 )then
            ibl2 = 1
            ibl = i0bl2 + ibl2
          endif
          ibl2 = ibl2 + 2
        enddo
c 
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2PUT
c
c
      entry BROY2PUT( ibroyfl, nveclocl, ulocl,dlocl, blsav )
c---------------------------------------------------------------
c Purpose: save u[it] and dft[it] to disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibroyfl .gt. 0 )then
c       write to disk ...
        call BROYFIT( ibroyfl,
     $   ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
        call WRITBIG( ibroyfl, nveclocl, ulocl )
        call WRITBIG( ibroyfl, nveclocl, dlocl )
      else
c       save to memory ...
        ibl = i0bl2 + ibl2
        if( ibl .gt. nbl2 )then
c         Off end of storage, reset index
          ibl2 = 1
          ibl = i0bl2 + ibl2
        endif
        iblu = ibl
        ibld = ibl + 1
        call STASHVEC( iblu, nveclocl, ulocl, blsav )
        call STASHVEC( ibld, nveclocl, dlocl, blsav )
        ibl2 = ibl2 + 2
       endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2GETU
c
c
      entry BROY2GETU( ibroyfl, nveclocl, ulocl, blsav )
c---------------------------------------------------------------
c Purpose: retrieve u[it] from disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibroyfl .gt. 0 )then
c       Get u[it] from disk, and skip dft[it] ...
        call BROYFIT( ibroyfl,
     $   ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
c     
        call READBIG( ibroyfl, nveclocl, ulocl )
        call READSKP( ibroyfl, nveclocl )
      else
c       retrieve from  memory ...
        ibl = i0bl2 + ibl2
        if( ibl .gt. nbl2 )then
c         Off end of storage, reset index
          ibl2 = 1
          ibl = i0bl2 + ibl2
        endif
        iblu = ibl
        call RETRVVEC( iblu, nveclocl, ulocl, blsav )
        ibl2 = ibl2 + 2
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BROY2GETD
c
c
      entry BROY2GETD( ibroyfl, nveclocl, dlocl, blsav )
c---------------------------------------------------------------
c Purpose: retrieve dft[it] from disk/memory
c Written: Peter A. Schultz, 30-January-2007
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( ibroyfl .gt. 0 )then
c       Get dft[it] from disk, skipping u[it]
        call BROYFIT( ibroyfl,
     $   ifl_broy,it_file, nit_file,nfl_broy,ibrfls )
c
        call READSKP( ibroyfl, nveclocl )
        call READBIG( ibroyfl, nveclocl, dlocl )
      else
c       retrieve from  memory ...
        ibl = i0bl2 + ibl2
        if( ibl .gt. nbl2 )then
          ibl2 = 1
          ibl = i0bl2 + ibl2
        endif
        ibld = ibl + 1
        call RETRVVEC( ibld, nveclocl, dlocl, blsav )
        ibl2 = ibl2 + 2
      endif
c
      RETURN
      END
