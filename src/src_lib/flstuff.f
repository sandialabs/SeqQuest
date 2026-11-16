c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLSTUFF
c
c
c#######################################################################
c Module FLSTUFF
c
c Purpose: take care of file/unit/name handling for fortran code
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  29Dec19-PAS/2.68a: reconcile FLNSTUFF/FLSPECS commons
c   6Jul06-APT/    : merged task parallel and image parallel
c   7Aug03-PAS/2.55a: cosmetics, and Cplant file synch workaround.
c  29Oct02-PAS/2.54: cleanup, document, and debug flag
c---------------------------------------------------------------
c
c The machinery in this package constructs fully qualified path names:
c     <filesp> = <path>/<filenm> = <path>/<jobnm>.<filetag>
c  The variable name in the calling arguments denotes what pieces of
c  this are assembled in the routine to give the ultimate <filesp>.
c
c The user does not specify unit numbers (except for standard in/out/err
c  in the initialization).  This package tracks available unit numbers,
c  allocates a unit number as an *output* when it opens a file, and
c  makes that unit number available again when the file is closed.  The
c  code then identifies all open units when the module is terminated
c  and closes all outstanding open files.
c
c External dependencies:
c  STOPXERR( msg ) which outputs <msg> and STOPs code when fatal error
c  MPNODE_G( node ) which returns global rank of node for debugging output
c  STRPARS( string,n1,n2,i1,i2 ) - first field (i1:i2) in string(n1:n2)
c
c >>>>> Module initialization (invoked before first file access):
c     subroutine FLINIT( IRD, IWR, IERRFL )
c        IRD,IWR,IERRFL(input) = set unit #'s for standard in/out/err
c        Initializes all unit numbers as unused.
c
c >>>>> Module termination (invoked when code completed):
c     subroutine FLFINI()
c        Closes all outstanding open file units.
c
c >>>>> File manipulation:
c  (NB: input character variables flstat and flform are shorthand for
c  fortran status and format keywords, and are translated in FL_CHKEYS)
c  Check for existence of file with name "filenm"
c     subroutine FLEXIST( iexist, filenm, flform )
c  Open ASCII file with tag "filetag" (input), and unit# ifl (output)
c     subroutine FLOPENA( ifl, filetag )
c  Open binary file with tag "filetag" (input), and unit# ifl (output)
c     subroutine FLOPENB( ifl, filetag )
c  Open file with name "filenm" (input), and unit# ifl (output)
c     subroutine FLOPENF( ifl, filenm, flstat, flform )
c  Open file with name "filenm", using specified unit# ifl (**input**)
c     entry FLOPENU( ifl, filenm, flstat, flform )
c  Open file with full spec "filesp" (input), and unit# ifl (output)
c     subroutine FLOPEN( ifl, filesp, flstat, flform, ierror )
c  Translation of fortran open keyword (flform/flstat):
c     subroutine FL_CHKEYS( flkey, flkeyword, ierror )
c  Close file with unit# ifl (input)
c     subroutine FLCLOSE( ifl )
c  Close and delete file with unit# ifl (input)
c     subroutine FLDELETE( ifl )
c
c >>>>> File name/path building routines
c  Extend a path name with nadd-length string "addpath", i.e.
c  <path> = <path> // <addpath> // '/' (i.e. note closing slash)
c     subroutine FLADDPATH( nadd, addpath )
c  Set job name "jobnm" (i.e. file names made as <jobnm>.<filetag> )
c     subroutine FLSETJOB( jobnm )
c  Take input string "filetag", output full file spec "filesp"
c     subroutine FLNAMESP( filetag, filesp )
c  Take input string "filetag", output filename "filenm"
c     subroutine FLNAME( filetag, filenm )
c  Take input filename "filenm", output file spec "filesp"
c     subroutine FLSPEC( filenm, filesp )
c
c >>>>> Unit number manipulation routines
c   subroutine FLGETU( ifl ) - return unused unit# ifl (input), taken
c   entry FLFINDU( ifl ) - returns unused unit# ifl (output), reserved
c   entry FLGETIRD( ifl ) - returns standard in unit# ifl (output)
c   entry FLGETIWR( ifl ) - returns standard out unit# ifl (output)
c   entry FLGETIERR( ifl ) - returns standard error unit# ifl (output)
c
c >>>>> Internal routines:
c   FLN_SUBS - primitives for path/filename manipulations
c   FLU_SUBS - primitives for unit number manipulation
c#######################################################################
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLINIT
c
c
      subroutine FLINIT( IRD, IWR, IERRFL )
c---------------------------------------------------------------
c Purpose: initialize file handling module
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  29Dec19-PAS/2.68a: get rid of unneeded&misnamed FLSPECS common
c---------------------------------------------------------------
c
c  IRD    = standard input unit number
c  IWR    = standard output unit number
c  IERRFL = standard error output unit number
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
c Set up units as all unused:
      call FLU_INIT
c Set up null directory structure:
      call FLN_INIT
c
c Set the standard i/o units, and reserve the unit numbers:
      call FLU_IRDSET( IRD )
      call FLU_RSV( IRD, ierror )
      call FLU_IWRSET( IWR )
      call FLU_RSV( IWR, ierror )
      call FLU_IERRSET( IERRFL )
      if( IERRFL .ne. IWR ) call FLU_RSV( IERRFL, ierror )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLFINI
c
c
      subroutine FLFINI
c---------------------------------------------------------------
c Purpose: close all open files
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( nfld = 100 )
      COMMON  /FLUSTUFF/ IRDU,IWRU,IERRU,IFLU1,IFLU2,iflstat(nfld)
c
c >>>> EXECUTABLE CODE:
c
c Check all unit numbers for open files, and close them:
      do  jfl=1,nfld
        if( iflstat(jfl) .eq. 2 )then
          call FLCLOSE( jfl )
          iflstat(jfl) = 0
        endif
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLEXIST
c
c
      subroutine FLEXIST( iexist, filenm, flform )
c---------------------------------------------------------------
c Purpose: check whether a file exists
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Returns iexist=0 if file does not exist
c         iexist=1 if file does exist
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: full file spec;  Output: format of file to be tested
      CHARACTER*(*)  filenm, flform
c Local: format of file to be tested
      CHARACTER*(12) form
      CHARACTER  filesp*128
c
c >>>> EXECUTABLE CODE:
c
      call FLU_NUM( ifl, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLEXIST1: no unit# avail' )
c
      call FL_CHKEYS( flform, form, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLEXIST2: bad format key' )
c
      call FLSPEC( filenm, filesp )
      OPEN( unit=ifl,file=filesp,status='OLD',form=form,err=1301 )
      CLOSE( unit=ifl )
c
      iexist = 1
      RETURN
 1301 iexist = 0
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLOPENA
c
c
      subroutine FLOPENA( ifl, filetag )
c---------------------------------------------------------------
c Purpose: open a standard ascii file
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: tag for file
      CHARACTER*(*) filetag
c Local: full file spec for file
      CHARACTER  filenm*128
c
c >>>> EXECUTABLE CODE:
c
      call FLNAME( filetag, filenm )
      call FLOPENF( ifl, filenm, 'UNK', 'FOR' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLOPENB
c
c
      subroutine FLOPENB( ifl, filetag )
c---------------------------------------------------------------
c Purpose: open a standard binary file
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: tag for file to be opened
      CHARACTER*(*) filetag
c Local: full file spec for file
      CHARACTER  filenm*128
c
c >>>> EXECUTABLE CODE:
c
      call FLNAME( filetag, filenm )
      call FLOPENF( ifl, filenm, 'UNK', 'UNF' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLOPENF
c
c
      subroutine FLOPENF( ifl, filenm, flstat, flform )
c---------------------------------------------------------------
c Purpose: open a file (in need of a path)
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c   6Feb02-PAS/2.51c: move some stuff around within FL.
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: full file spec, status of file, format of file
      CHARACTER*(*)  filenm, flstat,flform
c Local: (fortran) status of file, format of file
      CHARACTER*(12) stat,form
      CHARACTER      filesp*128
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Find a free unit number:
      call FLU_NUM( ifl, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPENF1: no unit# avail' )
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLOPENU
c
      entry FLOPENU( ifl, filenm, flstat, flform )
c
c Check that unit number is actually available:
      call FLU_CHK( ifl, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPENF2: unit# already open' )
c
c Set up the keywords for the fortran open statement:
      call FL_CHKEYS( flstat, stat, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPENF3: bad status key' )
      call FL_CHKEYS( flform, form, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPENF4: bad format key' )
c
c Open the file ...
      call FLSPEC( filenm, filesp )
      if( idebug_fl .gt. 0 )then
        call MPNODE_G( node )
        call STRPARS( filesp, 1,80, i1,i2 )
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPENF/node=',node,
     $   ', attempt unit=',ifl,   ', filesp=',filesp(i1:i2)
        write(*, '(a,i4, a,i2, a,a,a)' )  'FLOPENF/node=',node,
     $   ', attempt unit=',ifl,   ',         stat,form=',stat,form
      endif
c
cxxxcplant:start * * * * * * * * * * * * * * * * * * * * * * * * *
      if( flstat.eq.'OLD' )then
c       This section of code prompted by file system synch problem
c       on Cplant: code may not recognize existence of input file
c       on *first* attempt to open file.  Preliminary touch with
c       an INQUIRE (that might lie and say not-exist) is enough to
c       cause synch and have the open work. - 07Aug03-PAS
c        The ACCESS was Cplant suggestion - but fails on new OS.
c        if( flstat.eq.'OLD') call ACCESS( filesp, 0 )
        if( flstat.eq.'OLD' ) call FLEXIST( iexist, filenm, flform )
        if( iexist .eq. 0 )then
          call MPNODE_G( node )
          call STRPARS( filesp, 1,128, i1,i2 )
          write(*,'(a,i4,2x,a)') '***** OLD file missing on node',
     $     node,filesp(i1:i2)
        endif
      endif
cxxxcplant:end * * * * * * * * * * * * * * * * * * * * * * * * * *
      OPEN( unit=ifl,file=filesp,status=stat,form=form,err=1301 )
c  ... and record the unit number as being with an open file:
      call FLU_OPN( ifl, ierror )
c
      if( idebug_fl .gt. 0 )then
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPENF/node=',node,
     $   ', open ok unit=',ifl,   ', filesp=',filesp(i1:i2)
      endif
c
      RETURN
c
 1301 continue
      if( idebug_fl .gt. 0 )then
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPENF/node=',node,
     $   ', open FAILED unit=',ifl,   ', filesp=',filesp(i1:i2)
      endif
      call FLU_IERR( IWR )
      write(IWR,*) '***** ERROR *****'
      write(IWR,*) 'File open failed: unit number=',ifl
      write(IWR,*) '  file=',filesp
      call STOPXERR( 'FLOPENF: file open failed' )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLOPEN
c
c
      subroutine FLOPEN( ifl, filesp, flstat, flform, ierror )
c---------------------------------------------------------------
c Purpose: open a file with full file specification
c
c Written: Peter A. Schultz, 6-February-2002, for v2.51(c)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c ifl = output:  unit number for the opened file
c filesp = input: fully qualified path name for a file
c flstat = input: shorthand for fortran open "status"
c flform = input: shorthand for fortran open "format"
c ierror = output: 0 if successful open, otherwise non-zero
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: full file spec, status of file, format of file
      CHARACTER*(*)  filesp, flstat,flform
c Local: (fortran) status of file, format of file
      CHARACTER*(12) stat,form
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
      ierror = 0
c Find a free unit number:
      call FLU_NUM( ifl, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPEN1: no unit number avail' )
c
c Set up the keywords for the fortran open statement:
      call FL_CHKEYS( flstat, stat, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPEN3: bad status key' )
      call FL_CHKEYS( flform, form, ierror )
      if( ierror.ne.0 ) call STOPXERR( 'FLOPEN4: bad format key' )
c
c Open the file ...
      if( idebug_fl .gt. 0 )then
        call MPNODE_G( node )
        call STRPARS( filesp, 1,80, i1,i2 )
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPEN/node=',node,
     $   ', attempt unit=',ifl,   ', filesp=',filesp(i1:i2)
        write(*, '(a,i4, a,i2, a,a,a)' )  'FLOPEN/node=',node,
     $   ', attempt unit=',ifl,   ',         stat,form=',stat,form
      endif
cxxxcplant: could prep OLD open with touch, but not enfs here.
      OPEN( unit=ifl,file=filesp,status=stat,form=form,err=1301 )
c  ... and record the unit number as being with an open file:
      call FLU_OPN( ifl, ierror )
c
      if( idebug_fl .gt. 0 )then
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPEN/node=',node,
     $   ', open ok unit=',ifl,   ', filesp=',filesp(i1:i2)
      endif
c
      RETURN
c
 1301 ierror = 1
      if( idebug_fl .gt. 0 )then
        write(*, '(a,i4, a,i2, a,a)' )  'FLOPEN/node=',node,
     $   ', open FAILED unit=',ifl,   ', filesp=',filesp(i1:i2)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FL_CHKEYS
c
c
      subroutine FL_CHKEYS( flkey, flkeyword, ierror )
c---------------------------------------------------------------
c Purpose: find full correct keyword for fortran open statement
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Input: internal key for file open
      CHARACTER*(*)   flkey
c Output: full fortran keyword for open statement
      CHARACTER*(12)  flkeyword
c Local:
      CHARACTER*(3)   keyfl
c
c >>>> EXECUTABLE CODE:
c
      ierror = 0
      keyfl = flkey(1:3)
c
      if( keyfl .eq. 'UNK' )then
        flkeyword = 'unknown'
      elseif( keyfl .eq. 'NEW' )then
        flkeyword = 'new'
      elseif( keyfl .eq. 'OLD' )then
        flkeyword = 'old'
c
      elseif( keyfl .eq. 'UNF' )then
        flkeyword = 'unformatted'
      elseif( keyfl .eq. 'FOR' )then
        flkeyword = 'formatted'
c
      elseif( keyfl .eq. 'KEE' )then
        flkeyword = 'keep'
      elseif( keyfl .eq. 'DEL' )then
        flkeyword = 'delete'
c
      else
        flkeyword = ' '
        ierror = 1
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLCLOSE
c
c
      subroutine FLCLOSE( ifl )
c---------------------------------------------------------------
c Purpose: close (and keep) a file with unit number "ifl"
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Local:
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Close the file ...
      if( idebug_fl .gt.0 )then
        call MPNODE_G( node )
        write(*, '(a,i4, a,i2)' )  'FLCLOSE/node=',node,
     $   ', closed unit=',ifl
      endif
      CLOSE( unit=ifl, status ='KEEP' )
c  ... and record the unit number as being available again:
      call FLU_CLS( ifl, ierror )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLDELETE
c
c
      subroutine FLDELETE( ifl )
c---------------------------------------------------------------
c Purpose: close (and delete) a file with unit number "ifl"
c
c Written: Peter A. Schultz, 9-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Local:
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Close the file ...
      if( idebug_fl .gt. 0 )then
        call MPNODE_G( node )
          write(*, '(a,i4, a,i2)' )  'FLDELETE/node=',node,
     $     ', closed unit=',ifl
      endif
      CLOSE( unit=ifl, status ='DELETE' )
c  ... and record the unit number as being available again:
      call FLU_CLS( ifl, ierror )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLSETJOB
c
c
      subroutine FLSETJOB( jobnm )
c---------------------------------------------------------------
c Purpose: set a job name for the code
c
c Written: Peter A. Schultz, 7-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: job name
      CHARACTER*(*)  jobnm
c Local:
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
      call FLN_SETJ( jobnm )
c
      if( idebug_fl .gt. 0 )then
        call MPNODE_G( node )
        call STRPARS( jobnm, 1,80, i1,i2 )
        write(*, '(a,i4, a,a)' )  'FLSETJOB/node=', node,
     $   ', jobname=',jobnm(i1:i2)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLADDPATH
c
c
      subroutine FLADDPATH( nadd, addpath )
c---------------------------------------------------------------
c Purpose: add to path name
c
c Written: Peter A. Schultz, 7-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: addition to directory path for files
      CHARACTER*(*)  addpath
c Local:
      DATA  idebug_fl / 0 /
c
c >>>> EXECUTABLE CODE:
c
      call FLN_ADDP( nadd, addpath )
c
      if( idebug_fl .gt. 0 )then
        call MPNODE_G( node )
        write(*, '(a,i4, a,a)' )  'FLADDPATH/node=',node,
     $   ', addpath=',addpath(1:nadd)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLNAMESP
c
c
      subroutine FLNAMESP( filetag, filesp )
c---------------------------------------------------------------
c Purpose: construct full path name from file tag (add dir/job)
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: file tag;  Output: full file spec
      CHARACTER*(*)  filetag, filesp
c Local:
      CHARACTER*128  filenm
c
c >>>> EXECUTABLE CODE:
c
      call FLNAME( filetag, filenm )
      call FLSPEC( filenm, filesp )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLNAME
c
c
      subroutine FLNAME( filetag, filenm )
c---------------------------------------------------------------
c Purpose: connect job name to file tag to get file name
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revi  sion history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /FLNSTUFF/ ndirspec, njobspec, dirspec, jobspec
      CHARACTER  dirspec*128, jobspec*128
c
c Input: file tag;  Output: file name
      CHARACTER*(*)  filetag, filenm
c Local:
      CHARACTER*(1)  chdot
      DATA           chdot / '.' /
c
c >>>> EXECUTABLE CODE:
c
      filenm = jobspec(1:njobspec) // chdot // filetag
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLSPEC
c
c
      subroutine FLSPEC( filenm, filesp )
c---------------------------------------------------------------
c Purpose: connect file name to directory to get full file spec
c
c Written: Peter A. Schultz, 6-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /FLNSTUFF/ ndirspec, njobspec, dirspec, jobspec
      CHARACTER  dirspec*128, jobspec*128
c
c Input: file name;  Output: full file spec
      CHARACTER*(*)  filenm, filesp
c
c Local:
      CHARACTER*(128)  chtmp
c
c >>>> EXECUTABLE CODE:
c
      chtmp = filenm
      call STRPARS( chtmp, 1,128, n1,nd )
c
      if( n1 .gt. 0 )then
        filesp = dirspec(1:ndirspec) // filenm(n1:nd)
      else
        filesp = dirspec(1:ndirspec)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLGETU
c
c
      subroutine FLGETU( ifl )
c---------------------------------------------------------------
c Purpose: find a free i/o unit number "ifl"
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
c Get a free unit number ...
      call FLU_NUM( jfl, ierror )
c  ... and mark it as opened:
      call FLU_OPN( jfl, ierror )
      ifl = jfl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLFINDU
c
      entry FLFINDU( ifl )
c
c Get a free unit number ...
      call FLU_NUM( jfl, ierror )
c  ... and reserve it:
      call FLU_RSV( jfl, ierror )
      ifl = jfl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLGETIRD
c
      entry FLGETIRD( ifl )
c
c Get the standard input unit number:
      call FLU_IRD( jfl )
      ifl = jfl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLGETIWR
c
      entry FLGETIWR( ifl )
c
c Get the standard output unit number:
      call FLU_IWR( jfl )
      ifl = jfl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLGETIERR
c
      entry FLGETIERR( ifl )
c
c Get the standard output error unit number:
      call FLU_IERR( jfl )
      ifl = jfl
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLN_SUBS
c
c
c---------------------------------------------------------------
c Purpose: internal primitive routines to handle file/path names
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLN_INIT
c
c
      subroutine FLN_INIT
c---------------------------------------------------------------
c Purpose: initialized directory info to null
c
c Written: Peter A. Schultz, 6-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /FLNSTUFF/ ndirspec, njobspec, dirspec, jobspec
      CHARACTER  dirspec*128, jobspec*128
c
c >>>> EXECUTABLE CODE:
c
c Blank directory name ...
      dirspec = ' '
c  ... with zero length:
      ndirspec = 0
c Standard job name ...
      call FLN_SETJ( 'job' )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLN_ADDP
c
c
      subroutine FLN_ADDP( nch, chfl )
c---------------------------------------------------------------
c Purpose: extend a path name
c
c Written: Peter A. Schultz, 6-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /FLNSTUFF/ ndirspec, njobspec, dirspec, jobspec
      CHARACTER  dirspec*128, jobspec*128
c
      CHARACTER*(*)  chfl
c Local:
      CHARACTER*128  chtmp
      CHARACTER*1    slash
      DATA           slash / '/' /
c
c >>>> EXECUTABLE CODE:
c
      chtmp = chfl
      call STRPARS( chtmp, 1,nch, n1,n2 )
      if( n1 .eq. 0 ) RETURN
c
      len2 = n2 - n1 + 1
      len1 = ndirspec
      nd1 = len1 + 1
      nd2 = len1 + len2 + 1
      if( nd2.gt.128 ) call STOPXERR( 'dir spec too long >128' )
c
c Augment directory specification:
      ndirspec = nd2
      dirspec(nd1:nd2) = chfl(n1:n2) // slash
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLN_SETJ
c
c
      subroutine FLN_SETJ( chfl )
c---------------------------------------------------------------
c Purpose: overwrite a job name
c
c Written: Peter A. Schultz, 6-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /FLNSTUFF/ ndirspec, njobspec, dirspec, jobspec
      CHARACTER  dirspec*128, jobspec*128
c
      CHARACTER*(*)  chfl
c Local:
      CHARACTER*(128)  chtmp
c
c >>>> EXECUTABLE CODE:
c
      chtmp = chfl
      call STRPARS( chtmp, 1,128, n1,n2 )
      if( n1 .eq. 0 ) RETURN
c
      len2 = n2 - n1 + 1
      len1 = 0
      nd1 = len1 + 1
      nd2 = len1 + len2
c
      njobspec = nd2
      jobspec(nd1:nd2) = chfl(n1:n2)
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_SUBS
c
c
c---------------------------------------------------------------
c Purpose: internal primitive routines to handle files/units
c
c Written: Peter A. Schultz, 5-Decemeber-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c In the common block FLUSTUFF:
c   IRDU,IWRU,IERRU = standard in/out/err unit numbers
c   IFLU1,IFLU2 (not used) = two other standard units
c   ifstat(ifl) - status of unit number: 0=free,1=reserved,2=opened
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_INIT
c Initialize all unit numbers as free:
      subroutine FLU_INIT
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( nfld = 100 )
      COMMON  /FLUSTUFF/ IRDU,IWRU,IERRU,IFLU1,IFLU2,iflstat(nfld)
c
      do  jfl=1,nfld
        iflstat(jfl) = 0
      enddo
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IRDSET
c Set standard input unit
      entry FLU_IRDSET( ifl )
      IRDU = ifl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IWRSET
c Set standard output unit number
      entry FLU_IWRSET( ifl )
      IWRU = ifl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IERRSET
c Set standard error output unit number
      entry FLU_IERRSET( ifl )
      IERRU = ifl
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IRD
c Return standard input unit number
      entry FLU_IRD( ifl )
      ifl = IRDU
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IWR
c Return standard output unit number
      entry FLU_IWR( ifl )
      ifl = IWRU
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_IERR
c Return standard error output unit number
      entry FLU_IERR( ifl )
      ifl = IERRU
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_NUM
c Return available unit number
      entry FLU_NUM( ifl, ierror )
      ierror = 0
c I intentionally leave the single digit units untouched
      do  jfl=10,nfld
        if( iflstat(jfl) .eq. 0 )then
c         This unit number is available
          ifl = jfl
          RETURN
        endif
      enddo
c     We did not find a free unit number, very bad:
      ifl = -1
      ierror = 5
      call STOPXERR( 'FLU_NUM: no free unit numbers for file' )
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_CHK
c Check if file is already open
      entry FLU_CHK( ifl, ierror )
      ierror = 0
c     Check if unit number attached to opened file:
      if( iflstat(ifl) .eq. 2 ) ierror = 1
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_RSV
c Set unit number as reserved
      entry FLU_RSV( ifl, ierror )
      if( iflstat(ifl) .eq. 0 )then
        iflstat(ifl) = 1
        ierror = 0
      else
        ierror = 1
      endif
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_OPN
c Set unit number as in-use and opened, and unavailable
      entry FLU_OPN( ifl, ierror )
      if( iflstat(ifl) .ne. 2 )then
        ierror = 0
      else
        ierror = 1
      endif
      iflstat(ifl) = 2
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FLU_CLS
c Set unit number as closed file, making it available
      entry FLU_CLS( ifl, ierror )
      if( iflstat(ifl) .eq. 2 )then
        ierror = 0
      else
        ierror = 1
      endif
      iflstat(ifl) = 0
      RETURN
c
c    That's all Folks!
c
      END
