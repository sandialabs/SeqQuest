c 
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPTOOLS
c  
c
c=======================================================================
c Provides wrappers and simple tools for the mpi code
c
c  MPBARRIER    - an MPI barrier
c  MPSENDR4     - send an r4 message
c  MPRECVR4     - recv an r4 message
c  MPSENDR8     - send an r8 message
c  MPRECVR8     - recv an r8 message
c  MPSENDI      - send an i4 message
c  MPRECVI      - recv an i4 message
c  MPSENDTICKET - send a ticket
c  MPRECVTICKET - recv a ticket
c  MPBCASTI     - bcast an integer vector
c  MPBCAST8     - bcast an r8 vector
c  MPBCAST4     - bcast an r4 vector
c  MPREDUCI     - reduce(sum) an integer vector
c  MPREDUC8     - reduce(sum) an r8 vector
c  MPREDUC4     - reduce(sum) an r4 vector
c  MPSPLITR8    - split vector from master over processors
c  MPMERGER8    - merge vector from processors to master
c  MPSPLITGCOMM  - split global communicator into local communicators
c  MPSPLITICOMM  - split global communicator into local communicators
c  MPFREECOMMK   - free MPCOMM_K communicator
c
c=======================================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSYNCH
c
c
      subroutine MPBARRIER( icomm )
c---------------------------------------------------------------
c Purpose: synchronize processors
c---------------------------------------------------------------
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      INTEGER  icomm, ierror
cmpi
      call MPI_Barrier( icomm, ierror  )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR4
c
c
      subroutine MPSENDR4( node_to, lenmsg, r4msg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of R*4 data.
c
c Written: Peter A. Schultz, 11-February-2007, for 2.60
c          Adapted from r8 version
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in R*4 units
c  r4msg()   - the R*4 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      REAL       r4msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( r4msg, lenmsg,
     $     MPI_REAL            , node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR4
c
      subroutine MPRECVR4( node_from, lenmsg, r4msg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      REAL       r4msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c >>>> EXECUTABLE CODE:
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( r4msg, lenmsg,
     $     MPI_REAL            , node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, lenmsg,
     $ ierror )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR8
c
c
      subroutine MPSENDR8( node_to, lenmsg, r8msg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of R*8 data.
c
c Written: Peter A. Schultz, 21-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c  14Oct04-PAS/2.59-tp: add return msg length in recv
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in R*8 units
c  r8msg()   - the R*8 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( r8msg, lenmsg,
     $     MPI_DOUBLE_PRECISION, node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR8
c
      subroutine MPRECVR8( node_from, lenmsg, r8msg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( r8msg, lenmsg,
     $     MPI_DOUBLE_PRECISION, node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, lenmsg,
     $ ierror )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDI
c
c
      subroutine MPSENDI( node_to, lenmsg, imsg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of integer data.
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in integer units
c  imsg()    - the integer data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( imsg, lenmsg,
     $     MPI_INTEGER, node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVI
c
      subroutine MPRECVI( node_from, lenmsg, imsg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( imsg, lenmsg,
     $     MPI_INTEGER, node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_INTEGER, lenmsg,
     $ ierror )
c
      RETURN
      END
c 
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDTICKET
c  
c
      subroutine MPSENDTICKET( itarget, icomm )
c---------------------------------------------------------------
c Purpose: send/recv a ticket to a processor
c Written: Peter A. Schultz, 6-July-2006, for v2.60
c---------------------------------------------------------------
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (lenrcpt=1)
c
      lenmsg = lenrcpt
      rcpt = 1.d0
c
      call MPSENDR8( itarget, lenmsg , rcpt, icomm )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVTICKET
c  
      subroutine MPRECVTICKET( isource, icomm )
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (lenrcpt=1)
c
      lenmsg = lenrcpt
      call MPRECVR8( isource, lenmsg , rcpt, icomm )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCASTI
c
c
      subroutine MPBCASTI( node_from, lenmsg, imsg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for integer broadcast
c
c Written: Andrew C Pineda, 12-November-2009, for v2.63
c
c Revision history:
c  12Nov09-ACP/2.63: adapted from MPBCASTR8
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in integer units
c  imsg()    - the integer data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
c
cmpi
      call MPI_Bcast( imsg, lenmsg, MPI_INTEGER,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCASTI failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST8
c
c
      subroutine MPBCAST8( node_from, lenmsg, r8msg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r8 broadcast
c
c Written: Peter A. Schultz, 21-January-2002, for v2.59-tp
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in R*8 units
c  r8msg()   - the R*8 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
c
cmpi
      call MPI_Bcast( r8msg, lenmsg, MPI_DOUBLE_PRECISION,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCAST8 failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST4
c
c
      subroutine MPBCAST4( node_from, lenmsg, r4msg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r4 broadcast
c
c Written: Peter A. Schultz, 11-January-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in R*8 units
c  r4msg()   - the r*4 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL     r4msg
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r4msg(*)
c
cmpi
      call MPI_Bcast( r4msg, lenmsg, MPI_REAL,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCAST4 failure' )
c
      RETURN
      END
c     
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUCI
c
c
      subroutine MPREDUCI( node_to, nvec, ivecloc, ivecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for integer reduction/sum
c
c Written: Andrew C. Pineda,  23-November-2009, for v2.63
c
c Revision history:
c   none
c---------------------------------------------------------------
c       Use standard mpi reduction as the default scheme
c
      include  'mpif.h'
      itckt = 1
      lenmsg = 1
cmpi
      call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     &  node_to, icomm, ierr )
      call MPI_Reduce( ivecloc, ivecred, nvec, MPI_INTEGER,
     &  MPI_SUM, node_to, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPREDUCI failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC8
c
c
      subroutine MPREDUC8( node_to, nvec, vecloc, vecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r8 reduction/sum
c
c Written: Peter A. Schultz, 30-January-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_to   - target of reduction
c  nvec      - length of vector, in R*8 units
c  vecloc()  - the local R*8 data being sent
c  vecred()  - the reduced vector
c  icomm     - the communicator
c
c Choices for reduction (iopt):
c   default = native mpi reduce/sum
c         1 = reduce direct to target with tickets
c         2 = simple binary tree reduction
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      PARAMETER  ( iopt = 2 )
      PARAMETER  ( r8one = 1.d0 )
c Message space:
      DIMENSION  vecloc(*), vecred(*)
c
c >>>> EXECUTABLE CODE:
c
C      call MPNODES( nprocs )
C      call MPNODE( iproc )
c Use native rather than wrapper (do not know if k-parallel!)
cmpi
      call MPI_comm_rank( icomm, iproc, ierror )
      call MPI_comm_size( icomm, nprocs, ierror )
c
      if( nprocs .eq. 1 )then
c       Nothing to reduce, reduction is simply the local vector
        call DCOPY( nvec, vecloc,1, vecred,1 )
        RETURN
      endif
c
      if( iopt .eq. 1 )then
c       Assumimg reduction of long vec over modest procs, try to reduce
c       mpi contention by reduction directly to target using ticketing:
c
        if( iproc .eq. node_to )then
c         Target proc sends tickets, accumulates reduced vector
c
c         Start off with target's data in reduced vector:
          call DCOPY( nvec, vecloc,1, vecred,1 )
          do 100 ip=1,nprocs
            jproc = ip - 1
            if( node_to .ne. jproc )then
c             Send ticket to proc, get its data, accumulate vector
              call MPSENDTICKET( jproc, icomm )
              nvecrecv = nvec
              call MPRECVR8( jproc, nvecrecv, vecloc, icomm )
              call DAXPY( nvec, r8one, vecloc,1, vecred,1 )
            endif
  100     continue
c
        else
c         Wait for ticket, then send vector to target:
          call MPRECVTICKET( node_to, icomm )
          call MPSENDR8( node_to, nvec, vecloc, icomm )
        endif
c
      elseif( iopt .eq. 2 )then
c       Try a simple by-hand binary tree reduce, to node 0
c
c       Start off copy of local vector in reduced vector:
        call DCOPY( nvec, vecloc,1, vecred,1 )
c
        iprdel = 1
  200   continue
        if( iprdel .ge. nprocs ) goto 260
        nprset = 2*iprdel
c       Cycle over sending process:
        do 250 ipr=iprdel,nprocs-1,nprset
          iprecv = ipr - iprdel
          ipsend = ipr
          if( iproc .eq. iprecv )then
c           Send ticket to source, get its data, accumulate vector
            call MPSENDTICKET( ipsend, icomm )
            nvecrecv = nvec
            call MPRECVR8( ipsend, nvecrecv, vecloc, icomm )
            call DAXPY( nvec, r8one, vecloc,1, vecred,1 )
          elseif( iproc .eq. ipsend )then
c           Wait for ticket from target, then send vector to target:
            call MPRECVTICKET( iprecv, icomm )
            call MPSENDR8( iprecv, nvec, vecred, icomm )
c           This proc is done, get out ...
            goto 260
          endif
  250   continue
        iprdel = 2*iprdel
        goto 200
  260   continue
c
        if( node_to .ne. 0 )then
c         Send final result to right node and right location
          if( iproc .eq. 0 )then
            call MPSENDR8( node_to, nvec, vecred, icomm )
          elseif( iproc .eq. node_to )then
            nvecrecv = nvec
            call MPRECVR8( 0     , nvecrecv, vecred, icomm )
          endif
        endif
c
      else
c       Use standard mpi reduction as the default scheme
c
        itckt = 1
        lenmsg = 1
cmpi
        call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     $   node_to, icomm, ierr )
        call MPI_Reduce( vecloc, vecred, nvec, MPI_DOUBLE_PRECISION,
     $   MPI_SUM, node_to, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc3' )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC4
c
c
      subroutine MPREDUC4( node_to, nvec, vecloc, vecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r4 reduction/sum
c
c Written: Peter A. Schultz, 11-February-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_to   - target of reduction
c  nvec      - length of vector, in R*8 units
c  vecloc()  - the local R*4 data being sent
c  vecred()  - the reduced vector
c  icomm     - the communicator
c
c Choices for reduction (iopt):
c   default = native mpi reduce/sum
c         1 = reduce direct to target with tickets
c         2 = simple binary tree reduction
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      PARAMETER  ( iopt = 2 )
      REAL       r4one
      PARAMETER  ( r4one = 1.e0 )
c Message space:
      REAL       vecloc, vecred
      DIMENSION  vecloc(*), vecred(*)
c
c >>>> EXECUTABLE CODE:
c
C      call MPNODES( nprocs )
C      call MPNODE( iproc )
c Use native rather than wrapper (do not know if k-parallel!)
cmpi
      call MPI_comm_rank( icomm, iproc, ierror )
      call MPI_comm_size( icomm, nprocs, ierror )
c
      if( nprocs .eq. 1 )then
c       Nothing to reduce, reduction is simply the local vector
        call SCOPY( nvec, vecloc,1, vecred,1 )
        RETURN
      endif
c
      if( iopt .eq. 1 )then
c       Assumimg reduction of long vec over modest procs, try to reduce
c       mpi contention by reduction directly to target using ticketing:
c
        if( iproc .eq. node_to )then
c         Target proc sends tickets, accumulates reduced vector
c
c         Start off with target's data in reduced vector:
          call SCOPY( nvec, vecloc,1, vecred,1 )
          do 100 ip=1,nprocs
            jproc = ip - 1
            if( node_to .ne. jproc )then
c             Send ticket to proc, get its data, accumulate vector
              call MPSENDTICKET( jproc, icomm )
              nvecrecv = nvec
              call MPRECVR4( jproc, nvecrecv, vecloc, icomm )
              call SAXPY( nvec, r4one, vecloc,1, vecred,1 )
            endif
  100     continue
c
        else
c         Wait for ticket, then send vector to target:
          call MPRECVTICKET( node_to, icomm )
          call MPSENDR4( node_to, nvec, vecloc, icomm )
        endif
c
      elseif( iopt .eq. 2 )then
c       Try a simple by-hand binary tree reduce, to node 0
c
c       Start off copy of local vector in reduced vector:
        call SCOPY( nvec, vecloc,1, vecred,1 )
c
        iprdel = 1
  200   continue
        if( iprdel .ge. nprocs ) goto 260
        nprset = 2*iprdel
c       Cycle over sending process:
        do 250 ipr=iprdel,nprocs-1,nprset
          iprecv = ipr - iprdel
          ipsend = ipr
          if( iproc .eq. iprecv )then
c           Send ticket to source, get its data, accumulate vector
            call MPSENDTICKET( ipsend, icomm )
            nvecrecv = nvec
            call MPRECVR4( ipsend, nvecrecv, vecloc, icomm )
            call SAXPY( nvec, r4one, vecloc,1, vecred,1 )
          elseif( iproc .eq. ipsend )then
c           Wait for ticket from target, then send vector to target:
            call MPRECVTICKET( iprecv, icomm )
            call MPSENDR4( iprecv, nvec, vecred, icomm )
c           This proc is done, get out ...
            goto 260
          endif
  250   continue
        iprdel = 2*iprdel
        goto 200
  260   continue
c
        if( node_to .ne. 0 )then
c         Send final result to right node and right location
          if( iproc .eq. 0 )then
            call MPSENDR4( node_to, nvec, vecred, icomm )
          elseif( iproc .eq. node_to )then
            nvecrecv = nvec
            call MPRECVR4( 0     , nvecrecv, vecred, icomm )
          endif
        endif
c
      else
c       Use standard mpi reduction as the default scheme
c
        itckt = 1
        lenmsg = 1
cmpi
        call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     $   node_to, icomm, ierr )
        call MPI_Reduce( vecloc, vecred, nvec, MPI_REAL,
     $   MPI_SUM, node_to, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc3' )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITR8
c
c
      subroutine MPSPLITR8( nprocl, iprocl0, iprocl,
     $     nvec,nveclocl, vec, icomm)
c---------------------------------------------------------------
c Purpose: split over a set of processors a long vector existing
c          on the master processor
c Written: Peter A. Schultz, for task-parallel on 2.59
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
c Notes: the scheme used to split the vector to the master node
c   is the simplest possible.  The master sends pieces to each
c   node in sequence.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output:
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      nveclocl = nvec
      if( nprocl .le. 1 ) RETURN
c
c Set up the split:
c  Assume nodes are contiguously sequenced, starting with
c  the source (presumably master) node (iprocl0):
      iprocmx = iprocl0 + nprocl - 1
c  Length of vector to be on first processor:
      nveclocl = nvec / nprocl
      if( nveclocl*nprocl .ne. nvec ) nveclocl = nveclocl + 1
c  Index start of local vector:
      iv1 = 1
c
c Perform the split:
c
      iv2 = iv1 + nveclocl
      nprocleft = nprocl
      nvecleft = nvecleft - nveclocl
      do 100 inode=iprocl0+1,iprocmx
        itarget = inode
        nprocleft = nprocleft - 1
        nvecleft = nvec - iv2 + 1
c       Length of vector segment to be sent:
        nvtarget = nvecleft / nprocleft
        if( nvtarget*nprocleft .ne. nvecleft ) nvtarget = nvtarget + 1
c
        if( iprocl .eq. iprocl0 )then
c         Master sends to the target proc its part of the vector
          call MPSENDR8( itarget, nvtarget, vec(iv2), icomm )
c
        elseif( iprocl .eq. itarget )then
c         Proc receives its part of vector from master
          isource = iprocl0
          nvecrecv = nvtarget
          call MPRECVR8( isource, nvecrecv, vec(iv1), icomm )
          if( nvecrecv .ne. nvtarget ) call STOPXERR( 'MPSPLITR8 err' )
          nveclocl = nvecrecv
          goto 101
        endif
        iv2 = iv2 + nvtarget
c
  100 continue
  101 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPMERGER8
c
c
      subroutine MPMERGER8( nprocl, iprocl0, iprocl,
     $     nvec,ivec1,nveclocl, vec, icomm)
c---------------------------------------------------------------
c Purpose: merge onto master processor a full vector which
c          is broken up in blocks over a set of processors.
c
c Written: Peter A. Schultz, 27-October-2004, for 2.59 task-parallel
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
c Notes: the scheme used to merge the vector to the master node
c   is the simplest possible.  The master sends a "ticket" to
c   each node in sequence, which signals that node to send its
c   part of the vector to the master node.  The master appends
c   the incoming vector onto its existing vector, and proceeds
c   to ticket the next processor until all are done.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output:
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      if( nprocl.le.1 ) RETURN
      lenrcpt = 1
      rcpt = 1.d0
c
c Second, set up the merge:
c  Assume nodes are contiguously sequenced, starting with
c  the target (presumably master) node (iprocl0):
      iprocmx = iprocl0 + nprocl - 1
c  Length of vector (currently) on the local processor:
      nveclen = nveclocl
c  Index of start of existing local vector:
      iv1 = ivec1
c
c Finally, perform the merge to the root node:
c
      if( iprocl .eq. iprocl0 )then
c       Master node, collect from all nodes ...
c
        iv2 = iv1 + nveclen
        do 100 inode=iprocl0+1,iprocmx
          isource = inode
c
c         First, tell the source proc it is ok to send
          call MPSENDR8( isource, lenrcpt, rcpt, icomm )
c
c         Second, receive that proc's part of the vector
          nvecrecv = nvec - iv2 + 1
          call MPRECVR8( isource, nvecrecv, vec(iv2), icomm )
          iv2 = iv2 + nvecrecv
  100   continue
c
      else
c       Non-master, merging full vector to master ...
c
c       First, get ticket from master:
        call MPRECVR8( iprocl0, lenrcpt, rcpt, icomm )
c
c       Send local segment of vector to master:
        call MPSENDR8( iprocl0, nveclen, vec(iv1), icomm )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITGCOMM
c
c
      subroutine MPSPLITGCOMM( mygroup )
c---------------------------------------------------------------
c Purpose: generate local communicators from global communicator
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c >>>> EXECUTABLE CODE:
c
      call MPCOMM_G( icomm_g )
cseq:start
cseq      icomm_l = icomm_g
cseq      nprocs = 1
cseq      node = 0
cseq      node0 = 0
cseq:end
c This sets the master to the 0-node in each group:
      mykey = 0
cmpi:start
      call mpi_comm_split( icomm_g, mygroup, mykey, icomm_l, ierror )
      call mpi_comm_rank( icomm_l, node, ierror )
      call mpi_comm_size( icomm_l, nprocs, ierror )
cmpi:end
c
c Load the comm wrapper (local/image level)
      icommlvl = 2
      node0 = 0
      call MPLOADCOMM( nprocs, node, node0, icomm_l, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITICOMM
c
c
      subroutine MPSPLITICOMM( mygroup )
c---------------------------------------------------------------
c Purpose: generate local communicators from global communicator
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c >>>> EXECUTABLE CODE:
c
      call MPCOMM( icomm_l )
cseq:start
cseq      icomm_k = icomm_l
cseq      nprocs = 1
cseq      node = 0
cseq      node0 = 0
cseq:end
      mykey = 0
cmpi:start
      call mpi_comm_split( icomm_l, mygroup, mykey, icomm_k, ierror )
      call mpi_comm_rank( icomm_k, node, ierror )
      call mpi_comm_size( icomm_k, nprocs, ierror )
cmpi:end
c
c Load the comm wrapper (sub-image, k-parallel level):
      icommlvl = 3
      node0 = 0
      call MPLOADCOMM( nprocs, node, node0, icomm_k, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPFREECOMMK
c
c
      subroutine MPFREECOMMK()
c
      INTEGER istatus
      INCLUDE 'mpif.h'
c
      INTEGER  kparopt, icomm_k
c
      call KPFLAG_GET( kparopt )
      if( kparopt .eq. 2 )then
        call MPCOMM_K( icomm_k )
        call MPI_COMM_FREE( icomm_k, istatus )
      endif
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
