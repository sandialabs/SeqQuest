c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSTUFF
c
c
c###############################################################
c Purpose: manage MPI handling in the code
c
c Written: Peter A. Schultz, 7-December-2001, for 2.51
c
c Revision history:
c  19Nov08-PAS/2.62: to accommodate k-parallel communicator layer
c   6Jul06-APT/2.60: to merge task and image parallel
c###############################################################
c
c External dependencies:
c  none
c
c Public routines:
c  Initialize communicator for communicator level icommlvl:
c    MPLOADCOMM( node,node0,nodes,icomm, icommlvl )
c  Global communicator, icommlvl=1:
c    MPNODES_G( nodes ) - returns number of nodes in global communicator
c    MPNODE_G( node )   - returns rank of this node in global communicator
c    MPNODE0_G( node0 ) - returns rank of master in global communicator
c    MPCOMM_G( icomm )  - returns communicator id for global communicator
c  Local communicator (image layer), icommlvl=2:
c    MPSTAMP( IWR )  - print out communicator info for local communicator
c    MPNODES( nodes )   - returns number of nodes in local communicator
c    MPNODE( node )     - returns rank of this node in local communicator
c    MPNODE0( node0 )   - returns rank of master in local communicator
c    MPCOMM( icomm )    - returns communicator id for local communicator
c  "k-local" communicator (sub-image "k-parallel" layer), icommlvl=3:
c    MPNODES_K( nodes ) - returns number of nodes in k-local communicator
c    MPNODE_K( node )   - returns rank of this node in k-local communicator
c    MPNODE0_K( node0 ) - returns rank of master in k-local communicator
c    MPCOMM_K( icomm )  - returns communicator id k-local communicator
c  Note: MPSTUFF data structure is set 4-deep, allowing for another level
c
c Private routines:
c    MPSET{N,I,0,C}( arg, icommlvl ) - sets arg for icommlvl
c    MPGET{N,I,0,C}( arg, icommlvl ) - returns arg for icommlvl
c     ... where arg = nodes, rank of node, master rank, communicator
c         for {N,I,0,C} routines respectively
c
c###############################################################
c Notes on the MPSTUFF module:
c
c  This is not a particularly clean, general implementation,
c   but a minimalist module just to make wrapping the mpi possible.
c  We will refine this module as the parallelism in the code matures,
c   and the specific needs for better parallelism mgmt become clearer.
c  The "parallel" global/local/k-local accessors are notably clunky.
c  The code itself must keep track of a lot of information
c   about processor organization, very awkward. E.g., ...
c  An ad hoc machinery is necessary to keep track of node-id's
c   between levels of communicators, e.g., to collect in upper
c   communicator results from masters of lower communicators.
c  However, this does hide direct mpi from most of the core code.
c
c###############################################################
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSTAMP
c
c
      subroutine MPSTAMP( IWR )
c---------------------------------------------------------------
c Purpose: print out local processor allocation
c
c Written: Aidan P. Thompson, 6-July-2006, for v2.60
c
c Revision history:
c  21Jun07-PAS/2.60: pushed into subroutine
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER  IWR
      INTEGER  nodes,node,node0,icomm
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nodes )
      call MPNODE( node )
      call MPNODE0( node0 )
      call MPCOMM( icomm )
c
      if( IWR .ge. 0 )then
        write(IWR,'(1x,a)')
     $       'MPINFO::: Local Communicator         :::'
        write(IWR,'(1x,a,i4,a)')
     $       'MPINFO::: Local Context = ',icomm  ,'       :::'
        write(IWR,'(1x,a,i4,a)')
     $       'MPINFO::: Local Size =    ',nodes  ,'       :::'
        write(IWR,'(1x,a,i4,a)')
     $       'MPINFO::: Local Root =    ',node0  ,'       :::'
        write(IWR,'(1x,a,i4,a)')
     $       'MPINFO::: Local Rank =    ',node   ,'       :::'
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPLOADCOMM
c
c
      subroutine MPLOADCOMM( nprocs, node, node0, icomm, icommlvl )
c---------------------------------------------------------------
c Purpose: initialize wrapper layers
c
c Written: Peter A. Schultz, 19_November-2008, for 2.62
c
c Revision history:
c   none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
c  Currently max out at 4 levels in MPSTUFF:
      ncommlvl = 4
c  Initialize this and all subordinate communicators:
      do  icommset=icommlvl,ncommlvl
        call MP_SETN( nprocs, icommset )
        call MP_SETI( node  , icommset )
        call MP_SET0( node0 , icommset )
        call MP_SETC( icomm , icommset )
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODES
c
c
      subroutine MPNODES( nprocs )
c---------------------------------------------------------------
c Purpose: return the local number of processors/nodes
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 2
      call MP_GETN( nprocs, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE
c
c
      subroutine MPNODE( node )
c---------------------------------------------------------------
c Purpose: return the local node number of this processor
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 2
      call MP_GETI( node, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE0
c
c
      subroutine MPNODE0( node0 )
c---------------------------------------------------------------
c Purpose: return the local root node "node0"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 2
      call MP_GET0( node0, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPCOMM
c
c
      subroutine MPCOMM( icomm )
c---------------------------------------------------------------
c Purpose: return the local communicator/context
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 2
      call MP_GETC( icomm,icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODES_K
c
c
      subroutine MPNODES_K( nprocs )
c---------------------------------------------------------------
c Purpose: return the local number of processors/nodes
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 3
      call MP_GETN( nprocs, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE_K
c
c
      subroutine MPNODE_K( node )
c---------------------------------------------------------------
c Purpose: return the local node number of this processor
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 3
      call MP_GETI( node, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE0_K
c
c
      subroutine MPNODE0_K( node0 )
c---------------------------------------------------------------
c Purpose: return the k-local root node rank "node0"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 3
      call MP_GET0( node0, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPCOMM_K
c
c
      subroutine MPCOMM_K( icomm )
c---------------------------------------------------------------
c Purpose: return the k-local communicator/context
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 3
      call MP_GETC( icomm,icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODES_G
c
c
      subroutine MPNODES_G( nprocs )
c---------------------------------------------------------------
c Purpose: return the global number of processors/nodes
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 1
      call MP_GETN( nprocs, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE_G
c
c
      subroutine MPNODE_G( node )
c---------------------------------------------------------------
c Purpose: return the global node number of this processor
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 1
      call MP_GETI( node, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPNODE0_G
c
c
      subroutine MPNODE0_G( node0 )
c---------------------------------------------------------------
c Purpose: return the global root node "node0"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 1
      call MP_GET0( node0, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPCOMM_G
c
c
      subroutine MPCOMM_G( icomm )
c---------------------------------------------------------------
c Purpose: return the global communicator/context
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      icommlvl = 1
      call MP_GETC( icomm, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SUBS
c
c
c#######################################################################
c     Internal routines that manipulate MPSTUFF module data
c---------------------------------------------------------------
c Purpose: internal primitive routines to manipulate MPSTUFF data
c
c Written: Peter A. Schultz, 7-December-2001, for v2.51
c
c Revision history:
c  19Nov08-PAS/2.62: unify MPSTUFF data, expand
c   6Jul06-APT/2.60: merged task parallel and image parallel
c#######################################################################
c
      subroutine MP_SUBS( icommlvl )
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SUBS
c Clear to a non-parallel default
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /MPSTUFF/ nodetot(4),nodeid(4),node0id(4),icommid(4)
c
      nodetot(icommlvl) = 1
      nodeid(icommlvl)  = 0
      node0id(icommlvl) = 0
      icommid(icommlvl) = 0
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SETN
c Set the number of nodes
      entry MP_SETN( nprocs, icommlvl )
      nodetot(icommlvl) = nprocs
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SETI
c Set id of current node
      entry MP_SETI( node, icommlvl )
      nodeid(icommlvl) = node
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SET0
c Set node0
      entry MP_SET0( node0, icommlvl )
      node0id(icommlvl) = node0
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_SETC
c Set communicator/context
      entry MP_SETC( icomm, icommlvl )
      icommid(icommlvl) = icomm
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_GETN
c Get the number of nodes
      entry MP_GETN( nprocs, icommlvl )
      nprocs = nodetot(icommlvl)
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_GETI
c Get id of current node
      entry MP_GETI( node, icommlvl )
      node = nodeid(icommlvl)
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_GET0
c Get node0
      entry MP_GET0( node0, icommlvl )
      node0 = node0id(icommlvl)
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MP_GETC
c Get communicator/context
      entry MP_GETC( icomm, icommlvl )
      icomm = icommid(icommlvl)
      RETURN
c
c    That's all Folks!
c
      END
