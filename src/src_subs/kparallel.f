c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPARALLEL
c
c
c###############################################################
c Purpose: manage the k-parallel stuff
c
c Written: Peter A. Schultz, 22-November-2008, for 2.63
c
c Revision history:
c  Apr12-ACP&PAS/2.63: add stuff to enable k-parallel in scf
c###############################################################
c Notes on KPARALLEL
c  Processor allocations in each k-point are contiguous ranks,
c  as are the sets of ranks through the k-points
c
c  The strategy implemented here is to do as even a distribution
c  as possible. If perfect balance is not possible, allocate
c  the lighter load to the early processors, so that the image
c  master that sits in the early block is assured of doing
c  no more work (and perhaps less) than any other processor,
c  and thereby hide some of its extra tasks behind others' compute.
c
c
c  We have two possibilities nproc<nk, and nproc>nk.
c  With nproc<nk, if nk is not exactly divisible by nproc,
c    then each later proc gets one assigned one more k-point than
c    each early proc.
c  With nk<nproc, if nproc is not exactly divisible by nk,
c    then each later k-point gets one FEWER proc assigned to it
c    then each earlier k-point (i.e., with more procs/k in the
c    early k-points and procs, they finish more quickly).
c
c Module data: KPSTUFF - manages nature of k-parallelism
c  kparopt - flag that specifies nature of k-parallelims
c  Note: (1) KPSTUFF should probably also record the paramemters that
c            define the k-parallelism (e.g., nkloc, nk0, kpmaster ...)
c        (2) KPSTUFF might be better located folded into MPSTUFF, which
c            already manages the (more coarse) aspects of parallelism
c
c Public routines:
c  KPFLAG_SET - set flags controlling k-parallel operation
c  KPFLAG_GET - get flags controlling k-parallel operation
c  KPBETTER   - evaluate whether k-parallel is good enough
c  KPEFFIC    - computes computational efficiency of k-parallel (ideal)
c  KPDEFINE   - assignments defining k-parallel distribution
c  KPMASTER   - returns rank (local/image) of master of given k-point
c  KPWHERE    - returns range of processor ranks owning given k-point
c  KPSYNC()   - a sync to allow k-masters to catch up to their clients
c
c Public routines for development purposes
c  KPARMERGE - merge all k-hams onto master from k-masters
c  KPARDISTR - distribute k-hams onto k-masters from master
c
c###############################################################
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPFLAG_SET
c
c
      subroutine KPFLAG_SET( kparoptp )
      IMPLICIT NONE
      integer kparoptp
      integer kparopt
c
      COMMON  /KPSTUFF/ kparopt
c
      kparopt = kparoptp
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPFLAG_GET
c
c
      
      subroutine KPFLAG_GET( kparoptp )
      IMPLICIT NONE
      integer kparoptp
      integer kparopt
c
      COMMON  /KPSTUFF/ kparopt
c
      kparoptp = kparopt
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPBETTER
c
c
      subroutine KPBETTER( nk, kparopt )
c---------------------------------------------------------------
c Purpose: unless dictated, decide k-parallel cf. efficiency
c---------------------------------------------------------------
c
      IMPLICIT NONE
c Input:
      INTEGER  nk
c Input/Output:
      INTEGER  kparopt
c Local:
      DOUBLE PRECISION  effic
      INTEGER  nprocs
      INTEGER  idokpar
      DATA     idokpar / -2 /
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
C      call KPFLAG_GET( kparopt )
      if( nprocs.eq.1 .or. nk.eq.1 )then
c       There is no option, must be standard task-parallel
        kparopt = -1
C        call KPFLAG_SET( kparopt )
        goto 900
      endif
      if( kparopt .ne. 0 )then
c       This has been decided/dictated earlier
        goto 900
      endif
c
c Return efficiency of a k-parallel scheme vs. all task-parallel
      call KPEFFIC( effic, nk,nprocs )
c
c Now decide whether k-parallel is likely to outperform
c  a standard all-parallel, based on expected efficiencies.
c
      if( nprocs.le.4 )then
c       On 4 procs, typically get >75% from all-par, ask more from kpar
        if( effic .gt. 0.89d0 ) kparopt = idokpar
c
      elseif( nprocs.le.8 )then
c       On 8 procs, can get as much as 75% from all-par, ask more from kpar
        if( effic .gt. 0.79d0 ) kparopt = idokpar
c
      elseif( nprocs .le. 12 )then
c       On 12 procs, can get as much as 66% from all-par, ask more from kpar
        if( effic .gt. 0.74d0 ) kparopt = idokpar
c
      elseif( nprocs .le. 16 )then
c       On 2^4=16 proc, can get as much as 50% from all-par, ask more from kpar
        if( effic .gt. 0.66d0 ) kparopt = idokpar
c
      else
c       More than 16 processors, compel k-parallel
        kparopt = idokpar
c
      endif
C      call KPFLAG_SET( kparopt )
c
  900 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPEFFIC
c
c
      subroutine KPEFFIC( effickp, nk,nprocs )
c---------------------------------------------------------------
c Purpose: compute ideal efficiency of k-parallel scheme
c---------------------------------------------------------------
c
      IMPLICIT NONE
c Input:
      INTEGER  nk,nprocs
c Output:
      DOUBLE PRECISION  effickp
c Local:
      INTEGER  kpernbig, nperksml
      DOUBLE PRECISION  effic, wkldallp,wkldkpar, speedap,speedkp
c
c >>>> EXECUTABLE CODE:
c
      if( nk .ge. nprocs )then
c       Nprocs < nk, everyone gets their own k-points
c
c       Compare workload (#k/proc) of all-parallel proc vs. k-parallel busiest proc
        kpernbig = ( nk + nprocs - 1 ) / nprocs
        wkldallp = DBLE( nk ) / DBLE( nprocs )
        wkldkpar = DBLE( kpernbig )
c
        effic = wkldallp / wkldkpar
      else
c       Nprocs > nk, multiple nprocs per k-point
        kpernbig = 1
c
c       Compare speed (# procs/k) of all-p kpts vs. k-parallel fewest-node k-point
        nperksml = nprocs / nk
        speedap = DBLE( nprocs ) / DBLE( nk )
        speedkp = DBLE( nperksml )
c
        effic = speedkp / speedap
      endif
c
      effickp = effic
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPSETUP
c
c
      subroutine KPSETUP( nk, kparopt )
c---------------------------------------------------------------
c Purpose: set up k-parallel communicators
c---------------------------------------------------------------
c
c From input parameter ...
c  nk      = number of nk-points for this problem
c  kparopt = flag that determines type of k-parallism
c          = -2 => turn k-parallel ON
c          = -1 => keep k-parallel OFF
c          =  0 => defer choice to code
c
c   ... set the following output parameters:
c  nk0      = 0-index of k-pts assigned to this proc, in global nk
c  nkloc    = number of k-points assigned to this node
c  nklocmax = maximum number of k-points assigned to any processor
c
c   ... and set up the communicators
c
      IMPLICIT NONE
c Input:
      INTEGER  nk, kparopt
c Output:
      INTEGER  nk0,nkloc,nklocmax
c Local:
      INTEGER  nprocs, iproc, icomm
      INTEGER  icolorkp
c
c >>>> EXECUTABLE CODE:
c
c Check that KP has not already been initialized
      if( kparopt .gt. 0 ) RETURN
c
      call KPBETTER( nk, kparopt )
c
      if( kparopt .eq. -2 )then
        call MPNODES( nprocs )
        call MPNODE( iproc )
c       Define processor groups for k-parallel distribution
        call KPDEFINE( nk,nprocs,iproc,
     $   nklocmax,nkloc,nk0, icolorkp )
c
c       Create split communicators along "colors" defined by k-points
        call MPSPLITICOMM( icolorkp )
c
c       Now that kparallel is set up, set kparopt = 2 (k-parallel)
        kparopt = 2
      elseif( kparopt .eq. 0 .or. kparopt .eq. -1 )then
c       No k-parallel, set kparopt = 1 (not k-parallel)
        kparopt = 1
      else
        call STOPXERR( 'KPSETUP: kparopt is wrong' )
      endif
c
      call KPFLAG_SET( kparopt )
c        
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPDEFINE
c
c
      subroutine KPMINE( nk,nkloc,nk0 )
c---------------------------------------------------------------
c Purpose: for this processor, get assignment
c---------------------------------------------------------------
c
c From input parameters ...
c  nk     = number of nk-points for this problem
c  nprocs = total number of processors for this image
c  iproc  = rank of this particular processor in the image
c
c ... determine following output parameters:
c  nkloc    = number of k-points assigned to this node
c  nk0      = 0-index of k-pts assigned to this proc, in global nk
c  
      IMPLICIT NONE
c Input:
      INTEGER  kparopt, nk
c Output:
      INTEGER  nklocmax,nkloc,nk0, icolorkp
c Local:
      INTEGER  nprocs, iproc
      INTEGER  kpernbig,kpernsml, nsmlprocs
      INTEGER  nperkbig,nperksml, lastbigk,lastbign
c
c >>>> EXECUTABLE CODE:
c
c Default is standard all-parallel:
      nk0 = 0
      nkloc = nk
      call KPFLAG_GET( kparopt )
c
      if( kparopt.eq.2 )then
c       Get id of this proc in image
        call MPNODES( nprocs )
        call MPNODE( iproc )
c       Get the k-parallel distribution
        call KPDEFINE( nk,nprocs,iproc,
     $   nklocmax,nkloc,nk0, icolorkp )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPDEFINE
c
c
      subroutine KPDEFINE( nk,nprocs,iproc,
     $ nklocmax,nkloc,nk0, icolorkp )
c---------------------------------------------------------------
c Purpose: define k-parallel distribution
c---------------------------------------------------------------
c
c From input parameters ...
c  nk     = number of nk-points for this problem
c  nprocs = total number ofprocessors for this problem
c  iproc  = rank of this particular processor
c
c ... determine following output parameters:
c  nklocmax = maximum number of k-points assigned to any processor
c  nkloc    = number of k-points assigned to this node
c  nk0      = 0-index of k-pts assigned to this proc, in global nk
c  icolorkp = communicator grouping for this processor
c  
      IMPLICIT NONE
c Input:
      INTEGER  nk, nprocs, iproc
c Output:
      INTEGER  nklocmax,nkloc,nk0, icolorkp
c Local:
      INTEGER  kpernbig,kpernsml, nsmlprocs
      INTEGER  nperkbig,nperksml, lastbigk,lastbign
c
c >>>> EXECUTABLE CODE:
c
      if( nk .ge. nprocs )then
c       Nodes < nk, everyone gets their own k-points
c
c       Every proc in their own group (k-local comm):
        icolorkp = iproc + 1
c
        kpernbig = ( nk + nprocs - 1 ) / nprocs
        kpernsml = kpernbig - 1
        nsmlprocs = nprocs*kpernbig - nk
c
        if( iproc .lt. nsmlprocs )then
c         Proc is small-nk, use 0-index iproc:
          nk0 = kpernsml*iproc
          nkloc = kpernsml
        else
          nk0 = kpernbig*iproc - nsmlprocs*(kpernbig-kpernsml)
          nkloc = kpernbig
        endif
        nklocmax = kpernbig
c
      else
c       Nprocs > nk, every proc only 1 k-point, each kpt maybe multiple procs
c
        nkloc = 1
        nklocmax = 1
        nperkbig = ( nprocs + nk - 1 ) / nk
        nperksml = nperkbig  - 1
        lastbigk = nprocs - nperksml*nk
        lastbign = lastbigk * nperkbig
c
        if( iproc .lt. lastbign )then
c         Proc is in first block of k-points with more procs/k
          nk0 = iproc / nperkbig
        else
c         Proc is in second block of k-points with less procs/k
          nk0 = lastbigk + ( iproc - lastbign ) / nperksml
        endif
c
c       The proc groups (k-local comms) segregated by k-point:
        icolorkp = nk0 + 1
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPMASTER
c
c
      subroutine KPMASTER( nk,nprocs, k,kmasterpr )
c---------------------------------------------------------------
c Purpose: identify parallel distribution
c---------------------------------------------------------------
c
c     nk = number of global k-points
c nprocs = number of (local comm) procs serving this image with nk
c      k = k-point for which the k-local master is desired
c kmasterpr = rank in local comm for k-local master of k-point
c
      IMPLICIT NONE
c Input:
      INTEGER  nk, nprocs
      INTEGER  k
c Output: 
      INTEGER  kmasterpr
c
c Local:
      INTEGER  kparopt
      INTEGER  kproc
      INTEGER  nkloc
      INTEGER  kpernbig,kpernsml, nsmlprocs
      INTEGER  nperkbig,nperksml, lastbigk
c
c >>>> EXECUTABLE CODE:
c
c Note: the following assumes proc ranks (iproc) are 0-indexed
c
      call KPFLAG_GET( kparopt )
c
      if( kparopt.eq.2 )then
c       If k-parallel is ON ...
c
      if( nk .ge. nprocs )then
c       Nodes < nk, everyone got their own k-points
c
        kpernbig = ( nk + nprocs - 1 ) / nprocs
        kpernsml = kpernbig - 1
        nsmlprocs = nprocs*kpernbig - nk
c
        if( k .le. nsmlprocs*kpernsml )then
c         Proc owning this k-point is in small-nkloc early block of procs
          nkloc = kpernsml
          kproc = ( k - 1 ) / kpernsml
        else
c         Proc owning this k-point is in big-nkloc later block of procs
          nkloc = kpernbig
          kproc = ( k + nsmlprocs - 1 ) / kpernbig
        endif
c
      else
c       Nprocs > nk, every proc only 1 k-point, each kpt maybe multiple procs
c
        nkloc = 1
c
        nperkbig = ( nprocs + nk - 1 ) / nk
        nperksml = nperkbig - 1
        lastbigk = nprocs - nperksml*nk
c
        if( k .le. lastbigk )then
c         This k-point is in big-proc early block of k-points
          kproc = ( k - 1 ) * nperkbig
        else
c         This k-point is in small-proc later block of k-points
          kproc = lastbigk*nperkbig  + ( k - lastbigk - 1 ) * nperksml
        endif
c
      endif
c
      kmasterpr = kproc
c
      else
c       If k-parallel is OFF, the image master is k-master
        CALL MPNODE0( kmasterpr )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPWHERE
c
c
      subroutine KPWHERE( nk, ik, iproc_min, iproc_max )
c---------------------------------------------------------------
c Purpose: Find ranks of processors in the current image that
c          hold the k-point ik.
c---------------------------------------------------------------
c
c From input parameters ...
c  nk     = number of nk-points for this problem
c  ik     = ith k-point (1,...,nk)
c Output parameters
c  Determine rank(s) of processes
c    iproc_min  = lowest rank process
c    iproc_max  = highest rank process
c  where iproc_min.le.iproc_max.
c
      IMPLICIT NONE
c
c Input:
      INTEGER  nk, ik
c Output:
      INTEGER  iproc_min, iproc_max
c Local:
c     nprocs = total number of processors in this image
      INTEGER  kparopt
      INTEGER  nprocs, iproc
      INTEGER  kpernbig, kpernsml, nsmlprocs
      INTEGER  nperkbig, nperksml, lastbigk
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call KPFLAG_GET( kparopt )

      if( kparopt.eq.2 )then
c       If k-parallel is ON, then determine where ik is located
c
      if( nk .ge. nprocs )then
c
c       Nodes < nk
c       Oonly one process for each k.
c
        kpernbig = ( nk + nprocs - 1 ) / nprocs
        kpernsml = kpernbig - 1
        nsmlprocs = nprocs*kpernbig - nk
c
        if( ik .lt. nsmlprocs*kpernsml )then
          iproc = (ik-1) / kpernsml
        else
          iproc = nsmlprocs + (ik-nsmlprocs*kpernsml) / kpernbig
        endif
        iproc_min = iproc
        iproc_max = iproc
c
      else
c
c       Nprocs > nk
c       Every proc only 1 k-point, each k-pt might be on multiple procs.
c
c       let a= int(nprocs/nk). b=remainder(nprocs/nk)
c       nperkbig = a if b=0, a+1 if b>0.
        nperkbig = ( nprocs + nk - 1 ) / nk
c       nperksml = a-1 if b=0, a if b>0
        nperksml = nperkbig  - 1
c       lastbigk = nk if b=0, b if b>0
        lastbigk = nprocs - nperksml*nk
c
        if( ik .le. lastbigk )then
c         The lower k have an extra process assigned to them.
          iproc_max = ik*nperkbig - 1
          iproc_min = iproc_max - nperkbig + 1
        else
          iproc_max = lastbigk*nperkbig + (ik-lastbigk)*nperksml - 1
          iproc_min = iproc_max - nperksml + 1
        endif
c
      endif
c
      else
c       If k-parallel if OFF, then ik resides on all processors of image
        call MPNODE0( iproc_min )
        call MPNODES( iproc_max )
        iproc_max = ( iproc_max - 1 ) + iproc_min
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPSYNC
c
c
      subroutine KPSYNC()
c---------------------------------------------------------------
c Purpose: synch clients to their KPMASTER
c
c Written: Peter A. Schultz, 18-April-2012 f0r 2.63 (k-parallel)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER  kparopt
c Local:
      INTEGER  nprocl
      INTEGER  nproc_k, master_k, iproc_k, icomm_k
      INTEGER  r4msg
c
      call MPNODES( nprocl )
      if( nprocl .eq. 1 ) RETURN
c
      call KPFLAG_GET( kparopt )
      if( kparopt .eq. 2 )then
c       Get k-local MP info
        call MPNODES_K( nproc_k )
        if( nproc_k .gt. 1 )then
          call MPNODE0_K( master_k )
          call MPNODE_K( iproc_k ) 
          call MPCOMM_K( icomm_k ) 
c
c         Use a bcast from the master to synch clients to k-master
          r4msg = 0.0
          call MPBCAST4( master_k, 1, r4msg, icomm_k )
        endif
      endif
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPARMERGE
c
c
      subroutine KPARMERGE( norb,nk, vijk )
c---------------------------------------------------------------
c Purpose: merge k-parallel Hams from kmasters onto image master
c
c Written: Peter A. Schultz,  4-December-2008, for 2.63 (k-parallel)
c
c Revision history:
c   4Dec08-PAS/2.63: pushed down from vslomat
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input: nk=global k-pts, nkloc=local number of k-pts
      INTEGER    norb, nk,nkloc
c Input/Output
      DOUBLE PRECISION  vijk
      DIMENSION  vijk(norb*norb,*)
c
c Local declarations:
      INTEGER    kparopt
      INTEGER    nprocl,masterl,iprocl,icomml, master,nprocs,iproc
      INTEGER    nk0, k,ikloc,  kmasterpr
      INTEGER    mat, lenmsg
c
c  For k-parallel, send all the k-resolved matrices to image master
c
      call KPFLAG_GET( kparopt )
      if( kparopt .eq. 2 )then
c       Get image/local MP info
        call MPNODES( nprocl )
        call MPNODE0( masterl )
        call MPNODE( iprocl )
        call MPCOMM( icomml )
c       Get k-local MP info
        call MPNODES_K( nprocs )
        call MPNODE0_K( master )
        call MPNODE_K( iproc ) 
c  Set up k-parallel
c        call KPMINE( kparopt, nk,nkloc,nk0 )
        call KPMINE( nk,nkloc,nk0 )
c
        mat = norb*norb
c  
        if( iprocl.eq.masterl )then
c         Image master, collect matrices from k-local masters
          do  k=1,nk
            call KPMASTER( nk,nprocl, k,kmasterpr )
            if( kmasterpr .ne. masterl )then
              call MPSENDTICKET( kmasterpr, icomml )
              lenmsg = mat
              call MPRECVR8( kmasterpr, lenmsg, vijk(1,k), icomml )
            endif
          enddo
        elseif( iproc .eq. master .and. nkloc .gt. 0 )then
c         Master of this k-point, send image/local master the matrix
          do  ikloc=1,nkloc
            call MPRECVTICKET( masterl, icomml )
            lenmsg = mat
            call MPSENDR8( masterl, lenmsg, vijk(1,ikloc), icomml )
          enddo
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPARDISTR
c
c
      subroutine KPARDISTR( norb,nk,vijk )
c---------------------------------------------------------------
c Purpose: merge k-parallel Hams from kmasters onto image master
c
c Written: Peter A. Schultz,  4-December-2008, for 2.63 (k-parallel)
c
c Revision history:
c   4Dec08-PAS/2.63: pushed down from vslomat
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input: nk=global k-pts, nkloc=local number of k-pts
      INTEGER    norb, nk,nkloc
c Input/Output
      DOUBLE PRECISION  vijk
      DIMENSION  vijk(norb*norb,*)
c
c Local declarations:
      INTEGER    kparopt
      INTEGER    nprocl,masterl,iprocl,icomml, master,nprocs,iproc
      INTEGER    nk0, k,ikloc,  kmasterpr
      INTEGER    mat, lenmsg
c
c  For k-parallel, send all the k-resolved matrices to image master
c
      call KPFLAG_GET( kparopt )
      if( kparopt .eq. 2 )then
c       Get image/local MP info
        call MPNODES( nprocl )
        call MPNODE0( masterl )
        call MPNODE( iprocl )
        call MPCOMM( icomml )
c       Get k-local MP info
        call MPNODES_K( nprocs )
        call MPNODE0_K( master )
        call MPNODE_K( iproc ) 
c  Set up k-parallel
        call KPMINE( nk,nkloc,nk0 )
c
        mat = norb*norb
c  
        if( iprocl .eq. masterl )then
c         Image master, collect matrices from k-local masters
          do  k=1,nk
            call KPMASTER( nk,nprocl, k, kmasterpr )
            if( kmasterpr .ne. masterl )then
              call MPRECVTICKET( kmasterpr, icomml )
              lenmsg = mat
              call MPSENDR8( kmasterpr, lenmsg, vijk(1,k), icomml )
            endif
          enddo
        elseif( iproc .eq. master .and. nkloc .gt. 0 )then
c         Master of this k-point, send image/local master the matrix
          do  ikloc=1,nkloc
            call MPSENDTICKET( masterl, icomml )
            lenmsg = mat
            call MPRECVR8( masterl, lenmsg, vijk(1,ikloc), icomml )
          enddo
        endif
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
