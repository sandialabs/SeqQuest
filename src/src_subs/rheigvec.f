c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHEIGVEC
c
c
      subroutine RHEIGVEC( ivecfl, ncplx,norb,nstate,numpop,
     $ nk,nk0,nkloc,
     $ wtk,eigpop,rtpop, eigvec,eigvrd )
c---------------------------------------------------------------
c Purpose: read eigenvectors and sort them into an order most
c          convenient for use in the density routines, while
c          reducing precision, too.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
c  26Nov08-PAS/2.63: k-parallel eigvec distribution
c  24Sep08-PAS/2.62: better diagnostics on occup<0 error
c   6Jul06-APT/    : merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  27Oct04-PAS/2.59: task-parallelism
c---------------------------------------------------------------
c
c Note: have "nstate" eigenvectors that must be read through,
c but only "numpop" of those states are occupied and are needed
c for the density calculation, hence the different dimensions
c of the eigenvectors as read in (eigvrd) and output (eigvec).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL       eigvec
c
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  wtk(nk),eigpop(nstate,nk), rtpop(nstate)
      DIMENSION  eigvec(ncplx,numpop,nkloc,norb)
      DIMENSION  eigvrd(ncplx,norb,nstate)
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES( nprocl )
      call MPNODE( iprocl )
      call MPNODE0( masterl )
      call MPCOMM( icomml )
c
      call MPNODES_K( nproc_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
      nvecrd = ncplx*norb*nstate
c
      ikloc = 0
      do 100 ik=1,nk
c
        lenmsg = nvecrd
        call KPMASTER( nk,nprocl, ik,kmasterpr )
        if( nkloc .eq. nk ) kmasterpr = masterl
c
        if( iprocl .eq. masterl )then
c         Image/local master reads eigvecs, sends to correct k-master proc
          call READBIG( ivecfl, nvecrd, eigvrd )
          if( kmasterpr .ne. masterl )then
c           This k-point goes off to kmasterpr
            call MPSENDR8( kmasterpr, lenmsg, eigvrd, icomml )
            goto 100
          endif
        elseif( iprocl .eq. kmasterpr )then
c         Proc is master of this k-point, get eigvecrd from masterl
          call MPRECVR8( masterl, lenmsg, eigvrd, icomml )
        else
c         This proc is not leading this k-point
          goto 100
        endif
        ikloc = ikloc + 1
c
c       Set up occupancies for states (must be strictly non-negative).
        do  n=1,numpop
          rtpop(n) = wtk(ik)*eigpop(n,ik)
          if( rtpop(n).gt.zero )then
            rtpop(n) = SQRT( rtpop(n) )
          elseif( rtpop(n).lt.zero )then
c           We have an error - do some documentation
            call FLGETIERR( IERR )
            write(IERR,*) 'RHEIGVEC: state occups mangled, k=',ik
            write(IERR,'(1x,1p8d10.3)') (eigpop(i,k),i=1,numpop)
            write(IERR,*) 'If pop<0 intended, must use DM density'
            call STOPXERR( 'RHEIGVEC: negative state populations' )
          endif
        enddo
c
        if( ncplx.eq.2 )then
c         Complex eigenvector matrix
          do  i=1,norb
            do  n=1,numpop
              eigvec(1,n,ikloc,i) = rtpop(n)*eigvrd(1,i,n)
              eigvec(2,n,ikloc,i) = rtpop(n)*eigvrd(2,i,n)
            enddo
          enddo
        else
c         Real eigenvector matrix
          do  i=1,norb
            do  n=1,numpop
              eigvec(1,n,1,i) = rtpop(n)*eigvrd(1,i,n)
            enddo
          enddo
        endif
c
  100 continue
c
C      if( nkloc. ne. nk .and. nkloc .gt.0 .and. nproc_k .gt. 1 )then
      if( nproc_k .gt. 1 .and. nkloc .gt. 0 )then
c       Send around each k-group its r4 eigvec from kmaster
        lenmsg = ncplx*nkloc*norb*numpop
        call MPBCAST4( master_k, lenmsg, eigvec, icomm_k )
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
