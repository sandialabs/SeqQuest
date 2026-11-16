c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHEIGVEC
c
c
      subroutine PRHEIGVEC( ivecfl, ncplx,norb,nstate,numpop,nk,
     $ wtk,eigpop,rtpop, eigvec,eigvrd )
c---------------------------------------------------------------
c Purpose: read eigenvectors and sort them into an order most
c          convenient for use in the density routines, while
c          reducing precision, too.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
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
      REAL  eigvec
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  wtk(nk),eigpop(nstate,nk), rtpop(nstate)
      DIMENSION  eigvec(ncplx,numpop,nk,norb)
      DIMENSION  eigvrd(ncplx,norb,nstate)
c
c >>>> EXECUTABLE CODE:
c
      nvecrd = ncplx*norb*nstate
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      do 100 k=1,nk
c
        if( iproc.eq.master ) call READBIG( ivecfl, nvecrd, eigvrd )
        if( nprocs.gt.1 ) call MPBCAST8( master, nvecrd, eigvrd, icomm )
c
c       Set up occupancies for states (strictly non-negative).
        do  n=1,numpop
          rtpop(n) = wtk(k)*eigpop(n,k)
          if( rtpop(n).gt.zero )then
            rtpop(n) = SQRT( rtpop(n) )
          elseif( rtpop(n).lt.zero )then
            call STOPXERR( 'EF:pop<0 - must use DM density' )
          endif
        enddo
c
        if( ncplx.eq.2 )then
c         Complex eigenvector matrix
          do  i=1,norb
            do  n=1,numpop
              eigvec(1,n,k,i) = rtpop(n)*eigvrd(1,i,n)
              eigvec(2,n,k,i) = rtpop(n)*eigvrd(2,i,n)
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
c    That's all Folks!
c
      RETURN
      END
