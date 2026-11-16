c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SCHROED
c
c
      subroutine SCHROED( do_eigvecs,
     $ idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $ norb,nstate,nk, eigval, rwork, seval, scut,
     $ vmat, wmat, iwork )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c  10Jan13-PAS/2.64: pass do_eigvecs and cosmentics
c  25Jun07-PAS/2.60: merge tp and serial
c  11Sep05-PAS/2.59: remove lin.dep.-based dynamic memory use
c  14Aug05-PAS/2.59: d-Ham from disk + reduce vmat mem req.
c  08Jun05-RPM/2.59: added ability to do a canonical orthogonalization
c                    to remove linear dependence in the basis set
c   8Aug03-PAS/2.56: cosmetics; k-vector cleanup
c  24Jul01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract writes
c  23Aug99-PAS/2.38: install zhegvx expert driver eigensolve
c  10May99-PAS/2.35: eigvals cut from ivecfl output
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   5Sep98-PAS/2.25: print out spectral limits
c  14Jul98-PAS/2.23: option to use unchol'd delta-Ham
c  12Aug94-PAS/2.11: cleanup
c---------------------------------------------------------------
c
c When this routine is done
c   (1) it writes out eigenvectors to ivecfl
c   (2) it returns eigenvalues in "eigval()"
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input flags:
      LOGICAL    do_eigvecs
c Output: eigenvalue spectrum
      DIMENSION  eigval(nstate,nk)
c Scratch arrays:
      DIMENSION  rwork(7*norb)
      DIMENSION  seval(norb)
      DIMENSION  iwork(norb,*)
      DIMENSION  vmat(norb*norb)
      DIMENSION  wmat(norb*norb,*)
c
c Local declarations:
      CHARACTER  jobz*1,uplo*1
      COMPLEX*16 cone,czero
      DOUBLE PRECISION  small,scut
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
c
c Do eigenvectors with jobz='V' (jobz="N" turns off eigenvectors):
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Non-scf, compute eigenvalues only
        jobz = 'N'
        if( idosolv.eq.3 )
     $   call STOPXERR( 'idosolv=3 incompatible with non-scf' )
      endif
      uplo = 'L'
c
      REWIND( unit=iham0fl )
      REWIND( unit=iovlpfl )
c
c    Bloch vector loop:
c
      do 1000 k=1,nk
c
c       Assemble Hamiltonian(spin,k)
c
c       Retrieve delta-Hamiltonian from idhamfl
        call READBIG( idhamfl, mat, vmat )
c       Retrieve reference Hamiltonian from iham0fl
        call READBIG( iham0fl, mat, wmat(1,1) )
c       Add del-Ham to ref-Ham to get total Ham
        do  ij=1,mat
          vmat(ij) = wmat(ij,1) + vmat(ij)
        enddo
c
c       Retrieve overlap matrix from iovlpfl
        call READBIG( iovlpfl , mat, wmat(1,2) )
c
c       Load overlap into lower triangle full complex matrix wmat(3:4)
        call C16LOW( norb, wmat(1,2), wmat(1,3) )
c       Load Ham into lower triangle of full complex matrix wmat(1:2)
        call C16LOW( norb, vmat     , wmat(1,1) )
c
c       Need 3*norb-2 for r*8 workspace "rwork"
c       lwork is (complex) length for workspace "l"
        lwzwork = MAX( 4*norb, mat )
        lwork = lwzwork / 2
C        lwork = mat / 2
c
c       Ham is in wmat(,1:2), overlap in wmat(,3:4)
c       Eigenvectors will be returned in wmat(,1:2)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Solve Schroedinger equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if( idosolv .eq. 1 )then
c
c         ***** Full Hermitian eigenproblem *****
c
          call ZHEGV( 1,jobz, 'L',
     $     norb, wmat(1,1),norb,         wmat(1,3),norb,
c -->            H(2,n,n)-i/vec(2,n,n)-o S(2,n,n)-i
     $     eigval(1,k), vmat     ,lwork, rwork, INFO )
c -->                   work(2,lwork)    (3n-2)
c
        elseif( idosolv .eq. 2 )then
c
c         *****  Selected eigenvalue Hermitian eigenproblem *****
c
          call ZHEGVX( 1,jobz, 'I', 'L',
     $     norb, wmat(1,1),norb,         wmat(1,3),norb,
c -->            H(2,n,n)-i              S(2,n,n)-i
     $     -100d0,-20d0, 1,nstate, 2*DLAMCH('S'), nstout,eigval(1,k),
     $     wmat(1,5),norb,    vmat     ,lwork, rwork,
c -->      eigvec(2,n,nstate) work(2,lwork)    (7n)
     $     iwork(1,2),iwork, INFO )
c -->      int(5n)    int(n)
c
c         Put eigenvectors in the right place:
          if( do_eigvecs )then
            call DCOPY( 2*norb*nstate, wmat(1,5),1, wmat(1,1),1 )
          endif
c
        elseif( idosolv .eq. 3 )then
c
c         ***** Deal with linear dependence issues *****
c
c         Parameters that should be soft-coded
c         small - cutoff for filtering small vals (HARDWIRE)
c         scut - cutoff for eigenvalues of S (from input)
          small = 1d-12
c          scut = 1d-6
c
c         Form eigenvectors/vals of S, and from these form
c         orthogonal transformation matrix X. Overwrites wmat(1,3:4)
          cone = DCMPLX( 1d0,0d0 )
          czero = DCMPLX( 0d0,0d0 )
          call FILTERZMAT( wmat(1,3),norb,norb,small )
          call ZHEEV( 'V','L',norb,wmat(1,3),norb,
     $         seval,vmat     ,lwork,rwork,info )
c          print*, 'lowest few S evs'
c          write(*,'(8f10.4)') (seval(i),i=1,8)
          call CANON_XMAT( wmat(1,3),norb,norb,neff,seval,scut,small )
c
c         Square out the upper triangle of H
          call SYMZMAT( wmat(1,1),norb,norb,small )
          call ZGEMM( 'N','N',norb,neff,norb,cone,wmat(1,1),norb,
     $         wmat(1,3),norb,czero,wmat(1,5),norb )
          call ZGEMM( 'C','N',neff,neff,norb,cone,wmat(1,3),norb,
     $         wmat(1,5),norb,czero,wmat(1,1),norb )
c
c          Now solve the transformed matrix
c          print*, 'H',k,nk
c          call printmat( norb,wmat(1,1),3,3 )
           call ZHEEV( 'V','L',neff,wmat(1,1),norb,
     $         eigval(1,k),vmat     ,lwork,rwork,info )
c          print*, 'cprime'
c          call printmat(norb,wmat(1,1),3,3)
c          print*, "eigvals "
c          write(*,'(8f10.4)') (eigval(j,k),j=1,norb)
c
c         The evals are correct, but we need to mult the evecs by
c         X to get the right values
          call DCOPY( 2*norb*norb, wmat(1,1),1, wmat(1,5),1 )
          call ZGEMM( 'N','N',norb,neff,neff,cone,wmat(1,3),norb,
     $         wmat(1,5),norb,czero,wmat(1,1),norb )
c          print*, 'c'
c          call printmat( norb,wmat(1,1),3,3 )
c
        else
c
          call STOPXERR( 'solvopt /invalid SCHROED eigensolver option' )
c
        endif
c
        if( INFO .ne. 0 )then
          call STOPXERR( 'diag-err/ SCHROED diagonalization error' )
        endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Record eigenvalues and eigenvectors
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if( do_eigvecs )then
c         Save eigenfunction coefficients:
          call WRITBIG( ivecfl, 2*norb*nstate, wmat(1,1) )
        endif
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c===============================================================
c    Utility routines
c===============================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CANON_XMAT
c
      subroutine CANON_XMAT( u,lda,n,neff,seval,scut,small )
c
c     Form the canonical orthogonalization transformation 
c     matrix from the eigenvalues (seval) and eigenvectors
c     (u) of the overlap matrix.
c
c     We also need to insure that there aren't any very small
c     values that result from this
c
c     Also returns neff, the number of effective states, i.e.
c     the states with a stable S eigenvalue
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      INTEGER  lda,n,i,j,neff
      COMPLEX*16  u(lda,n)
      DOUBLE PRECISION  seval(n),scut,small
c
      neff = 0
      do  j=1,n
        if( seval(j) .gt. scut )then
          neff = neff + 1
          do  i=1,n
            u(i,neff) = u(i,j) / DSQRT( seval(j) )
          enddo
        endif
      enddo
      call FILTERZMAT( u,lda,n,small )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMZMAT
c
      subroutine SYMZMAT( a,lda,n,small )
c
c     symmetrize the upper triangle of complex*16 matrix a
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      INTEGER  lda,n,i,j
      COMPLEX*16  a(lda,n)
      DOUBLE PRECISION  small
c
      call FILTERZMAT( a,lda,n,small )
      do  i=1,n
        do  j=1,i
          a(j,i) = DCONJG( a(i,j) )
        enddo
      enddo
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERZMAT
c
      subroutine FILTERZMAT( a,lda,n,small )
c
c     Set all values in complex matrix *a* below *small* to zero
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero = 0.d0 )
c
      INTEGER  lda,n,i,j
      COMPLEX*16  a(lda,n)
      DOUBLE PRECISION  a_re, a_im, small
c
      do  j=1,n
        do  i=1,n
          a_re = DBLE( a(i,j) )
          a_im = DIMAG( a(i,j) )
          if( DABS( a_re ) .lt. small ) a_re = zero
          if( DABS( a_im ) .lt. small ) a_im = zero
          a(i,j) = DCMPLX( a_re, a_im )
        enddo
      enddo
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000
characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
