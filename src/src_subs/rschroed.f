      subroutine RSCHROED(idosolv,scut,ispin,nspin,
     $     idhamfl,iham0fl,iovlpfl,ivecfl,
     $     norb,nstate,eigval,sefp,iwork,
     $     ifail,lwork4,vmat)

c Purpose: solve Schroedinger equation given input ham and ovlp
c
c Written: Rick Muller,  June 4, 2007

      implicit NONE
      integer info,i,j,neff,nstout
      integer idosolv,idhamfl,iham0fl,iovlpfl,ivecfl
      integer ispin,nspin
      integer lwork4,lwork
      integer norb,nstate,mat
      double precision eigval(nstate),scut
      double precision vmat(norb,norb,*)
c     sefp first used to store eigenvalues of S in algorithm
c     to remove linear dependence from basis, and then as
c     occupation factor to set up density matrix construction
      double precision  sefp(nstate)
      integer  iwork(5*norb),ifail(norb)
      CHARACTER  jobz*1
      double precision DLAMCH
      double precision zero,half,one,two
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
c
c  Do eigenvectors with jobz='V' (jobz='N' turns off eigvecs)
      jobz = 'V'
      mat = norb*norb

c     Retrieve delta-Ham from disk if spin (if not, already in mem)
      if( nspin .eq. 2 ) call READBIG( idhamfl, mat, vmat(1,1,1) )
c
c     Retrieve iteration-independent Hamiltonian, overlap matrices:
      REWIND( unit=iham0fl )
      call READBIG( iham0fl, mat, vmat(1,1,2) )
      REWIND( unit=iovlpfl )
      call READBIG( iovlpfl , mat, vmat(1,1,3) )
c
c  Combine iteration-independent Ham and del-Ham to get total Ham:
      call DAXPY( mat, one, vmat(1,1,1),1, vmat(1,1,2),1 )
c
c  Total hamiltonian is now in vmat(,,2), overlap in vmat(,,3)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Solve Schroedinger Equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c lwork is workspace available within vmat(,,4), and normally you
c would set lwork=mat and that would be more than ample.
c But for small problems (e.g. H-atom), we crash as the routine wants
c more than the small mat it gets.  AAAArrrrrgggghhhh!
c Hence I pass in the full space available ...
      lwork = lwork4
c
      if( idosolv.eq.1 )then
c
c     ***** Full real symmetric eigenproblem *****
c
         call DSYGV( 1,jobz, 'L',
     $        norb, vmat(1,1,2),norb,vmat(1,1,3),norb, eigval,
c     -->          H(n,n)-i/vec(n,n)-o S(n,n)-i
     $        vmat(1,1,4),lwork, INFO )
c     -->      work(lwork)
c
c     Put eigenvectors into vmat(,,1)
         call DCOPY( norb*nstate, vmat(1,1,2),1,vmat(1,1,1),1 )
c
      elseif( idosolv.eq.2 )then
c
c     ***** Selected eigvals/vecs real symmetric eigenproblem *****
c
         call DSYGVX( 1,jobz, 'I','L',
     $        norb, vmat(1,1,2),norb, vmat(1,1,3),norb,
c     -->            H(n,n)-i          S(n,n)-i
     $        -100d0,-20d0,1,nstate,2*DLAMCH('S'),
     $        nstout,eigval,
     $        vmat(1,1,1),norb, vmat(1,1,4), lwork,
c     -->      vec(n,nstate)     work(lwork)
     $        iwork,    ifail, INFO )
c     -->      int(5nt)-s in(n)-s
c
      elseif( idosolv .eq. 3 )then
c
c     ***** Deal with linear dependence issues *****
c
c     Rick's new code for the solver that removes 
c     linearly dependent basis sets
c     
c     
c     Compute the eigenvalues/vectors of the overlap matrix
c2         call DSYEV( 'V','L',norb,vmat(1,1,3),norb,
c2     $        sefp ,vmat(1,1,4),lwork,info )
c
c     Use the eigenvectors to do form the orthogonalizing matrix X
c2         neff = 0
c2         do  j=1,norb
c2            if( sefp(j) .gt. scut )then
c2               neff = neff + 1
c2               do  i=1,norb
c2                  vmat(i,neff,3) = vmat(i,j,3) / DSQRT( sefp(j) )
c2               enddo
c2            endif
c2         enddo
         call symorth(norb,vmat(1,1,3),norb,sefp,vmat(1,1,4),
     $        lwork,scut,neff)
c     
c     First square H out, since only the lower triangle is stored
         do  j=1,norb
            do  i=1,j
               vmat(i,j,2) = vmat(j,i,2)
            enddo
         enddo
c
c     Use X for similarity transformation of H, stored in vmat(1,1,2)
c     Overwrite vmat(1,1,2)
c1         call DGEMM( 'N','N',norb,neff,norb,one,vmat(1,1,2),norb,
c1     $        vmat(1,1,3),norb,zero,vmat(1,1,1),norb )
c1         call DGEMM( 'T','N',neff,neff,norb,one,vmat(1,1,3),norb,
c1     $        vmat(1,1,1),norb,zero,vmat(1,1,2),norb )
         call simx(norb,neff,vmat(1,1,2),norb,vmat(1,1,3),norb,
     $        vmat(1,1,1),norb)

c     Now solve the transformed matrix
         call DSYEV( 'V','L',neff,vmat(1,1,2),norb,
     $        eigval,vmat(1,1,4),lwork,INFO )
c     
c     The eigvals are now correct. Transform the eigvecs using X
c     and put them in vmat(1,1,1)
         call DGEMM( 'N','N',norb,neff,neff,one,vmat(1,1,3),norb,
     $        vmat(1,1,2),norb,zero,vmat(1,1,1),norb )
c
      elseif( idosolv .eq. 4 )then
         
         

       else
c
          call STOPXERR( 'solv-opt/ ERLSOLV has invalid solver option' )
c     
       endif
c
       if( INFO .ne. 0 )then
          call STOPXERR( 'diag-err/ ERLSOLV diagonalization error' )
       endif
c
c     Save the the eigenvectors in vmat(,,1) to ivecfl
       call WRITBIG( ivecfl, norb*nstate, vmat(1,1,1) )
c
        
       return
       end

      subroutine simx(n,m,A,lda,B,ldb,AB,ldab)
c     Compute the similarity tranformation A=BtAB
      implicit none
      integer n,m,lda,ldb,ldab
      double precision A(lda,n),B(ldb,n),AB(ldab,n)

      call DGEMM('N','N',n,m,n,1d0,A,lda,B,ldb,0d0,AB,ldab)
      call DGEMM('T','N',m,m,n,1d0,B,ldb,AB,ldab,0d0,A,lda)

      return
      end

      subroutine simxt(n,m,A,lda,B,ldb,AB,ldab)
c     Compute the similarity tranformation A=BABt
      implicit none
      integer n,m,lda,ldb,ldab
      double precision A(lda,n),B(ldb,n),AB(ldab,n)

      call DGEMM('N','T',n,m,n,1d0,A,lda,B,ldb,0d0,AB,ldab)
      call DGEMM('N','N',m,m,n,1d0,B,ldb,AB,ldab,0d0,A,lda)

      return
      end

      subroutine symorth(n,S,lds,seval,work,lwork,scut,neff)
c     Compute the orthogonalizing transform based on S
      implicit none
      integer n,lds,info,i,j,neff,lwork
      double precision S(lds,n),work(lwork),seval(n),scut

      call DSYEV('V','L',n,S,lds,seval,work,lwork,info)
      neff = 0
      do j=1,n
         if (seval(j) .gt. scut) then
            neff = neff + 1
            do i=1,n
               S(i,neff) = S(i,j)/SQRT(seval(j))
            enddo
         endif
      enddo
      return
      end
