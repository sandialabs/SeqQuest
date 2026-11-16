c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOIJ
c
c
      subroutine RHOIJ( ivecfl,idmatfl,iematfl, norb,nstate,nk,
     $ eigval,eigpop,fpop,npop, wtk,
     $ dmat,emat, eigvec )
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis
c
c Written: Peter A. Schultz,  9-April-1999, for v2.34
c
c Revision history:
c  31Dec19-PAS/2.68b: only read/use occupied eigvecs (mem reduce)
c  18Jan02-PAS/2.51b: replace code that led to segmentation flt.
c  20Jul01-PAS/2.49: spin-polarized dft, fpop made scratch, dmat
c  10May99-PAS/2.35: eigenfunction grid density
c   9Apr99-PAS/2.34: complex lapack eigensolver implemented
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Oct97-PAS/2.21: change fpop/wtk convention
c  12Aug94-PAS/2.11: cleanup/purge complex/clean args
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  eigval(nstate,nk),eigpop(nstate,nk),wtk(nk)
c scratch arrays:
      DIMENSION  fpop(nstate)
      DIMENSION  dmat(norb,*),emat(norb,norb)
      DIMENSION  eigvec(2,npop  ,norb)
c NB: dmat used temp as 2*nstate*norb for eigvec transpose
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
c
c    Bloch vector loop:
c
      do 100 k=1,nk
c
c Retrieve (only occupied numpop, not all nstate) eigenvectors ...
        call READBIG( ivecfl, 2*norb*npop  , dmat )
c -->                                        2*norb*numpop-o
c NB: dmat used as 2*norb*numpop scratch space here
c
c  ... and transpose for use below:
        call TREIGVEC( norb,npop  , dmat,eigvec )
c
        numpop = npop
c
        do  n=1,numpop
          fpop(n) = eigpop(n,k)
C          do  i=1,norb
C            eigvec(1,n,i) = dmat(1,i,n)
C            eigvec(2,n,i) = dmat(2,i,n)
C          enddo
        enddo
c
c Clear out dmat, emat:
        call MKZERO( mat, dmat )
        call MKZERO( mat, emat )
c
c Done with fpop, scale by weight for use below.
        wt = wtk(k)
        do  n=1,numpop
          fpop(n) = wt*fpop(n)
        enddo
c
c Load one triangle ...
        do 32 i=1,norb
c
c         Do *off*-diagonals first ...
          do 31 j=1,i-1
            dmr = zero
            dmi = zero
            emr = zero
            emi = zero
            do 20 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,j)*eigvec(1,n,i)
     $                   + eigvec(2,n,j)*eigvec(2,n,i) )
              ti = occnk*( eigvec(2,n,j)*eigvec(1,n,i)
     $                   - eigvec(1,n,j)*eigvec(2,n,i) )
              dmr = dmr + tr
              dmi = dmi + ti
              emr = emr + tr*eigval(n,k)
              emi = emi + ti*eigval(n,k)
   20       continue
            dmat(j,i) = dmi
            dmat(i,j) = dmr
            emat(j,i) = emi
            emat(i,j) = emr
   31     continue
c
c          ... and do the diagonal element by itself:
c           (if not split off, compiler might seg-fault)
            dmr = zero
            emr = zero
            do 22 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,i)*eigvec(1,n,i)
     $                   + eigvec(2,n,i)*eigvec(2,n,i) )
              dmr = dmr + tr
              emr = emr + tr*eigval(n,k)
   22       continue
            dmat(i,i) = dmr
            emat(i,i) = emr
c
   32   continue
c
        call WRITBIG( idmatfl, mat, dmat )
        call WRITBIG( iematfl, mat, emat )
c
  100 continue
c
c    That's all Folks!
c
      RETURN
      END
