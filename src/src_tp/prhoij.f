c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOSTRIP
c
c
      subroutine RHOSTRIP( jstripe, norb, nprocs )
c---------------------------------------------------------------
c Purpose: width of mat-stripe (larger or "equitable" widths)
c Written: Peter A. Schultz, 31Dec19, for 2.68b
c---------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  norb,nprocs,    jstripe
c
      jstripe = (norb+nprocs-1) / nprocs
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOIJ
c
c
      subroutine PRHOIJ( ivecfl,idmatfl,iematfl, norb,nstate,nk,
     $ eigval,eigpop,fpop,npop, wtk,
     $ dmat,emat, eigvec )
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis, parallel
c
c Written: Peter A. Schultz,  9-April-1999, for v2.34
c
c Revision history:
c  29Mar20-PAS/2.68e: READB1 to read partial ivecfl(npop<nstate)
c  31Dec19-PAS/2.68b: reduced-mem eigvec(npop)
c  11May07-PAS/2.60: merged tp and serial
c  12Jul06-PAS/2.60: stripe-distributed dmat generation
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  DMC: task parallel
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
c mem/cpu optimization notes:
c  MUST use transpose despite (2*norb*npop) mem cost
c    - stride in dm-gen loops leads to 4-5 X cost increase
c  pre-weight eigvec(n;) wtih sqrt(wtk*pop) in transpose not helpful
c    - eliminates 1 (second-order) op in dm-gen, but only negligible gain
c    - code risk, might want neg wtk in future for some unforeseen purpose
c Possible future ideas:
c   skip emat in scf, only do emat after scf (complicated, extra solve)
c   block or block-cyclic dm-mat gen (complicated code)
c Now, this mem cost below core solvers, cpu cost is marginal, leave alone
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  eigval(nstate,nk),eigpop(nstate,nk),wtk(nk)
c Scratch arrays:
      DIMENSION  fpop(nstate)
c
c     *** NB: dmat used as 2*norb*npop   scratch below, in addition
c         to norb^2 output.  Currently, the idle emat is contiguous
c         to dmat, so this is safe, even for npop  =norb.
      DIMENSION  dmat(norb,*),emat(norb,norb)
      DIMENSION  eigvec(2,npop  ,norb)
c
c Local declarations:
      DATA       lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      call FLGETIWR( IWR )
      if( lstout.gt.0 .and. nprocs.gt.1 .and. iproc.eq.masterp )then
        write(IWR,*) 'DMAT: all-parallel stripe/merge dmat(k)'
        write(IWR,*) 'PRHOIJ(06): compact stripe+transpose'
      endif
c
      j1 = 1
      jn = norb
      jlen = norb
c Take care of arithmetic for stripe mgmt, using (1):
c   (1) most equitable, rather than
c   (2) blocked (i.e. BLACS compatible), or (3) 0-favored
c  All assign mem allocations of wider stripe (everyone has same larger size)
c  Thinner work allocations (shorter width by 1) will go to first processors
c  First, determine width of (wider) stripe in matrix
C      jstripe = (norb+nprocs-1) / nprocs ! replaced with call (so
C      peigsolv can use
      call RHOSTRIP( jstripe, norb, nprocs )
c  Allocate memory locations (not used - idea for combined mem later)
C      j0dmat = 0
C      j0emat = j0dmat + jstripe
C      j0vec  = j0emat + jstripe
c
      jlen = jstripe - 1
      nshort = ( jstripe * nprocs ) - norb
      j1 = 1 + iproc*jlen
      if( iproc.ge.nshort )then
c       Last nlong procs get bigger stripe (by 1), assume 0:(nprocs-1) id
        jlen = jlen + 1
        j1 = j1 + iproc - nshort
      endif
      jn = j1 - 1 + jlen
c
c    Bloch vector loop:
c
      do 1000 k=1,nk
c
c       Retrieve (only occupied npop, not all nstate) eigenvectors ...
c       NB: dmat used as 2*norb*npop   scratch space here
        if( iproc .eq. masterp )then
c         A read of the partial vector, npop out of nstate saved
          call READB1 ( ivecfl, 2*norb*nstate, 2*norb*npop, dmat )
        endif
        call MPBCAST8( masterp, 2*norb*npop  , dmat  , icomm )
c
c       Every processor works on every k-point
c
c       Transpose (only occupied!) eigvecs for use below:
        call TREIGVEC( norb,npop  , dmat,eigvec )
        numpop = npop
c
c       Scale fpop by k-point sampling weight for use in dmat loops
        wt = wtk(k)
        do  n=1,numpop
          fpop(n) = wt*eigpop(n,k)
        enddo 
c       
c       Clear out dmat, emat stripes
        nmatstr = jlen*norb
        call MKZERO( nmatstr, dmat )
        call MKZERO( nmatstr, emat )
c       
        jdm = 0
        jem = 0
        do 400 j=j1,jn
          jdm = jdm + 1
          jem = jem + 1
c
c         Do imaginary part first ...
          do 210 i=1,j-1
            dmi = zero
            emi = zero
            do 200 n=1,numpop 
              occnk = fpop(n)
              ti = occnk*( eigvec(2,n,i)*eigvec(1,n,j)
     $                   - eigvec(1,n,i)*eigvec(2,n,j) )
              dmi = dmi + ti 
              emi = emi + ti*eigval(n,k) 
  200       continue
            dmat(i,jdm) = dmi 
            emat(i,jem) = emi
  210     continue
c           
c          ... and then switch to real part:
          do 310 i=j,norb
            dmr = zero
            emr = zero
            do 300 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,i)*eigvec(1,n,j)
     $                   + eigvec(2,n,i)*eigvec(2,n,j) )
              dmr = dmr + tr
              emr = emr + tr*eigval(n,k)
  300       continue
            dmat(i,jdm) = dmr
            emat(i,jem) = emr
  310     continue
c
c         Close loop over j/norb
  400   continue
c
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat compute' )
c
c >>>>> Store dmat and emat away
c
c Here, dmat&emat in distributed stripes and we are done with eigvecs()
c This version, we merge to masterp and write to disk to save
c
c       Must do emat() first, as it is in the way of the full dmat
        if( nprocs.gt.1 )then
c         Merge emat to masterp
          nvec = mat
          nmatstr = jlen*norb
c         Now stripe frontloaded locally, so  ...
c          ivec1 = (j1-1)*norb + 1
          ivec1 = 1
c
          call MPMERGER8( nprocs, masterp,iproc,
     $                    nvec,ivec1,nmatstr, emat,
     $                    icomm )
          if( lstout.gt.3 ) call TIMER( 'prhoij/emat merged' )
        endif
c
c       Save full emat to iematfl for use in the force section
        if( iproc .eq. masterp )then
          call WRITBIG( iematfl, mat, emat )
          if( lstout.gt.3 ) call TIMER( 'prhoij/emat writ' )
        endif
c       And now done with the emat(), so its space is free, switch to dmat ...
c
        if( nprocs.gt.1 )then
c         Merge dmat to masterp
          nvec = mat
          nmatstr = jlen*norb
c         Now stripe frontloaded locally, so  ...
c          ivec1 = (j1-1)*norb + 1
          ivec1 = 1
c
          call MPMERGER8( nprocs, masterp,iproc,
     $                    nvec,ivec1,nmatstr, dmat,
     $                    icomm )
          if( lstout.gt.3 ) call TIMER( 'prhoij/dmat merged' )
        endif
c
c       Save full dmat to idmatfl for use later
        if( iproc .eq. masterp )then
          call WRITBIG( idmatfl, mat, dmat )
          if( lstout.gt.3 ) call TIMER( 'prhoij/dmat writ' )
        endif
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOIJ2
c
c
      subroutine PRHOIJ2( ivecfl,idmatfl,iematfl, norb,nstate,nk,
     $ eigval,eigpop,fpop,npop, wtk,
     $ dmat,emat, eigvec )
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis, k-parallel
c
c Written: Peter A. Schultz,  9-April-1999, for v2.34
c
c Revision history:
c  10Jul06-PAS/2.60: k-parallel version prhoij2 written
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  DMC: task parallel
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
c Scratch arrays:
      DIMENSION  fpop(nstate)
c     *** NB: dmat used as 2*norb*npop   scratch below, in addition
c         to norb^2 output.  Currently, the idle emat is contiguous
c         to dmat, so this is safe, even for npop  =norb.
      DIMENSION  dmat(norb,*),emat(norb,norb)
      DIMENSION  eigvec(2,npop  ,norb)
c
c Local declarations:
      DATA       lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
cdev:
      call FLGETIWR( IWR )
      write(IWR,*) 'DEV/PRHOIJ2-pij03: npop+pull up stripe-arith'
ctp:
c  Get local MP info:
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      if( lstout.gt.0 .and. nprocs.gt.1 .and. iproc.eq.masterp )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DMAT: k-parallel/idle-master version'
      endif
c
c    Bloch vector loop:
c
      do 1000 kblk=1,nk, nprocs
c       We will start with master(=0) doing first k-pt in set:
        kproc0 = -1
c
        kfirst = kblk
        klast = kblk + nprocs - 1
        if( klast.gt.nk )then
c         More procs than k-pts. Truncate loop
          klast = nk
c         Make masterp=0 an idle proc, so as to hide some i/o
          kproc0 = 0
        endif
c
c       >>>> Retrieve eigenvectors and send to proper worker proc
c
        kproc = kproc0
        kdmat = 0
        do 100 k=kfirst,klast
          kproc = kproc + 1
          if( kproc.eq.nprocs ) kproc = 0
c
          nvecs = 2*norb*npop
          if( iproc.eq.masterp )then
c           Retrieve (only occupied, not all nstate) eigvecs for this k-pt
            call READBIG( ivecfl, nvecs, eigvec )
            if( iproc.eq.kproc )then
c              ... and keep on master to generate dmat(k)
              call DCOPY( nvecs, eigvec,1, dmat,1 )
              kdmat = k
            else
c              ... or send to worker to generate dmat(k)
              lenmsg = nvecs
              call MPSENDR8( kproc, lenmsg, eigvec, icomm )
            endif
          elseif( iproc.eq.kproc )then
c            ... receive eigvec on worker
            lenmsg = nvecs
            call MPRECVR8( masterp, lenmsg, dmat  , icomm )
            kdmat = k
          endif
c
  100   continue
c
        if( kdmat.le.0 ) goto 500
c
c       >>>> Compute density matrix dmat, energy-weighted, emat
c
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat assignd' )
c       We are working on dmat(k):
        k = kdmat
c
c       Transpose (only occupied) eigvecs for use below:
        call TREIGVEC( norb,npop  , dmat,eigvec )
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
c       Clear out dmat, emat:
        call MKZERO( mat, dmat )
        call MKZERO( mat, emat )
c
c       Scale fpop by k-point sampling weight for use below.
        wt = wtk(k)
        do  n=1,numpop
          fpop(n) = wt*fpop(n)
        enddo
c
        do 400 j=1,norb
c
c         Do imaginary part first ...
          do 210 i=1,j-1
            dmi = zero
            emi = zero
            do 200 n=1,numpop
              occnk = fpop(n)
              ti = occnk*( eigvec(2,n,i)*eigvec(1,n,j)
     $                   - eigvec(1,n,i)*eigvec(2,n,j) )
              dmi = dmi + ti
              emi = emi + ti*eigval(n,k)
  200       continue
            dmat(i,j) = dmi
            emat(i,j) = emi
  210     continue
c
c          ... and then switch to real part:
          do 310 i=j,norb
            dmr = zero
            emr = zero
            do 300 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,i)*eigvec(1,n,j)
     $                   + eigvec(2,n,i)*eigvec(2,n,j) )
              dmr = dmr + tr
              emr = emr + tr*eigval(n,k)
  300       continue
            dmat(i,j) = dmr
            emat(i,j) = emr
  310     continue
c
c         Close loop over j/norb
  400   continue
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat compute' )
c
c       Drop here if no dmat to work on ...
  500   continue
c
c       >>>> Put dmat and emat into storage for use later
c
c       Currently, all sent to masterp and written to disk from there
        kproc = kproc0
        do 600 k=kfirst,klast
          kproc = kproc + 1
          if( kproc.eq.nprocs ) kproc = 0
c
          if( iproc .eq. masterp )then
c           Get dmats and emats from workers, and store
c           Use ticketing to keep from overwhelming master/commun.
c
c           Save dmat to idmatfl for use later
            if( iproc.ne.kproc )then
c             First, send kproc a ticket ... ok to send
              call MPSENDTICKET( kproc, icomm )
              lenmsg = mat
              call MPRECVR8( kproc, lenmsg, dmat, icomm )
            endif
            call WRITBIG( idmatfl, mat, dmat )
c
c           Save emat to iematfl for use in the force section
            if( iproc.ne.kproc )then
c             First, send kproc another ticket ... ok to send
              call MPSENDTICKET( kproc, icomm )
              lenmsg = mat
              call MPRECVR8( kproc, lenmsg, emat, icomm )
            endif
            call WRITBIG( iematfl, mat, emat )
c
            if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat writes' )
          elseif( iproc.eq.kproc )then
c           Send dmat, emat to master for storage
c           First, wait for ticket from master:
            call MPRECVTICKET( masterp, icomm )
c           Send dmat(k) to master:
            call MPSENDR8( masterp, mat, dmat, icomm )
c           Get another ticket from master:
            call MPRECVTICKET( masterp, icomm )
c           Send emat(k) to master:
            call MPSENDR8( masterp, mat, emat, icomm )
            if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat sent' )
          endif
  600   continue
c
c       End loop over blocks of k-points
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
