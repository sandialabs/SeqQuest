c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0STOREKP
c
c
      subroutine H0STOREKP( imatfile, nprocstr, master, iproc, icomm,
     $ nmatstr, nmatpr, mymatpr, mat,
     $ norb, nk, ncplx, natm, ntyp, itypa, norba,
     $ hmatstr,  hmatijk, matstr, hmatkbuf )
c -->  nmatpr    mat              ncplx*matstr
c---------------------------------------------------------------
c Purpose: Put hamiltonians to disk, converting distributed striped matrices
c          into square matrics.  Like H0STORE, but without the need to hold
c          all-k striped matrices at once on the master.
c          Requires an auxiliary buffer, hmatkbuf, of size ncplx*matstr.
c
c Written: Andrew C. Pineda,  23-November-2009, for 2.63
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(nk,*), hmatijk(norb,norb)
c NB: hmatstr only needs ~(1/nprocstr) the size of hmatstr in H0STORE.
c     It does not need to hold all the full merged (striped) matrices.
c
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Scratch array:
      DIMENSION  hmatkbuf(*)
c      hmatkbuf = Buffer space to hold unsorted striped matrix for a single k-point.
c      Its size needs to be matstr*ncplx.
c
c >>>> EXECUTABLE CODE:
c
      if( iproc .eq. master )then
        REWIND( unit=imatfile )
      endif
c
      do  ik=1, nk, 1
c
c       Collect pieces for current k-point, ik, onto master in hmatkbuf
        call STRMERGEKP( ik, nk, nprocstr, master, iproc,
     $                    icomm, ncplx, mymatpr, matstr,
     $                    hmatstr, hmatkbuf )
c
c       Sort current k-point on the master, packing the imaginary part
c       into the lower left half of the matrix if present.
c
        if( iproc .eq. master )then
          call MSORT2IJ( norb, ncplx, natm, ntyp, itypa, norba,
     $                   hmatkbuf, hmatijk )
c         Save matrix to imatfile from master only:
          call H0WRIT( imatfile, mat, 1, hmatijk )
        endif
c
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRMERGEKP
c
c
      subroutine STRMERGEKP( ik, nk, nprocstr, master, iproc,
     $                       icomm, ncplx, mymatpr, matstr,
     $                       hmat, hbuffer )
c---------------------------------------------------------------
c Purpose: Merge one k-point from the stripe-distributed matrix
c          onto master
c
c Written: Andrew C. Pineda, 23-November-2009 for 2.62.
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      INTEGER    mymatprb
      DIMENSION  hmat(*)
      DIMENSION  hbuffer(*)
c       NB: hbuffer needs to hold one full k-point, size matstr*ncplx.
c
cxkp    ***** NB: code assumes master=0 ***** where does this assumption act?
c
c  Merge the stripes for this k-point onto the master.
c  Loop over processes holding the stripes, imaginary parts, if present, get
c  packed after the real parts.
c
c >>>> EXECUTABLE CODE:
c
      jmaster = 1
      jloc = ik
c
c  Each process sends its stripe for this k-point to master process.
c  Note on each stripe the array indices/dimensions are (nk,mymatpr,ncplx).
c  So it is strided by nk.
c
      do 400 jproc = 0, nprocstr-1
         do 300 icplx=1, ncplx, 1
c
            if( jproc .eq. master )then
               if( iproc .eq. master )then
c                 This stripe computed on the master, just load it
                  do  i=jmaster, jmaster+mymatpr-1, 1
                     hbuffer(i) = hmat(jloc)
                     jloc = jloc + nk
                  enddo
                  jmaster = jmaster + mymatpr
               endif
            else
c              Non-master process did this, send stripe to the Master
               if( jproc .eq. iproc )then
                  do i=1, mymatpr,1
                     hbuffer(i) = hmat(jloc)
                     jloc = jloc + nk
                  end do
c
                  call MPSENDI( master,  1,       mymatpr, icomm ) 
                  call MPSENDR8( master, mymatpr, hbuffer, icomm )
               endif
c
c              Master process receives every stripe.
               if( iproc .eq. master )then
                  lenmsg = 1
                  call MPRECVI( jproc, lenmsg,
     $                 mymatprb,         icomm ) 
                  call MPRECVR8( jproc, mymatprb,
     $                 hbuffer(jmaster), icomm )
                  jmaster = jmaster + mymatprb
               endif
c
            endif
c
  300    continue
  400 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSORT2IJ
c
c
      subroutine MSORT2IJ( norb, ncplx, natm, ntyp, itypa, norba,
     $                     hmatstr, hmatij )
c---------------------------------------------------------------
c Purpose: Sort single k-point from the Hamiltonian matrix in striped/atom-blocked
c          into i-j form.
c
c Written: Andrew C. Pineda, 23-November-2009
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(*), hmatij(norb,norb)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declarations:
      LOGICAL    diagatm
c
c >>>> EXECUTABLE CODE:
c
c Clear out square matrix target of sort:
      nmatijk = norb*norb
      call MKZERO( nmatijk, hmatij )
c
c     Atom j loop:
c
      jjnext = 0
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        j0 = jstart
c
c       Atom i loop:
c
        ijnext = jjstart
        inext = jstart
        do 1000 iatm=jatm,natm
          ityp = itypa(iatm)
          norbati = norba(ityp)
          istart = inext
          inext = istart + norbati
c
          matij = norbati*norbatj
          ijstart = ijnext
          ijnext = ijstart + matij*ncplx
c
          i0 = istart
c
          diagatm = .false.
          if( iatm.eq.jatm ) diagatm = .true.
c
          ij0 = ijstart
          ij0 = ij0 - norbati
c
c         Atom orbital loops:
c
          iorbat1 = 1
          do 720 jorbat=1,norbatj
            j = j0 + jorbat
c
            ij0 = ij0 + norbati
            if( diagatm ) iorbat1 = jorbat
            do 710  iorbat=iorbat1,norbati
              i = i0 + iorbat
              ij = ij0 + iorbat
              ijc = ij + matij
c
                hmatij(i,j) = hmatstr(ij)
c
              if( ncplx.eq.2 .and. i.ne.j )then
                  hmatij(j,i) = hmatstr(ijc)
              endif
c
c             Close Atom orbital loops
  710       continue
  720     continue
c
c         Close Atom i loop:
 1000   continue
c
c       Close Atom j loop:
 2000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0STORE
c
c
      subroutine H0STORE( imatfile, nprocstr, master, iproc, icomm,
     $ nmatstr, nmatpr, mat,
     $ norb, nk, ncplx, natm, ntyp, itypa, norba,
     $ hmatstr, hmatijk )
c -->  nmatstr  mat
c---------------------------------------------------------------
c Purpose: Put striped/atom-blocked matrix into storage
c          Here, into k-resolved full H(n,n) onto disk
c
c Written: Peter A. Schultz,  9-October-2008, for 2.62
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Note: I pass ik into this routine simply to avoid a problem
c that had occurred (in a different place) where an optimizing
c compiler purged a loop-index from a call to H0SORT1
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(nk,*), hmatijk(norb,norb)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c >>>> EXECUTABLE CODE:
c
      call STRMERGE( nprocstr, master,iproc,icomm, nmatstr,nmatpr,
     $ hmatstr )
c -->  nmatpr-i/nmatstr-o
c
      if( iproc .eq. master )then
        REWIND( unit=imatfile )
        ik = 0
  100   continue
        ik = ik + 1
        call MSORT1IJK( norb,ik,nk,ncplx, natm,ntyp,itypa, norba,
     $   hmatstr,   hmatijk )
c -->    hmatstr-i  hmatijk-o
c -->    nmatstr    mat
c       Save matrix to imatfile:
        call H0WRIT( imatfile, mat, 1, hmatijk )
        if( ik .lt. nk ) goto 100
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRMERGE
c
c
      subroutine STRMERGE( nprocstr, master,iproc,icomm, nmatstr,nmatpr,
     $ hmat )
c---------------------------------------------------------------
c Purpose: merge distributed matrix onto master
c Written: Peter A. Schultz, 4-October-2008, for 2.62
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  hmat(*)
c
      if( nprocstr .gt. 1 .and. iproc .lt. nprocstr )then
c       Merge the stripes to full set of stripes on master:
        nmatpr1 = 1
        call MPMERGER8( nprocstr, master, iproc,
     $   nmatstr,nmatpr1,nmatpr, hmat, icomm )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSORTIJK
c
c
      subroutine MSORTIJK( norb,nk,ncplx, natm,ntyp,itypa, norba,
     $ hmatstr,hmatijk )
c---------------------------------------------------------------
c Purpose: sort Hamiltonian matrix from striped/atom-blocked
c          into k-resolved i-j form.
c
c Written: Peter A. Schultz, 25-Septemeber-2008
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(nk,*), hmatijk(norb,norb,nk)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declarations:
      LOGICAL    diagatm
c
c >>>> EXECUTABLE CODE:
c
c Clear out square matrix target of sort:
      nmatijk = norb*norb*nk
      call MKZERO( nmatijk, hmatijk )
c
c     Atom j loop:
c
      jjnext = 0
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        j0 = jstart
c
c       Atom i loop:
c
        ijnext = jjstart
        inext = jstart
        do 1000 iatm=jatm,natm
          ityp = itypa(iatm)
          norbati = norba(ityp)
          istart = inext
          inext = istart + norbati
c
          matij = norbati*norbatj
          ijstart = ijnext
          ijnext = ijstart + matij*ncplx
c
          i0 = istart
c
          diagatm = .false.
          if( iatm.eq.jatm ) diagatm = .true.
c
          ij0 = ijstart
          ij0 = ij0 - norbati
c
c         Atom orbital loops:
c
          iorbat1 = 1
          do 720 jorbat=1,norbatj
            j = j0 + jorbat
c
            ij0 = ij0 + norbati
            if( diagatm ) iorbat1 = jorbat
            do 710  iorbat=iorbat1,norbati
              i = i0 + iorbat
              ij = ij0 + iorbat
              ijc = ij + matij
c
              do  k=1,nk
                hmatijk(i,j,k) = hmatstr(k,ij)
              enddo
c
              if( ncplx.eq.2 .and. i.ne.j )then
                do  k=1,nk
                  hmatijk(j,i,k) = hmatstr(k,ijc)
                enddo
              endif
c
c             Close Atom orbital loops
  710       continue
  720     continue
c
c         Close Atom i loop:
 1000   continue
c
c       Close Atom j loop:
 2000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSORT1IJK
c
c
      subroutine MSORT1IJK( norb,ik,nk,ncplx, natm,ntyp,itypa, norba,
     $ hmatstr,hmatijk )
c---------------------------------------------------------------
c Purpose: sort Hamiltonian matrix from striped/atom-blocked
c          into k-resolved i-j form.
c
c Written: Peter A. Schultz, 25-Septemeber-2008
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(nk,*), hmatijk(norb,norb)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declarations:
      LOGICAL    diagatm
c
c >>>> EXECUTABLE CODE:
c
      k = ik
      if( k.lt.1 .or. k.gt.nk )then
c       The compiler killed the passed k
        call FLGETIERR( IERR )
        write(IERR,*) 'MSORT1 receives corrupted k-index'
        write(IERR,*) 'Likely opting compiler purged argument in call'
        call STOPXERR( 'MSORT1IJK - trouble with compiler?' )
      endif
c
c Clear out square matrix target of sort:
      nmatijk = norb*norb
      call MKZERO( nmatijk, hmatijk )
c
c     Atom j loop:
c
      jjnext = 0
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        j0 = jstart
c
c       Atom i loop:
c
        ijnext = jjstart
        inext = jstart
        do 1000 iatm=jatm,natm
          ityp = itypa(iatm)
          norbati = norba(ityp)
          istart = inext
          inext = istart + norbati
c
          matij = norbati*norbatj
          ijstart = ijnext
          ijnext = ijstart + matij*ncplx
c
          i0 = istart
c
          diagatm = .false.
          if( iatm.eq.jatm ) diagatm = .true.
c
          ij0 = ijstart
          ij0 = ij0 - norbati
c
c         Atom orbital loops:
c
          iorbat1 = 1
          do 720 jorbat=1,norbatj
            j = j0 + jorbat
c
            ij0 = ij0 + norbati
            if( diagatm ) iorbat1 = jorbat
            do 710  iorbat=iorbat1,norbati
              i = i0 + iorbat
              ij = ij0 + iorbat
              ijc = ij + matij
c
                hmatijk(i,j) = hmatstr(k,ij)
c
              if( ncplx.eq.2 .and. i.ne.j )then
                  hmatijk(j,i) = hmatstr(k,ijc)
              endif
c
c             Close Atom orbital loops
  710       continue
  720     continue
c
c         Close Atom i loop:
 1000   continue
c
c       Close Atom j loop:
 2000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSORT1STR
c
c
      subroutine MSORT1STR( norb,ik,nk,ncplx, natm,ntyp,itypa, norba,
     $ hmatstr,hmatijk )
c---------------------------------------------------------------
c Purpose: sort Hamiltonian matrix into striped/atom-blocked
c          from k-resolved i-j form.
c
c Written: Peter A. Schultz, 14-October-2008
c
c Revision history:
c  14Oct08-PAS/2.62: adapted from msort1ijk
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  hmatstr(nk,*), hmatijk(norb,norb)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declarations:
      LOGICAL    diagatm
c
c >>>> EXECUTABLE CODE:
c
      k = ik
      if( k.lt.1 .or. k.gt.nk )then
c       The compiler killed the passed k
        call FLGETIERR( IERR )
        write(IERR,*) 'MSORT1 receives corrupted k-index'
        write(IERR,*) 'Likely opting compiler purged argument in call'
        call STOPXERR( 'MSORT1STR - trouble with compiler?' )
      endif
c
c     Atom j loop:
c
      jjnext = 0
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        j0 = jstart
c
c       Atom i loop:
c
        ijnext = jjstart
        inext = jstart
        do 1000 iatm=jatm,natm
          ityp = itypa(iatm)
          norbati = norba(ityp)
          istart = inext
          inext = istart + norbati
c
          matij = norbati*norbatj
          ijstart = ijnext
          ijnext = ijstart + matij*ncplx
c
          i0 = istart
c
          diagatm = .false.
          if( iatm.eq.jatm ) diagatm = .true.
c
          ij0 = ijstart
          ij0 = ij0 - norbati
c
c         Atom orbital loops:
c
          iorbat1 = 1
          do 720 jorbat=1,norbatj
            j = j0 + jorbat
c
            ij0 = ij0 + norbati
            if( diagatm ) iorbat1 = jorbat
            do 710  iorbat=iorbat1,norbati
              i = i0 + iorbat
              ij = ij0 + iorbat
              ijc = ij + matij
c
                hmatstr(k,ij) = hmatijk(i,j)
c
              if( ncplx.eq.2 .and. i.ne.j )then
                  hmatstr(k,ijc) = hmatijk(j,i)
              endif
c
c             Close Atom orbital loops
  710       continue
  720     continue
c
c         Close Atom i loop:
 1000   continue
c
c       Close Atom j loop:
 2000 continue
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
