c
c 
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DM0FETCH
c
c     
      subroutine DM0FETCH( idmatfl,idmatsfl, nspin,
     $ norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ dmatstr,dmatij )
c---------------------------------------------------------------
c Purpose: retrieve full density (or energy-weighted) matrices,
c          and send to processors in striped form.
c     
c Written: Peter A. Schultz, 13-October-2008, for v2.62 (stripe-par)
c     
c Revision history:
c  22Jan12-PAS/2.63: harmonized with ACP/dm0fetchkp, installed below
c---------------------------------------------------------------
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c     
      PARAMETER  ( zero = 0.d0, one = 1.d0 )
c
c Input arrays:
      DIMENSION  itypa(natm),norba(ntyp)
c Output array:  dmastr[nmatstr]  spin&k=1: max(nmatstr,mat)
      DIMENSION  dmatstr(*)
c
c Scratch array: dmatij[mat(k=1)|nspin*mat(k>1)]
      DIMENSION  dmatij(*)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      mat = norb*norb
      if( iproc .eq. master )then
        REWIND( unit=idmatfl )
        if( nspin.eq.2 ) REWIND( unit=idmatsfl )
      endif
c
      ik = 0
  100 continue
      ik = ik + 1
c
c  The master retrieves amd sorts the matrices (up&dn spin summed)
      mat = norb*norb
      if( iproc .eq. master )then
        call READBIG( idmatfl, mat, dmatij )
        if( nspin .eq. 2 )then
c         Read other spin and sum into this:
          if( nk .gt. 1 )then
c           We need to use an extra buffer to avoid wiping dmatstr
            call READBIG( idmatsfl, mat, dmatij(1+mat) )
            call DAXPY( mat, one, dmatij(1+mat),1, dmatij,1 )
          else
c           With nk=1, dmatstr is idle, overlap use as buffer here
            call READBIG( idmatsfl, mat, dmatstr )
            call DAXPY( mat, one, dmatstr,1, dmatij,1 )
          endif
        endif
c        ... and inserts this dmat(ik) into full nk-striped matrix:
        call MSORT1STR( norb,ik,nk,ncplx, natm,ntyp,itypa,norba,
     $   dmatstr,dmatij )
      endif
c
      if( ik .lt. nk ) goto 100
c
c  Master sends the correct stripe to the all processors:
      call MSENDSTR( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ nprocs,iproc,master,icomm,
     $ dmatstr )
c
c    That's all Folks!
c
      RETURN
      END
c
c    
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSENDSTR
c
c
      subroutine MSENDSTR( norb,nk,ncplx, natm,ntyp,itypa,norba,
     $ nprocs,iproc,master,icomm,
     $ dmatstr )
c---------------------------------------------------------------
c Purpose: Send matrix stripes out to their respective owners
c     
c Written: Peter A. Schultz, 14-October-2008
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
      DIMENSION  dmatstr(*)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c >>>> EXECUTABLE CODE:
c     
c Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
      nmatstr = ncplx*nk*matstr
c
c Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
      nmatpr = mymatpr*ncplx*nk
c
      if( nprocstr.le.1 ) goto 900
c
c Master will start sending AFTER its own allotment:
      jjnext = nmatpr
c
      do 1000 jproc=1,nprocstr-1
c
        if( iproc.eq.master )then
c         Find length to be sent, and from what starting point:
          call MATLOCSTR( jproc,nprocs,nprocstr,
     $      jatm1,  jatm2, jpj0, maxmatpr, jmatpr ,
     $     matstr,
     $     norb,nk,ncplx, natm,ntyp,itypa, norba )
          nmatjpr = jmatpr*ncplx*nk
          ij1 = jjnext + 1
          jjnext = jjnext + nmatjpr
          lenmsg = nmatjpr
          call MPSENDR8( jproc, lenmsg, dmatstr(ij1), icomm )
c
        elseif( iproc.eq.jproc )then
c         Receive this stripe from the master:
          ij1 = 1
          lenmsg = nmatpr
          call MPRECVR8( master, lenmsg, dmatstr(ij1), icomm )
c
        endif
c
 1000 continue
c
  900 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c 
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DM0FETCHKP
c
c
      subroutine DM0FETCHKP( idmatfl, idmatsfl, nspin,
     $ norb, nk, ncplx, natm, ntyp, itypa, norba,
     $ dmatstr,              dmatij,       dmatstrbuf )
c      nk*mymatpr*ncplx-o    mat-s         nc*matstr-s
c      = nmatpr                            spin: & mat-s
c---------------------------------------------------------------
c Purpose: retrieve full density (or energy-weighted) matrices,
c          and send to processors in striped form.
c          (Inverse of H0STOREKP.)
c     
c Written: Andrew C. Pineda, 7-January-2010, for 2.62 (k-blocked,
c          stripe-par)
c
c Total space required:
c    no-spin:  nmatpr + mat + ncplx*matstr
c    spin:     nmatpr + mat + max( ncplx*matstr, mat )
c  note: could be packed tighter still (nmatpr less) for nk=1,
c        by overlapping needed scratch space with output space
c     
c Revision history:
c  7-Jan-2010/ACP: Derived from DM0FETCH.
c---------------------------------------------------------------
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c     
      PARAMETER  ( zero = 0.d0, one = 1.d0 )
c
c Input arrays:
      DIMENSION  itypa(natm), norba(ntyp)
c
c Output array:  dmatstr[nmatpr], nmatpr = ncplx*nk*mymatpr 
      DIMENSION  dmatstr(*)
c
c Scratch array: dmatij[mat]
      DIMENSION  dmatij(*)
c Scratch array: dmatstrbuf: ncplx*matstr, spin: max( buf, mat )
      DIMENSION  dmatstrbuf(*)
c    NB:  matstr <= mat so can use a buffer of ncplx*mat here to use 
c                dmatstrbuf for other purposes
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      mat = norb*norb
      if( iproc .eq. master )then
        REWIND( unit=idmatfl )
        if( nspin.eq.2 ) REWIND( unit=idmatsfl )
      endif
c
      ik = 0
  100 continue
      ik = ik + 1
c
c The master retrieves the matrix (up&dn spin summed) one k-point at a time.
c
      mat = norb*norb
      if( iproc .eq. master )then
        call READBIG( idmatfl, mat, dmatij )
        if( nspin .eq. 2 )then
c         Read other spin and sum into this
          call READBIG( idmatsfl, mat, dmatstrbuf )
          call DAXPY( mat, one, dmatstrbuf,1, dmatij,1 )
        endif
c
c Inserts this dmat(ik) into striped matrix for 1 k-point on the master.
c
        ikphony = 1
        nkphony = 1
        call MSORT1STR( norb, ikphony, nkphony, ncplx, natm,
     $   ntyp, itypa, norba,
     $   dmatstrbuf,       dmatij )
c -->    ncplx*matstr-o    mat-i
c
      endif
c
c Scatter the stripe across processors to the full distributed nk-striped matrix.
c
      call MSCATTERSTRKP( norb, ik, nk, ncplx, natm,
     $ ntyp, itypa, norba,
     $ nprocs, iproc, master, icomm,
     $ dmatstr,           dmatstrbuf )
c -->  ncplx*mymatpr*nk   ncplx*matstr-i
c -->  = nmatpr
c
      if( ik .lt. nk ) goto 100
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MSCATTERSTRKP
c
c
      subroutine MSCATTERSTRKP( norb, ik, nk, ncplx, natm,
     $     ntyp, itypa, norba,
     $     nprocs, iproc, master, icomm,
     $     dmatstr, dmatstrbuf )
c---------------------------------------------------------------
c Purpose: Send pieces of a matrix stripe out to their respective owners
c     
c Written: Andrew C. Pineda, 7-Jan-2010, for 2.62
c
c Total space required: (nk*mymatpr*ncplx + matstr*ncplx)
c
c Revision history:
c   7-Jan-2010/ACP:  Derived from MSENDSTR
c  21-Jan-2010/ACP:  Both implementations now work.
c  29-Jan-2010/ACP:  Cleaned out development code and commented out original code.
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output arrays:
c     dmatstr[nk*mymatpr*ncplx], mymatpr is local to processor
      DIMENSION  dmatstr(*)
c     Array buffer for a single stripe or partial stripe (1 k-point)
c     dmatstrbuf[matstr*ncplx] on master, dmatstrbuf[mymatpr*ncplx] elsewhere
      DIMENSION  dmatstrbuf(*)
c Input arrays:
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
      LOGICAL REFERENCE_IMPLEMENTATION
      DATA REFERENCE_IMPLEMENTATION /.false./
c
c >>>> EXECUTABLE CODE:
c     
c Determine size of a striped matrix for 1 k-point
c
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c      nmatstr = ncplx*nk*matstr
c
c Figure out how much belongs on local processor:
      call MATLOCSTR( iproc, nprocs, nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c      nmatpr = mymatpr*ncplx*nk
      matpr = mymatpr*ncplx
c
c      if( nprocstr.le.1 ) goto 900
c
c Master will start sending AFTER its own allotment:
c      jjnext = matpr
c
c      ij_send_offset=1
      ij_send_offset=1
      do 1000 jproc=0, nprocstr-1
         if(REFERENCE_IMPLEMENTATION) then

            if(iproc.eq.jproc) then
c Receiving process sends how much data it needs unless it is also the master process.
               lenmsg=1
               if(iproc.ne.master) then
                  call MPSENDI(master, lenmsg, matpr, icomm)
               else
                  ij_recv_offset=ij_send_offset
               endif
            endif

            if(iproc.eq.master) then
c Master process receives the length of the data it needs to send and sends it.
               lenmsg=1
               if(iproc.ne.jproc) then
                  call MPRECVI(jproc, lenmsg, matjpr, icomm)
                  lenmsg=matjpr
                  call MPSENDR8(jproc, lenmsg,
     $                 dmatstrbuf(ij_send_offset), icomm)
               endif
               ij_send_offset=ij_send_offset+matjpr
            endif

            if(iproc.eq.jproc) then
               lenmsg=matpr
               if(iproc.ne.master) then
c If receiving process is not the master get the data from the master
                  ij_recv_offset=1
                  call MPRECVR8( master, lenmsg,
     $                 dmatstrbuf(ij_recv_offset),
     $                 icomm )
               endif

c Receiving process puts the data in the correct place in local memory.
               ioffset=ik
               do ij3=0, matpr-1
                  dmatstr(ioffset)=dmatstrbuf(ij_recv_offset+ij3)
                  ioffset=ioffset+nk
               end do

            endif
         else
c     
            if( iproc.eq.master )then
c    Find length of stripe to be sent from master to jproc
               call MATLOCSTR( jproc, nprocs, nprocstr,
     $              jatm1,  jatm2, jpj0, maxmatpr, jmatpr ,
     $              matstr,
     $              norb, nk, ncplx, natm, ntyp, itypa, norba )
               matjpr = jmatpr*ncplx
               if(jproc.ne.master) then
c     Destination (jproc) is not master, do a send operation.
                  lenmsg = matjpr
                  call MPSENDR8( jproc, lenmsg,
     $                 dmatstrbuf(ij_send_offset),
     $                 icomm )
               else
c   Destination is on master, so the master only needs to know where the 
c   data starts to unpack it below. (jproc.eq.master)
                  ij_recv_offset=ij_send_offset
                  lenmsg = matjpr
               endif
c   Increment send offset to next block of data for next pass.
               ij_send_offset = ij_send_offset + matjpr
            endif
c
c   Receive this stripe if this is the destination
c
            if( iproc.eq.jproc )then
c   Do receive operation if I'm not the master. Can use all of dmatstrbuf.
               if(iproc.ne.master) then
                  ij_recv_offset=1
                  lenmsg = matpr
                  call MPRECVR8( master, lenmsg,
     $                 dmatstrbuf(ij_recv_offset),
     $                 icomm )
               endif
c     Put the stripe in its proper location on this processor
c     in dmatstr. This is just starting at offset = ik and strided by nk.

               ioffset=ik
               do ij3=0, matpr-1
                  dmatstr(ioffset)=dmatstrbuf(ij_recv_offset+ij3)
                  ioffset=ioffset+nk
               end do
c     
            endif
c     
         endif
 1000 continue
c
  900 continue
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
