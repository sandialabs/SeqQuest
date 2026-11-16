c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FFTDRHO
c
c
      subroutine FFTDRHO( idrhofl,  n1r,n2r,n3r,nptr,
     $ ws1,ws2,ws3,ftarray,
     $ rho,gvec, rhog, drhog )
c---------------------------------------------------------------
c Purpose: compute density derivatives from grid density via fft
c
c Written: Peter A. Schultz, 9-November-2000, for v2.46 (gga)
c
c Revision history:
c  26Jun07-PAS/2.60: task-parallel
c  20Jul01-PAS/2.49: spin-polarized dft, move rewind up.
c---------------------------------------------------------------
c
c First and second derivatives of the density are written to
c disk file unit "idrhofl" in order: X,Y,Z, XX,YY,ZZ, XY,YZ,ZX
c
c Can either keep, or destroy, rho() depending on "keeprh".
c rho() is preserved in this code.
c
c Transform the density p(r) into reciprocal space via FFT:
c    p(r) = SUM[G]: C[G] * exp(iG.r)
c C[G] are complex values obtained via forward fft of p(r)
c The real space derivatives of p(r) are taken in recip space:
c    d(p(r))/du = SUM[G]: (iGu*C[G]) * exp(iG.r)
c    d2(p)/dudv = SUM[G]: (-GuGvC[G]) * exp(iG.r)
c where u,v=x,y,z, and G=(Gx,Gy,Gz), a real three-vector.
c C'[G] = iGu*C[G] and C''[G] = -Gu*Gv*C[G] are the expansion
c coeffs in reciprocal space and are used to get derivatives in
c real space via back fft.  I do the forward fft once and save,
c and do repeated back fft's to collect all the derivatives.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c fft arrays:
      DIMENSION  ws1(*),ws2(*),ws3(*), ftarray(2,*)
c Input arrays: (density, recip space vectors)
      DIMENSION  rho(nptr), gvec(3,nptr)
c Scratch arrays: (fft of density, and fft derivs of density)
      DIMENSION  rhog(2,nptr)
      DIMENSION  drhog(2,nptr)
c
c Local declarations:
      LOGICAL    keeprh
c
c >>>> EXECUTABLE CODE:
c
c Do we destroy rho(), or keep rho() and compact drhog()?
      keeprh = .true.
c
c Load density into real part of complex array for the fft:
      call MKZERO( 2*nptr, rhog )
      call DCOPY( nptr, rho,1, rhog,2 )
c
c Set up parallelism
c   The derivatives and back fft's will be distributed over
c   processors, and the result shipped back to the master
c   for storage.
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
c  Assume master is first (0) and procs sequential ranked
      iprocu1 = 0
c  With enough procs, make master fft-idle and just do i/o
      if( nprocs .gt. 3 ) iprocu1 = 1
c
c  Fourier transform density:
c
      if( iproc.ne.master .or. iprocu1.eq.0 )
     $call CFFT3D( rhog, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c
c >>>> First derivatives: DX, DY, DZ
c
      iprocu = nprocs
      do 101 iu=1,3
        iprocu = iprocu + 1
        if( iprocu .ge. nprocs ) iprocu = iprocu1
c
        if( iproc .eq. iprocu )then
          do  ng=1,nptr
            drhog(1,ng) = -rhog(2,ng)*gvec(iu,ng)
            drhog(2,ng) =  rhog(1,ng)*gvec(iu,ng)
          enddo
          call CFFT3DB( drhog, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
        endif
c
        if( iproc .eq. iprocu  .or.  iproc .eq. master )then
          lenmsg = nptr
c
          if( keeprh )then
c           Compact transformed array to keep only the real part
            if( iproc .eq. iprocu )then
              iq = 2
              id = -1
              ir = 0
              do  nr=1,nptr
                ir = ir + iq - 1
                iq = iq + id
                id = -id
                drhog(iq,ir) = drhog(1,nr)
              enddo
              if( iproc .ne. master )then
                call MPSENDR8( master, lenmsg, drhog, icomm )
              endif
            elseif( iproc .eq. master )then
              call MPRECVR8( iprocu, lenmsg, drhog, icomm )
            endif
            if( iproc .eq. master )then
              call WRITBIG( idrhofl, nptr, drhog )
            endif
          else
c           Overwrite the density array
            if( iproc .eq. iprocu )then
              do  nr=1,nptr
                rho(nr) = drhog(1,nr)
              enddo
              if( iproc .ne. master )then
                call MPSENDR8( master, lenmsg, rho, icomm )
              endif
            elseif( iproc .eq. master )then
              call MPRECVR8( iprocu, lenmsg, rho, icomm )
            endif
            if( iproc .eq. master )then
              call WRITBIG( idrhofl, nptr, rho )
            endif
          endif
c
        endif
c
  101 continue
c
c >>>> Second derivatives: XX, YY, ZZ, XY, YZ, ZX
c
      do 202 idv=0,2
        do 201 iu=1,3-idv
          iv = iu + idv
          iprocu = iprocu + 1
          if( iprocu .ge. nprocs ) iprocu = iprocu1
c
          if( iproc .eq. iprocu )then
            do  ng=1,nptr
              drhog(1,ng) = -rhog(1,ng)*gvec(iu,ng)*gvec(iv,ng)
              drhog(2,ng) = -rhog(2,ng)*gvec(iu,ng)*gvec(iv,ng)
            enddo
            call CFFT3DB( drhog, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
          endif
c
          if( iproc .eq. iprocu  .or.  iproc .eq. master )then
            lenmsg = nptr
c
            if( keeprh )then
c             Compact transformed array to keep only the real part
              if( iproc .eq. iprocu )then
                iq = 2
                id = -1
                ir = 0
                do  nr=1,nptr
                  ir = ir + iq - 1
                  iq = iq + id
                  id = -id
                  drhog(iq,ir) = drhog(1,nr)
                enddo
                if( iproc .ne. master )then
                  call MPSENDR8( master, lenmsg, drhog, icomm )
                endif
              elseif( iproc .eq. master )then
                call MPRECVR8( iprocu, lenmsg, drhog, icomm )
              endif
              if( iproc .eq. master )then
                call WRITBIG( idrhofl, nptr, drhog )
              endif
            else
c             Overwrite density array
              if( iproc .eq. iprocu )then
                do  nr=1,nptr
                  rho(nr) = drhog(1,nr)
                enddo
                if( iproc .ne. master )then
                  call MPSENDR8( master, lenmsg, rho, icomm )
                endif
              elseif( iproc .eq. master )then
                call MPRECVR8( iprocu, lenmsg, rho, icomm )
              endif
              if( iproc .eq. master )then
                call WRITBIG( idrhofl, nptr, rho )
              endif
            endif
c
          endif
c
  201   continue
  202 continue
c
c    That's all Folks!
c
      RETURN
      END
