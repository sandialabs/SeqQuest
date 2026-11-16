c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOFIX
c
c
      subroutine VSLOFIX( near2c, s1atom, orbint,orbofk,
     $ norb,nk, natm,ntyp,nshld, noad, itypa,numshl,lshel,
     $ vslnear, ncplx,
     $ vmat )
c---------------------------------------------------------------
c Purpose: fix grid matrix elements with nearby construction
c
c Written: Peter A. Schultz, 29-April-1997, for v2.17
c
c Revision history:
c   4Dec08-PAS/2.63: k-parallel scf matrices
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   9Sep98-PAS/2.25: bug fix (argh), clean and match
c    Mar98-PAS/2.22: allow option to turn off 2ctr nearby
c   3Jun97-PAS/2.20: cut old grid cutoff stuff for boxes
c---------------------------------------------------------------
c
c Notes: Checking whether atoms overlap can be skipped - it is
c        a cheap n**2 problem here.
c
c     ACP Notes: vmat is not initialized to zero here, rather if needed
c     this done in the calling routine. Work is done at the image level.
c     (Processes within the same k-point communicator should compute the
c     exact same quantities. Also processes do not communicate as part
C     of the computation here, so if no rounding errors have accumulated
c     in the inputs, the results should agree exactly with the serial
C     computation.) 
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  orbint(2,norb,nk),orbofk(2,natm,norb,nk)
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  vslnear(natm)
c Input/Output arrays:
      DIMENSION  vmat(norb,norb,*)
c
c Local declarations:
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( near2c .eq. 0 ) RETURN
c
      call MPNODES( nprocl )
      call MPNODE( iprocl )
c
      call KPMINE( nk,nkloc,nk0 )
c     Check if there is nothing to do for this proc:
      if( nkloc .eq. 0 ) RETURN
c
      inext = 0
      do 200 iatm=1,natm
        ityp = itypa(iatm)
c
        iorbat00 = 0
        do 180 ishl=1,numshl(ityp)
          norbshi = 2*lshel(ishl,ityp) + 1
          i0 = inext
          inext = inext + norbshi
          iorbat0 = iorbat00
          iorbat00 = iorbat00 + norbshi
c
          jnext = 0
          do 100 jatm=1,iatm
            jtyp = itypa(jatm)
            nshlj = numshl(jtyp)
            if( jatm.eq.iatm ) nshlj = ishl
c
            vnri = vslnear(iatm)
            vnrj = vslnear(jatm)
c
            jorbat00 = 0
            do 80 jshl=1,nshlj
              norbshj = 2*lshel(jshl,jtyp) + 1
              j0 = jnext
              jnext = jnext + norbshj
              jorbat0 = jorbat00
              jorbat00 = jorbat00 + norbshj
c
              joff = 0
              if( jatm.eq.iatm .and. jshl.eq.ishl )then
                norbshj = 0
                joff = 1
              endif
c
              i = i0
              iorbat = iorbat0
              do 170 iorbsh=1,norbshi
                iorbat = iorbat + 1
                i = i + 1
c
                norbshj = norbshj + joff
c
                j = j0
                jorbat = jorbat0
                do 70 jorbsh=1,norbshj
                  jorbat = jorbat + 1
                  j = j + 1
c
                  do 50 ikloc=1,nkloc
                    ik = nk0 + ikloc
                    vintr = zero
                    vinti = zero
c
                    if( iatm.eq.jatm )then
                      vintr = vnri*s1atom(iorbat,jorbat,ik,iatm)
                      if( iorbat.ne.jorbat )then
                        vinti = vnri*s1atom(jorbat,iorbat,ik,iatm)
                      endif
                    elseif( near2c .eq. 2 )then
                      si2r = orbint(1,i,ik)*orbofk(1,iatm,j,ik)
     $                     + orbint(2,i,ik)*orbofk(2,iatm,j,ik)
                      si2i = orbint(2,i,ik)*orbofk(1,iatm,j,ik)
     $                     - orbint(1,i,ik)*orbofk(2,iatm,j,ik)
                      sj2r = orbofk(1,jatm,i,ik)*orbint(1,j,ik)
     $                     + orbofk(2,jatm,i,ik)*orbint(2,j,ik)
                      sj2i = orbofk(2,jatm,i,ik)*orbint(1,j,ik)
     $                     - orbofk(1,jatm,i,ik)*orbint(2,j,ik)
                      vintr = vnri*si2r + vnrj*sj2r
                      vinti = vnri*si2i + vnrj*sj2i
                    endif
c
                    vmat(i,j,ikloc) = vmat(i,j,ikloc) + vintr
                    if( j.ne.i .and. ncplx.eq.2 )then
                      vmat(j,i,ikloc) = vmat(j,i,ikloc) + vinti
                    endif
c
   50             continue
c
   70           continue
c
  170         continue
c
   80       continue
  100     continue
c
  180   continue
  200 continue
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
