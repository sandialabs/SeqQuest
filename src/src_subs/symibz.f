c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMIBZ
c
c
      subroutine SYMIBZ( ndim, nk0,vk0,wtk0, nkd,nk,vk,wtk,wtkscal,
     $ rvec, nsym, rmatsym )
c---------------------------------------------------------------
c Purpose: reduce k-points to IBZ using symmetry
c
c Written: Peter A. Schultz, 19-October-1999, for v2.40
c
c Revision history:
c   1Aug07-PAS/2.60: another little dimension thing
c   6Aug03-PAS/2.55a: extract k-write; install kcell generate
c  25Jul01-PAS/2.49: wtk norm convention
c  21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c Variables:
c   rvec(3,3) - primitive lattice vectors
c   nk0       - number of input k-points
c   vk0(3,*)  - input k-points
c   wtk0(3,*) - weights for input k-points (not used right now)
c   nk       - number of output k-points
c   vk(3,*)  - output k-points
c   wtk(3,*) - weights for output k-points
c   nsym     - number of (non-identity) symmetries
c   rmatsym  - cartesian representations of symmetries
c
c Note: must add time reversal/inversion to list of symmetries
c       to check against.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  rvec(3,3)
      DIMENSION  vk0(3,nk0), wtk0(nk0)
      DIMENSION  rmatsym(3,4,*)
c output arrays:
      DIMENSION  vk(3,*), wtk(*)
c local declarations:
      DIMENSION  ri(3),sri(3)
      DIMENSION  gvec(3,3), dk(3),cg(3)
      DATA  zero,one / 0.d0,1.d0 /
      DATA  extra / 1.d-6 /, bigdif / 1.d-10 /
c
c >>>> EXECUTABLE CODE:
c
c Produce primitive reciprocal vectors:
      call KCELL( ndim, rvec, gvec )
c
c Take the first k-point for sure:
      nk = 1
      wtk(1) = one
      vk(1,1) = vk0(1,1)
      vk(2,1) = vk0(2,1)
      vk(3,1) = vk0(3,1)
c
c Now cycle through rest of k-points ...
c
      do 1000 ik0=2,nk0
        ri(1) = vk0(1,ik0)
        ri(2) = vk0(2,ik0)
        ri(3) = vk0(3,ik0)
c
c  ... transform under symmetry operations ...
c (... adding a time-reversal/inversion ...)
c
        do 300 inv=1,2
c
          do 200 isym=0,nsym
c
c Compute transformed vector:
            if( isym.eq.0 )then
              sri(1) = ri(1)
              sri(2) = ri(2)
              sri(3) = ri(3)
            else
              sri(1) = rmatsym(1,1,isym)*ri(1) + rmatsym(2,1,isym)*ri(2)
     $               + rmatsym(3,1,isym)*ri(3)
              sri(2) = rmatsym(1,2,isym)*ri(1) + rmatsym(2,2,isym)*ri(2)
     $               + rmatsym(3,2,isym)*ri(3)
              sri(3) = rmatsym(1,3,isym)*ri(1) + rmatsym(2,3,isym)*ri(2)
     $               + rmatsym(3,3,isym)*ri(3)
            endif
c
c  ... and try to match to vectors already in set, modulo the gvecs
            do 100 ik=1,nk
              dk(1) = vk(1,ik) - sri(1)
              dk(2) = vk(2,ik) - sri(2)
              dk(3) = vk(3,ik) - sri(3)
              call PROJLV( dk, gvec(1,1),gvec(1,2),gvec(1,3), cg )
              difsq = zero
              do  jd=1,3
                if( cg(jd).lt.zero ) cg(jd) = -cg(jd)
                cg(jd) = cg(jd) + extra
                igfac = cg(jd)
                difsq = difsq + (cg(jd)-DBLE(igfac))**2
              enddo
              if( difsq.lt.bigdif )then
c               We have a match
                wtk(ik) = wtk(ik) + one
                goto 1000
              endif
  100       continue
c
  200     continue
c
c         Install time reversal:
          ri(1) = -ri(1)
          ri(2) = -ri(2)
          ri(3) = -ri(3)
c
  300   continue
c
c No match, hence new symmetry-unique k-point
        nk = nk + 1
        if( nk.gt.nkd ) goto 1301
        wtk(nk) = one
        vk(1,nk) = vk0(1,ik0)
        vk(2,nk) = vk0(2,ik0)
        vk(3,nk) = vk0(3,ik0)
 1000 continue
c
c Scale k-points weights to internal conventions:
      wtkscal = zero
      do  k=1,nk
        wtkscal = wtkscal + wtk(k)
      enddo
      scalek = one/wtkscal
      do  k=1,nk
        wtk(k) = scalek*wtk(k)
      enddo
c
      RETURN
c
 1301 continue
      call STOPXERR( 'k#ibz-d - max number kpt/ibz exceeded' )
c
c    That's all Folks!
c
      END
