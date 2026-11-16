c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ROTGRP
c
c
      subroutine ROTGRP( dmatwk, nuc,jatm,rij,
     $ norb,nk, natm,ntyp,nshld, noad,
     $ itypa,numshl,lshel, ioffa, vk,rl,
     $ nsyms,dsym,naofsym,
     $ dmat,cphi )
c---------------------------------------------------------------
c Purpose: take the complex input matrix "dmat", known only in
c   irreducible part of the surface BZ, and use input info
c   concerning the rotations and reflections under which slab is
c   symmetric, determine real configuration space matrix dmatwk
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output array:
      DIMENSION  dmatwk(noad,noad)
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  ioffa(natm)
      DIMENSION  rij(2),vk(3,nk),rl(2), dsym(5,5,3,*),naofsym(natm,*)
      DIMENSION  dmat(2,norb,norb,nk)
c scratch array:
      DIMENSION  cphi(2,nk)
complex dmat,cphi
c local declarations:
      DIMENSION  rrot(2)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      symops = one/ DBLE(nsyms)
      do 1000 isym=1,nsyms
c
c      compute phase (recall that dsym are inverse rotation matrices)
c
c  10/10/85 - the phase associated with each symmetry operator is:
c    exp(i*vk*(a lattice vector)), the lattice vector is the change
c    in rij due to action of inverse symmetry operator in question
c
        do  ir=1,2
          rrot(ir) = -(rij(ir)+rl(ir))
          do  jr=1,2
            rrot(ir) = rrot(ir) + dsym(ir,jr,2,isym)*rij(jr)
          enddo
        enddo
        do  k=1,nk
          arg = (vk(1,k)*rrot(1)+vk(2,k)*rrot(2))
c          cphi(k) = DCMPLX( COS(arg)*symops, SIN(arg)*symops )
          cphi(1,k) = symops*COS(arg)
          cphi(2,k) = symops*SIN(arg)
        enddo
c
c      specify atoms in basis
c
        nucsym = naofsym(nuc,isym)
        jatsym = naofsym(jatm,isym)
        ktyp = itypa(nuc)
        jtyp = itypa(jatm)
        nshli = numshl(ktyp)
        nshlj = numshl(jtyp)
c
        norbwki = 0
        norbi = ioffa(nucsym)
c
        do 85 ishl=1,nshli
          lorbi = lshel(ishl,ktyp)
          norbshi = 2*lorbi + 1
c
          norbwkj = 0
          norbj = ioffa(jatsym)
          do 80 jshl=1,nshlj
            lorbj = lshel(jshl,jtyp)
            norbshj = 2*lorbj + 1
c
      do 66 mmi=1,norbshi
        imm = norbi + mmi
        do 65 mmj=1,norbshj
          jmm = norbj + mmj
          arg = zero
          do 30 k=1,nk
C           arg = arg + dmat(imm,jmm,k)*cphi(k)
            arg = arg + dmat(1,imm,jmm,k)*cphi(1,k)
     $                - dmat(2,imm,jmm,k)*cphi(2,k)
   30     continue
c
          do 62 mi=1,norbshi
            im = norbwki + mi
            do 61 mj=1,norbshj
              jm= norbwkj + mj
              dmatwk(im,jm) = dmatwk(im,jm) + arg*
     $         dsym(mmi,mi,lorbi+1,isym)*dsym(mmj,mj,lorbj+1,isym)
   61       continue
   62     continue
c
   65   continue
   66 continue
c
c      note: there are cases, such as hcp(0001) and fcc(111)
c      surfaces, for which the point group operations are too
c      few to transform the "irreducible part" of the surface BZ
c      into the full SBZ. In these cases the slack is taken up
c      by the operation of time reversal, which need not take
c      atoms into atoms. Time reversal, however, requires that the
c      expression on the RHS of the preceding equation be replaced
c      by its real part. Since "dmatwk" has been defined as real in
c      this routine, the effects of time reversal are correctly
c      included here, and in the input data file, only the point
c      group operations of the film should be included.
c
c
            norbwkj = norbwkj + norbshj
            norbj = norbj + norbshj
   80     continue
          norbwki = norbwki + norbshi
          norbi = norbi + norbshi
   85   continue
c
 1000 continue
c
      RETURN
      END
