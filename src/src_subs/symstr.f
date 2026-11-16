c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMSTR
c
c
      subroutine SYMSTR( ndim, str1, nsym, rmatsym )
c---------------------------------------------------------------
c Purpose: symmetrize stress tensor
c
c Written: Peter A. Schultz, 28-September-1999, v2.40(3dstress)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output arrays:
      DIMENSION  str1(3,3),str2(3,3),str3(3,3)
      DIMENSION  rmatsym(3,4,*)
c local declarations
      DIMENSION  vi(3),svi(3)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( nsym.le.0 ) RETURN
c
c Load identity:
      do  jd=1,3
        do  id=1,3
          str2(id,jd) = str1(id,jd)
        enddo
      enddo
c
c  Cycle over non-ident syms, and accumulate transformed tensors:
c
      do 300 isym=1,nsym
c
        do  jdim=1,3
          vi(1) = str1(1,jdim)
          vi(2) = str1(2,jdim)
          vi(3) = str1(3,jdim)
c
          svi(1) = rmatsym(1,1,isym)*vi(1) + rmatsym(2,1,isym)*vi(2)
     $           + rmatsym(3,1,isym)*vi(3)
          svi(2) = rmatsym(1,2,isym)*vi(1) + rmatsym(2,2,isym)*vi(2)
     $           + rmatsym(3,2,isym)*vi(3)
          svi(3) = rmatsym(1,3,isym)*vi(1) + rmatsym(2,3,isym)*vi(2)
     $           + rmatsym(3,3,isym)*vi(3)
c
          str3(1,jdim) = svi(1)
          str3(2,jdim) = svi(2)
          str3(3,jdim) = svi(3)
        enddo
c
        do  idim=1,3
          vi(1) = str3(idim,1)
          vi(2) = str3(idim,2)
          vi(3) = str3(idim,3)
c
          svi(1) = rmatsym(1,1,isym)*vi(1) + rmatsym(2,1,isym)*vi(2)
     $           + rmatsym(3,1,isym)*vi(3)
          svi(2) = rmatsym(1,2,isym)*vi(1) + rmatsym(2,2,isym)*vi(2)
     $           + rmatsym(3,2,isym)*vi(3)
          svi(3) = rmatsym(1,3,isym)*vi(1) + rmatsym(2,3,isym)*vi(2)
     $           + rmatsym(3,3,isym)*vi(3)
c
          str2(idim,1) = str2(idim,1) + svi(1)
          str2(idim,2) = str2(idim,2) + svi(2)
          str2(idim,3) = str2(idim,3) + svi(3)
        enddo
c
  300 continue
c
c Load symmetrized stress tensor back in original space:
      symscal = one/ DBLE( nsym+1 )
      do  jd=1,3
        do  id=1,3
          str1(id,jd) = symscal*str2(id,jd)
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
