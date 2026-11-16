c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTSTR
c
c
      subroutine DFTSTR( ndim, strdxc, nptr,weight,
     $ dxc, drho )
c---------------------------------------------------------------
c Purpose: evaluate grid integral for 1st derivative GGA stress
c
c Written: Peter A. Schultz, 5-April-2002, for v2.53 (gga stress)
c
c Revision history:
c  10Oct03-PAS/2.57: accomodate spin-dependent GGA stress
c---------------------------------------------------------------
c
c Notes:
c   Evaluate the grid integral for GGA stress evaulation:
c
c     stress(i,j) = INT( dxc(r) * drh[i](r) * drh[j](r) )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0 )
c
c Input declarations:
      DIMENSION  dxc(nptr)
      DIMENSION  drho(nptr,3)
c Output stuff:
      DIMENSION  strdxc(3,3)
c
c >>>> EXECUTABLE CODE:
c
c If molecule, no stress to evaluate
      if( ndim .eq. 0 ) RETURN
c
      do  jd=1,ndim
        do  id=1,jd
          strij = zero
          do  nr=1,nptr
            strij = strij + dxc(nr)*drho(nr,id)*drho(nr,jd)
          enddo
          strdxc(id,jd) = strdxc(id,jd) + weight*strij
        enddo
      enddo
c
c Load the flip side
      do  jd=2,ndim-1
        do  id=1,jd-1
          strdxc(jd,id) = strdxc(id,jd)
        enddo
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
