c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDAVE
c
c
      subroutine GRIDAVE(IWR,n1,n2,n3,nd,arr,ave,plave,avelbl)
c---------------------------------------------------------------
c Purpose: compute planar averages of quantities on 3D grid
c
c Written: Peter A. Schultz
c
c Revision history:
c  26Jan11-PAS/2.62: cosmetic mods for output
c  22Oct98-PAS/2.27: communicate global averages
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c  Variables:
c   IWR      - output unit number (IWR.le.0 turns output off)
c   n1,n2,n3 - dimensions of grid
c   nd       - dimension of arr(), 1:real ; 2:complex
c   arr()    - input array for which averages are desired
c   ave()    - global averages
c   plave()  - output array of planar averages
c   avelbl   - arbitrary length character string describing average
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  arr(nd,n1*n2,n3), ave(nd),plave(nd,n3)
      CHARACTER*(*)  avelbl
c local declarations:
      LOGICAL  printon
      DATA  zero / 0.d0 /
c
      printon = (IWR.gt.0)
c
      if( printon )  write(IWR,9000)  avelbl
 9000 format(1x,'PLANAR GRID AVERAGE OF ',a)
 9100 format(1x,'Global average:',1pd16.8,
     $       1x,'Planar averages:'
     $       / (1x,1p5d15.7) )
c
      do 1000 id=1,nd
c
        ave(id) = zero
c
        do  i3=1,n3
c
          sum = zero
          do  i12=1,n1*n2
            sum = sum + arr(id,i12,i3)
          enddo
c
          sum = sum/ DBLE(n1*n2)
          plave(id,i3) = sum
          ave(id) = ave(id) + sum
c
        enddo
c
        ave(id) = ave(id)/ DBLE(n3)
c
        if( printon ) write(IWR,9100)  ave(id),(plave(id,i3),i3=1,n3)
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
