c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBOUT
c
c
      subroutine NEBOUT( iflpath,ea_neb,eb_neb,nimg_neb, pathout )
c---------------------------------------------------------------
c Purpose: output NEB data for post-processing
c
c Written: Renee M. Van Ginhoven, 11-May-2005, for 2.59
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Declare variables
      INTEGER iflpath,nimg_neb,img
      REAL*8 ea_neb, eb_neb
      REAL*8 totlen
c
      DIMENSION pathout(5,nimg_neb)
c
c Local declarations
      DATA zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      write(iflpath,*)'# Img, distance, energy, fdottan,'
      write(iflpath,3000) INT(zero), zero, ea_neb, zero
c
      do img = 1,nimg_neb
         write(iflpath,3000) img, 
     $        pathout(1,img), pathout(2,img), pathout(3,img),
     $        pathout(4,img), pathout(5,img)
      enddo
c
c Compute total path length
      totlen = pathout(1,nimg_neb) + pathout(4,nimg_neb)
      write(iflpath,3000)  nimg_neb+1, totlen, eb_neb, zero
 3000 FORMAT(i2,2x,f10.6,2x,f20.8,2x,3(f10.6,2x))
c
c
c
      RETURN
      END
