c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> REAPPUT
c
c
      subroutine REAPPUT( node_put, image, engytotl,
     $ ndata, data_img, ntmp, tmp, icomm )
c---------------------------------------------------------------
c Purpose: Put out results of assignment to master node
c
c Written: Peter A. Schultz, 28-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output (data for image)
      DIMENSION  data_img(*)
c Scratch space (buffer space for send/receive)
      DIMENSION  tmp(*)
c Local:
      DATA       zero / 0.d0 /
c converting int to real, add a little so back-conversion is clean:
      DATA       xtra / 0.1d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( image.gt.0 )then
c       Pack up the image id, and the data that goes along
        rimage = DBLE( image ) + xtra
        tmp(1) = rimage
        tmp(2) = engytotl
        call DCOPY( ndata, data_img,1, tmp(3),1 )
        lenmsg = 1 + 1 + ndata
      else
c       We did not have an assignment to work on and return
        rimage = xtra
        tmp(1) = rimage
        lenmsg = 1
      endif
c
      call MPSENDR8( node_put, lenmsg, tmp, icomm )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> REAPGET
c
c
      entry      REAPGET( node_get, image, engytotl,
     $ ndata, data_img, ntmp, tmp, icomm )
c---------------------------------------------------------------
c Purpose: Get results of assignment worked on
c
c Written: Peter A. Schultz, 28-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
      lenmsg = ndata + 2
      call MPRECVR8( node_get, lenmsg, tmp, icomm )
c
      rimage = tmp(1)
      image = rimage
c
      if( image.gt.0 )then
c       Unload data
        if( ndata.ne. (lenmsg-2) ) call STOPXERR( 'REAPGET: lenmsg' )
        engytotl = tmp(2)
        call DCOPY( ndata, tmp(3),1, data_img,1 )
      else
c       Nothing to get ... this node was idle
      endif
c
c    That's all Folks!
c
      RETURN
      END
