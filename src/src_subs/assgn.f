c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ASSGNPUT
c
c
      subroutine ASSGNPUT( node, node_img, neb_master, image,imgstep,
     $ ndata, data_img, ntmp, tmp, icomm )
c---------------------------------------------------------------
c Purpose: Send out assignment of image data to work on
c
c Written: Peter A. Schultz, 21-January-2002, for v2.54
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
c       Pack up the image id, step, number, and data that goes along
c
        rimage = DBLE( image ) + xtra
        tmp(1) = rimage
        rstep  = DBLE( imgstep ) + xtra
        tmp(2) = rstep
        i_data = 3
        call DCOPY( ndata, data_img,1, tmp(i_data),1 )
        lenmsg = 2 + ndata
c
      else
c       No image is assigned, construct null assignment
c
        rimage = xtra
        tmp(1) = rimage
        lenmsg = 1
c
      endif
c
      call MPSENDR8( node_img, lenmsg, tmp, icomm )
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ASSGNGET
c
c
      entry      ASSGNGET( node, node_img, neb_master, image,imgstep,
     $ ndata, data_img, ntmp, tmp, icomm )
c---------------------------------------------------------------
c Purpose: Get assignment of data to work on
c
c Written: Peter A. Schultz, 21-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
      lenmsg = ndata + 2
c
      call MPRECVR8( neb_master, lenmsg, tmp, icomm )
c
      rimage = tmp(1)
      image = rimage
c
      if( image.gt.0 )then
c       Unpack data received for assigned image
c
        if( ndata.ne. (lenmsg-2) ) call STOPXERR( 'ASSGNGET: lenmsg' )
        imgstep = tmp(2)
        i_data = 3
        call DCOPY( ndata, tmp(i_data),1, data_img,1 )
c
      else
c       Null assignment ... this node is done - stop
c
        if( node.ne.neb_master ) call STOPX
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
