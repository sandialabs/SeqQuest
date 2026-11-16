c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBREAP
c
c
      subroutine NEBREAP( IWR   , iflnebf, image, engytotl, ierror,
     $ node,neb_master,nodes,
     $ nimg_neb, natm,imgstat, engyneb, fatm,
     $ fatm_neb,  ntmp, tmp, igrouplist, ngroup, icomm )
c---------------------------------------------------------------
c Purpose: collect results from NEB images
c
c Written: Peter A. Schultz, 28-January-2002, for v2.54
c
c Revision history:
c  10Apr15-PAS/2.66: nebstat propagated for scf fail in image
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output quantities
      DIMENSION  imgstat(*)
      DIMENSION  engyneb(*)
      DIMENSION  fatm(3,natm)
c Scratch space:
      DIMENSION  fatm_neb(3,natm,*)
      DIMENSION  tmp(*)
c Local:
c    status value for image with completed force
c Input: List of first node id in each group of nodes.
      DIMENSION igrouplist(*)
c
      DATA  imgstdone / 2 /
c
c >>>> EXECUTABLE CODE:
c
      ierror = 0
c
      if( node.eq.neb_master )then
        ineberr = 0
c
        if( nodes.eq.1 )then
c         Serial calculation
c
c         Find next undone image:
          do 100 img=1,nimg_neb
            if( imgstat(img).eq.1 )then
c             Found unfinished image, check against input ...
              if( img.ne.image )then
                write(IWR,*) 'Completed image number =',image
                write(IWR,*) 'Apparent assigned number=',img
                write(IWR,'(a,20i2)') 'imgstat=',
     $           (imgstat(i),i=1,nimg_neb)
                call STOPXERR( 'NEBREAP: image.ne.image_neb' )
              endif
              eimage = engytotl
c              ... save to NEB status file, and to NEB output:
              call NEBPUTF( IWR    ,  image,eimage,  natm, fatm )
              call NEBPUTF( iflnebf,  image,eimage,  natm, fatm )
c              ... and transfer coordinates into work coords:
              engyneb(image) = eimage
              call DCOPY( 3*natm, fatm,1, fatm_neb(1,1,image),1 )
              goto 190
            endif
  100     continue
          ierror = -1
          call STOPXERR( 'NEBREAP: no unfinished images to reap' )
c
  190     continue
          write(IWR   ,*) '>>>>> Collected NEB image#',image
          imgstat(image) = imgstdone
c
c         End branch #nodes=1
        else
c         #nodes>1
c         Idiot-parallel NEB -> retrieve image results from nodes
c
          img0 = 0
          do 300 igroup=1,ngroup
c
c           Find next undone image:
            do 200 img=img0+1,nimg_neb
              if( imgstat(img).eq.1 )then
                imgtodo = img
                goto 210
              endif
  200       continue
c           Did not find an image needing work for this node
            goto 301
c
  210       continue
            img0 = imgtodo
c
            nodei = igrouplist(igroup)
            if( nodei.eq.neb_master )then
c             This image stayed on master NEB node
c
              if( image.ne.imgtodo )
     $         call STOPXERR( 'NEBREAP/neb_master: image.ne.imgtodo' )
              imgi = image
              eimage = engytotl
              call DCOPY( 3*natm, fatm,1, fatm_neb(1,1,imgi),1 )
c
            else
c             This image was done elsewhere ...
c
c             Receive/unpack assignment {image,E[img],F[img]} results
              call NEBSTAT( nodei, image, ierr, icomm )
              if( ierr .ne. 0 )then
c               This image failed for some reason, report error
                write(IWR,*)'***** ERROR: remote image failed=',img0
                write(IWR,*)'>>>>> Error code reported=',ierr,' (scf)'
                ineberr = ineberr + 1
c               Skip rest of processing this image
                goto 300
              endif
              call REAPGET( nodei, imgi, eimage,
     $         3*natm, fatm_neb(1,1,imgtodo), ntmp, tmp, icomm )
c
              if( imgi.ne.imgtodo )then
 9020           format(a,2i6)
                write(IWR   ,9020) 'Image assigned and image returned',
     $           ' disagree on node#',nodei
                write(IWR   ,9020) 'Image sent, returned=',imgtodo,imgi
                call STOPXERR( 'NEBREAP: image.ne.image_neb from nodei')
              endif
c
            endif
            write(IWR   ,'(a,i4,a,i6)')
     $       '>>>>> NEB image#',imgi,' results from node#',nodei
            imgstat(imgi) = imgstdone
c
            engyneb(imgi) = eimage
            call NEBPUTF( IWR    , imgi,eimage,natm,fatm_neb(1,1,imgi) )
            call NEBPUTF( iflnebf, imgi,eimage,natm,fatm_neb(1,1,imgi) )
  300     continue
  301     continue
c
c         End branch of #nodes>1.
        endif
c
        ierror = ineberr
c       End branch node=neb_master
      else
c       node.ne.nebmaster (and, by implication, #nodes>1)
c
c     Check if this node is the root node for the group
         call MPNODE( iproc )
         call MPNODE0( master )
         if( iproc.eq.master )then
c
c       Pack/send image result {image,E[image],F[image]} to NEB master
            call NEBSTATOK( neb_master, image, icomm )
            call REAPPUT( neb_master, image, engytotl,
     $           3*natm, fatm, ntmp, tmp, icomm )
         endif
c
      endif
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
c End
