c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBASSGN
c
c
      subroutine NEBASSGN( IWR, image,neb_step, ierror,
     $ node,neb_master,nodes,
     $ nimg_neb, natm,imgstat, ratm,
     $ ratm_neb,  ntmp, tmp, igrouplist, ngroup, icomm )
c---------------------------------------------------------------
c Purpose: assign NEB images to nodes
c
c Written: Peter A. Schultz, 21-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c  10Feb05-PAS/2.58b: fix NEBASTOP/neb_step misuse
c  11Nov02-PAS/2.54d: better stop
c---------------------------------------------------------------
c
c At this point, I have coordinates and completion status
c - in imgstat() - for each image in the NEB.  This routine
c assigns images to groups of nodes, sending each node in 
c each group the coordinates (and identity) for the image 
c it is assigned.
c Must anticipate the following possibilities:
c (0) This documentation was blindly extended to from nodes 
c     to groups of nodes. It may no longer be correct.     
c (1) Status is: all images complete (ie. with E and forces)
c     Noone has any work to do ... must update geometries.
c     This should have been done elsewhere -> error if happens
c (2) Each group gets an image to work on (the default).
c     Easy.
c (3) Some groups are busy through an NEB step ...
c     Easy.
c (4)  ... and some idle groups never get work to do.
c     This is non-trivial, because some nodes will get stopped
c     in the nebsetup while others will get through to the
c     nebupdate phase.  They will need to be made to sync up,
c     either to start the next NEB step, or at finish of NEB.
c
c     The available groups of nodes are stored in the array
c     igrouplist, giving the global node id of the first node
c     in each group. ngroup is the number of groups. nroup is
c     assumed to never exceed the number of images, although
c     there may be reasons to want to allow this.
c
c 
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output (coordinates for scf calculation)
      DIMENSION  ratm(3,natm)
c Input quantities
      DIMENSION  ratm_neb(3,natm,*)
      DIMENSION  imgstat(*)
c Input: List of first node id in each group of nodes.
      DIMENSION igrouplist(*)
c Scratch space:
      DIMENSION  tmp(*)
c
c >>>> EXECUTABLE CODE:
c
      ierror = 0
c
      if( node.eq.neb_master )then
c
        image = 0
        if( nodes.eq.1 )then
c         Serial calculation
c
c         Find next undone image:
          do 100 img=1,nimg_neb
            if( imgstat(img).eq.1 )then
c             Found free image ...
              image = img
c              ... and transfer coordinates into work coords:
              call DCOPY( 3*natm, ratm_neb(1,1,image),1, ratm,1 )
              goto 190
            endif
  100     continue
          ierror = 1
          call STOPXERR( 'NEBASSGN: no unfinished images to assign' )
c
  190     write(IWR,*) '>>>>> Working on NEB image#',image
c
c         End branch #nodes=1
        else
c         Non-Idiot-parallel NEB, #nodes>1 -> assign images to groups of nodes
c
          img0 = 0
          do 300 igroup = 1,ngroup
c
c           Find next undone image:
            do 200 img=img0+1,nimg_neb
              if( imgstat(img).eq.1 )then
                imgtodo = img
                goto 210
              endif
  200       continue
c           Did not find an image needing work for nodei
            goto 301
c
  210       img0 = imgtodo
c
            do nodei = igrouplist(igroup),igrouplist(igroup+1)-1
            if( nodei.eq.neb_master )then
c             This image assignment stays on neb_master node
              image = imgtodo
              call DCOPY( 3*natm, ratm_neb(1,1,imgtodo),1, ratm,1 )
            else
c             Pack and ship off image {image,R[image]} to nodei
              call ASSGNPUT( node, nodei, neb_master, imgtodo,neb_step,
     $         3*natm, ratm_neb(1,1,imgtodo), ntmp, tmp, icomm )
            endif
c
            write(IWR,'(a,i4,a,i6)')
     $       '>>>>> NEB image#',imgtodo,' assigned to node#',nodei
            end do
c
  300     continue
  301     continue
c
c         End branch of #nodes>1.
        endif
c
c       End branch node=neb_master
      else
c       node.ne.neb_master (and, by implication, #nodes>1)
c
c       Receive and unpack assignment {image,R[image]} from neb_master
        call ASSGNGET( node, neb_master, neb_master, image,neb_step,
     $   3*natm, ratm, ntmp, tmp, icomm )
c
c       If this node got a null assignment, signal to stop:
        if( image.eq.0 )then
          write(IWR,*) '>>>>> Null NEB assignment at node#',node
          call STOPXOK( 'end NEB' )
        endif
c
        write(IWR,'(a,i4,a,i6)')
     $   '>>>>> NEB image#',image  ,' assigned to this node#',node
c
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBASTOP
c
c
      entry      NEBASTOP( IWR, node,neb_master,nodes, ntmp, tmp, 
     $     igrouplist, ngroup, icomm )
c---------------------------------------------------------------
c Purpose: send null assignments to nodes to terminate NEB
c
c Written: Peter A. Schultz, 29-January-2002, for v2.54
c
c Revision history:
c  10Feb05-PAS/2.58b: fix NEBASTOP/neb_step misuse
c---------------------------------------------------------------
c
      if( node.eq.neb_master )then
c
        write(IWR,*) '>>>> Terminating NEB'
        if( nodes.eq.1 )then
c         Serial calculation: #nodes=1 -> do nothing
c
c         End branch #nodes=1
        else
c         Idiot-parallel NEB: #nodes>1 -> send null assignments
c
          i_data = 1
          n_data = 1
          i_tmp2 = i_data + n_data
          ntmp2 = ntmp - n_data
          lendata = 0
c
          do 600 igroup=1,ngroup
            do nodei = igrouplist(igroup),igrouplist(igroup+1)-1
               node_to = nodei
c     **** using "node_to=0,nodes-1" loop triggers a compiler
c     **** error (seg-fault).  Using an intermediate variable
c     **** prevent compiler mistake.  6Oct02-PAS
               imgtodo = 0
               neb_stop = 0
               if( node_to.ne.node )then
c     Pack up and ship off null assignment
                  call ASSGNPUT( node,node_to,neb_master, 
     $                 imgtodo,neb_stop,
     $                 lendata, tmp(i_data), ntmp, tmp(i_tmp2), icomm )
               endif
            end do
 600     continue
c
c         End #nodes>1
        endif
c
c       End node=neb_master
      else
c       node.ne.neb_master -> not allowed here.
        call STOPXERR( 'NEBASTOP: node.ne.neb_master' )
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
