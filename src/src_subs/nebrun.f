c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBRUN
c
c
      subroutine NEBRUN( IWR,IWRNEB, inebstat, ierrg,
     $ do_neb, image,engytotl,
     $ igstep,igstop,nhistg, gconv,gblend, tstep,
     $ nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
     $ Lantikink_neb,Loptall_neb,Lclimb_neb,Lnoninteract,
     $ itypa,atmass, nafrc, iforce,
     $ nslics,islics,islicats,vslics,    iatframe,
     $ ndim, dolvatm, rprim,orig,rscale, nsym,nsyma, rmatsym, naofsym,
     $ natm, imgstat, engyneb, fatm, ratm,ratma,ratmb, wksma,
     $ ntmp, tmp, igrouplist, ngroup )
c---------------------------------------------------------------
c Purpose: manage a NEB run
c
c Written: Peter A. Schultz, 28-January-2002, for v2.54
c          Based on original code by Aidan Thompson (SNL)
c
c Revision history:
c  30Mar20-PAS/2.68e: error write upon geom fail
c  31Jul17-PAS/2.68: unified slics, added slinky layers; 68d fix update
c   5Aug15-PAS/2.67: handle error in geom updates
c  11Apr15-PAS/2.66: handle errors in images
c  20Feb09-PAS/2.62: pass through vgfix option
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  11May05-RMV/2.59: climb/noninteracting/post-analysis data; patched
c  10Feb05-PAS/2.58b: fix status for non-masters
c  25Apr04-PAS/2.58: purge "igtype" passing
c  11Sep03-PAS/2.57: remove itmaxd-dependence; neb+broy bugfix
c  10Jan03-PAS/2.55: l-BFGS into GRELAX
c  16Oct02-PAS/2.54c: fix g77 problem (cxg77)
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
c   Option flags for NEB:
      LOGICAL    Lantikink_neb, Loptall_neb
      LOGICAL    Lclimb_neb, Lnoninteract
c   Atomic info: type of each atom, mass for each type
      DIMENSION  itypa(natm), atmass(*)
c   Info about lattice, for internal/external coord converts
      LOGICAL    dolvatm
      DIMENSION  rprim(3,4),orig(3),rscale(3)
c   Info defining relaxation/frame constraints; and symmetry
      DIMENSION  iforce(natm)
cRMV SLIC: projected coordinate variables and arrays
      DIMENSION  islics(4,*),islicats(*),vslics(3,*)
      DIMENSION  iatframe(3)
c
      DIMENSION  rmatsym(*), naofsym(*)
c   Status and energies for images
      DIMENSION  imgstat(*), engyneb(*)
c   Forces, coords for this local image calculation
      DIMENSION  fatm(3,natm), ratm(3,natm)
c   Coordinates for the initial-point, and end-product
      DIMENSION  ratma(3,natm),ratmb(3,natm)
c   List of first node id in each group of nodes.
      DIMENSION  igrouplist(*)
c
c Scratch space:
      DIMENSION  wksma(natm,*)
      DIMENSION  tmp(ntmp)
c
c Local declarations:
      LOGICAL    set_done
      DATA       zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Find out who is what:
      call MPNODES_G( nodes )
      call MPNODE0_G( neb_master )
      call MPNODE_G( node )
      call MPCOMM_G( icomm )
c
c Carve out memory to be used:
      if( node.eq.neb_master )then
        nebatm =   natm*nimg_neb
      else
        nebatm =   natm
      endif
      i_ratm = 1
      i_fatm = i_ratm + nebatm*3
      i_tmp2 = i_fatm + nebatm*3
      ntmp2  = ntmp - 2*nebatm*3
      if( node.eq.neb_master )then
c       Carve out more temp memory
        i_nx    = i_tmp2
        i_wrlx  = i_nx    + nebatm
        i_grlx  = i_wrlx  + nebatm
        i_frlx  = i_grlx  + nebatm*3
        i_vrlx  = i_frlx  + nebatm*3
        i_mrlx  = i_vrlx  + nebatm*3
        i_trlx1 = i_mrlx  + nebatm*3
        i_trlx2 = i_trlx1 + nebatm*3
        i_tmp3  = i_trlx2 + nebatm*3
        ntmp3 = ntmp - i_tmp3 + 1
      endif
c
      if( node.eq.neb_master )then
c       Retrieve R for all, and F for finished, images
c
c       Turn off writing to output files:
        IWRX = -1
c
c       Retrieve coords from saved file
c
        write(IWRNEB,*) '>>>>> retrieve all image coords',
     $                   ' from neb_geom'
        call FLOPENA( iflnebg, 'neb_geom' )
c
        call NEBRGET( IWRX  ,iflnebg,
     $   neb_step, nimg_neb, ea_neb,eb_neb,
cxg77     $   natm, tmp(i_fatm),   ratma,ratmb,
     $   natm, imgstat(nimg_neb+1),   ratma,ratmb,
c -->          imgstat(nimg_neb) - do not stomp true imgstat (1)
     $   tmp(i_ratm) )
c -->    ratm_neb(3an)-o
c
c       Retrieve forces from saved file, if they exist
c        (also points code to right place in file to append new frc's).
c
        write(IWRNEB,*) '>>>>> retrieve previous NEB image results',
     $                   ' from neb_efrc'
        call FLOPENA( iflnebf, 'neb_efrc' )
        REWIND( unit=iflnebf )
c
        call NEBFGET( IWRX  , iflnebf, nimg_neb,
     $   natm, imgstat, engyneb,
     $   wksma,       tmp(i_fatm) )
c -->    fatmx(3a)-s, fatm_neb(3an)-o
c
      endif
c
c Collect image results from nodes (and save to file):
c
      write(IWR,*) '>>>>> collect image results at master NEB node'
      if( node.eq.neb_master )
     $write(IWRNEB,*) '>>>>> collect image results at master NEB node'
      call NEBREAP( IWRNEB, iflnebf, image, engytotl, ierror,
     $ node,neb_master,nodes,
     $ nimg_neb, natm,imgstat, engyneb, fatm,
     $ tmp(i_fatm),  ntmp2, tmp(i_tmp2), 
     $         igrouplist, ngroup, icomm )
c -->  fatm-neb(3an)-io     tmp-s
c
       if( node .eq. neb_master .and. ierror .ne. 0 )then
c        We have a problem - one of the images (at least) went bad
c        It is a job problem, because this job has to die, and we don't
c        really have a clean way to tell everyone without a hang because
c        the bad image may not be responsive.
         write(IWR   ,*)'***** ERROR: scf fail in image, NEB terminated'
         write(IWRNEB,*)'***** ERROR: scf fail in image, NEB terminated'
         call FLUSH( IWR )
         call FLUSH( IWRNEB )
         call MPBCASTI( neb_master, 1, ierror, icomm )
         call STOPXERR( 'ERROR in NEB images' )
       else
c        This is everyone (neb_master detected no error) inc. neb_mr
c         .OR. neb_mr detected error above and tells everyone else here
c              (everyone except failed image(s), which must be caught
c               in nebstat routines that notified neb_master of problem.
         call MPBCASTI( neb_master, 1, ierror, icomm )
         if( ierror .ne. 0 )then
c          Error detected, stop the job
           call STOPXERR( 'ERROR in NEB images' )
         endif
c
       endif
c NEB Master node manages and works the results:
c
c Set default completion status to "continue"
c NEB Master will evaluate status, and terminate processes.
      inebstat = 0
      if( node.eq.neb_master )then
c
c       Find if all image forces done:
        set_done = .true.
        do  img=1,nimg_neb
          if( imgstat(img).lt.2 ) set_done = .false.
        enddo
c
c       inebstat=0(continue),1(out_of_steps),2(converged)
        if( set_done )then
          write(IWRNEB,*) '***** All NEB image forces done *****'
c
c         Retrieve all image coordinates (already have the forces):
          IWRX = -1
c
c         Convert input images into internal cartesian bohr coords
          call NEBNTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $     nimg_neb, imgstat,  ratma, ratmb,
     $     tmp(i_ratm),      wksma(1, 1),wksma(1, 4) )
c -->      ratm_neb(3an)-io  ratm(3a)-s  ratm0(3a)-s
c
c         Modify computed forces with NEB forces
          write(IWRNEB,*) '>>>>> Apply NEB forces to images'
          call NEBFORCE( IWRNEB,  Lantikink_neb, Loptall_neb,
     $     Lclimb_neb, Lnoninteract,
     $     nimg_neb, spring_neb,
     $     natm, imgstat, engyneb, ea_neb, eb_neb,
     $     tmp(i_fatm), tmp(i_ratm), ratma, ratmb,
c -->      fatm_neb(3an)-io ratm_neb(3an)-i
     $     wksma(1,1), wksma(1,4), wksma(1,7), tmp(i_tmp2) )
c -->      rilo(3a)-s  rihi(3a)-s  tmpa(9a)-s  pathout-o
c          RMV: pathdat(5,nimg) gets path info for post-process
c          (ndat=5)
c
c         Write record of NEB-modified forces to NEB output file
c
          write(IWRNEB,*) '********** NEB modified forces **********'
          call NEBFPUT( IWRNEB, nimg_neb,
     $     natm, img_stat, engyneb,
     $     tmp(i_fatm) )
c -->      fatm_neb(3an)-io
          write(IWRNEB,*)
c
c         Write to post-NEB file (RMV): fsystem dot tangent, and
c          distance to previous and next images (contents of pathdat)
c
          call FLOPENA( iflpath, 'neb_post' )
          REWIND( unit = iflpath )
          write(iflpath,*)  '#write out path data for post-processing'
          call NEBOUT( iflpath,ea_neb,eb_neb,nimg_neb,
     $     tmp(i_tmp2) )
c -->      pathout(5,nimg)-i
          call FLCLOSE( iflpath )  
c
c         Apply the various constraints on the system forces
c
          write(IWRNEB,*) '>>>>> Apply constraints to NEB forces'
          n3atm = 3*natm
          img_off = 0
          nx_off = 0
          nrlx_off = 0
          nrlxdim = 0
          do 200 img=1,nimg_neb
c
c           Load full coordinates and forces into temporary locations
            call DCOPY( n3atm, tmp(i_ratm+img_off),1, wksma(1,10),1 )
            call DCOPY( n3atm, tmp(i_fatm+img_off),1, wksma(1,13),1 )
            img_off = img_off + n3atm
c
c           Fix frame of reference - project out frame-violating forces
            call FIXFRAME( natm, iatframe,
     $       wksma(1,10), wksma(1,13) )
c -->        ratm(3a)-io  fatm(3a)-io
c
c           Select out atoms to be relaxed
            call GSELECT( natm, nafrc, iforce, nrlxg,
     $       nslics,islics,islicats,vslics,
     $       wksma(1,19),wksma(1,22),
c -->        ratmP(3a)-o ratmP(3a)-s
     $       itypa, atmass,  natmg,wksma(1,16),
c -->                              nxatm(a)-o
     $       wksma(1,10),wksma(1,13),
c -->        ratm0(3a)-i frc(3a)-i
     $       wksma(1, 1), wksma(1, 4), wksma(1, 7) )
c -->        ratmg(3a)-o  frcg(3a)-o   amassg(3a)-o
c
c           Load the various arrays into the right places ...
c            ... number of free components
c                Be careful: packing ints, from real into real space!
            call ICOPYOFF( natmg, wksma(1,16),0, tmp(i_nx),nx_off )
            nx_off = nx_off + natmg
            natmrlx = nx_off
c            ... coords, forces, and masses for each degree of freedom
            call DCOPY( nrlxg, wksma(1, 1),1, tmp(i_grlx+nrlx_off),1 )
            call DCOPY( nrlxg, wksma(1, 4),1, tmp(i_frlx+nrlx_off),1 )
            call DCOPY( nrlxg, wksma(1, 7),1, tmp(i_mrlx+nrlx_off),1 )
            nrlx_off = nrlx_off + nrlxg
            nrlxdim = nrlx_off
  200     continue
c
c         Check convergence over entire set of images simultaneously
          call CONVCHK( IWRNEB, 'nebfrc', igstat, igstep,igstop, gconv,
     $     chgmax,chgrms,  3,natmrlx, tmp(i_frlx) )
c -->                                 frcg(3an)-io
          inebstat = igstat
c
          if( inebstat.ne.2 )then
c           Not converged, propagate geoms for next NEB step
c
c           Propagate new configurations for next NEB step:
c
            write(IWRNEB,*) '>>>>> Blend to get new image coordinates'
            engypass = zero
            call GRELAX( IWRNEB, gblend, tstep, engypass,
     $       igstep,nhistg, nrlxdim, chgmax,
     $       natmrlx,tmp(i_wrlx),  tmp(i_nx),
c -->                watmg(an)-io  nxatm(an)-i
     $       tmp(i_grlx),  tmp(i_frlx), tmp(i_vrlx), tmp(i_mrlx),
c -->        ratmg(3an)-io frcg(3an)-i  vatmg(3an)-s amassg(3an)-i
     $       tmp(i_trlx1), tmp(i_trlx2), ntmp3,tmp(i_tmp3) )
c -->        wk3(3an)-s    wk4(3an)-s
c
c           Update master coordinate vectors, impose/check constraints
c
            write(IWRNEB,*)'>>>>> Update image coords for next NEB step'
            nrlx_off = 0
            img_off = 0
            do  img=1,nimg_neb
              call DCOPY( nrlxg, tmp(i_grlx+nrlx_off),1, wksma(1, 1),1 )
              call DCOPY( n3atm, tmp(i_ratm+img_off) ,1, wksma(1,13),1 )
c
              call GUPDATE( ierrg, wksma(1,13),
c -->                              ratm0(3a)-io
     $         natm, iforce, ndim,nsym,nsyma, rprim, rmatsym, naofsym,
     $         nslics,islics,islicats,vslics,  iatframe,
     $         wksma(1,19),
c -->          ratmP(3a)-i
     $         wksma(1, 1), wksma(1, 4),wksma(1, 7), wksma(1,10) )
c -->          ratmg(3a)-is ratmx(3a)-s dratm0(3a)-s dratmsy(3a)-s
c
c             Put updated coordinates into master coordinate vector
              call DCOPY( n3atm, wksma(1,13),1, tmp(i_ratm+img_off),1 )
              nrlx_off = nrlx_off + nrlxg
              img_off = img_off + n3atm
c
              if( ierrg .ne. 0 )then
                inebstat = -1
                goto 313
              endif
              imgstat(img) = 1
            enddo
c
c           Convert images into external user frame
            call NEBXTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $       nimg_neb, imgstat,  ratma, ratmb,
     $       tmp(i_ratm),      wksma(1, 1),wksma(1, 4) )
c -->        ratm_neb(3an)-io  ratm(3a)-s  ratm0(3a)-s
c
c           Reset force file to null results to reflect the new step
            call NEBFWIPE( iflnebf )
c
            write(IWRNEB,*) '>>>>> record new image coords to file'
            neb_step = igstep
c
c           Record new coords to NEB status file ...
            REWIND( unit=iflnebg )
            call NEBRSAV( iflnebg,
     $       neb_step, nimg_neb, ea_neb,eb_neb,
     $       natm, imgstat, ratma,ratmb,
     $       tmp(i_ratm) )
c -->        ratm_neb(3an)-io
c
c            ... and write new band coords to the NEB output file:
            call NEBRSAV( IWRNEB,
     $       neb_step, nimg_neb, ea_neb,eb_neb,
     $       natm, imgstat, ratma,ratmb,
     $       tmp(i_ratm) )
c -->        ratm_neb(3an)-io
c
c           End NEB set-done, propagated new images for next NEB step
          endif
c
c         End NEB set-done
        else
          write(IWRNEB,*) '>>>>> NEB step not done, do more images'
        endif
c
  313   continue
        call FLCLOSE( iflnebg )
        call FLCLOSE( iflnebf )
c
cxxx    Should have a status/sync here?
        if( inebstat.eq.0 )then
c         Send new image assignments
c
          write(IWRNEB,*) '>>>>> Assign NEB image(s)'
          call NEBASSGN( IWRNEB, image,neb_step, ierror,
     $     node,neb_master,nodes,
     $     nimg_neb, natm,imgstat, ratm,
     $     tmp(i_ratm),  ntmp2, tmp(i_tmp2), 
     $         igrouplist, ngroup, icomm )
c-->       ratm_neb(3an)-io     tmp-s
c
        else
          if( inebstat .lt. 0 )then
             write(IWR,*)    ' ***** ERROR in NEB geoms'
             write(IWRNEB,*) ' ***** ERROR in NEB geoms'
          endif
c         Send null assignments to terminate processes
c
          write(IWRNEB,*) '>>>>> Terminate NEB, null assignments'
          call NEBASTOP( IWRNEB, node,neb_master,nodes,
     $     ntmp2,tmp(i_tmp2), 
     $         igrouplist, ngroup, icomm )
c
        endif
c
c       End node.eq.neb_master branch
      else
c       Start node.ne.neb_master (implying nodes>1)
c
cxxx    Should have status/sync here?
c       Receive new image assignments
        call NEBASSGN( IWR, image,neb_step, ierror,
     $   node,neb_master,nodes,
     $   nimg_neb, natm,imgstat, ratm,
     $   tmp(i_ratm),  ntmp2, tmp(i_tmp2), 
     $        igrouplist, ngroup, icomm )
c-->     ratm_neb(3an)-io     tmp-s
c
c       Set step number to neb step number:
        igstep = neb_step
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
