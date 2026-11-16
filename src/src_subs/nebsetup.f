c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEBSETUP
c
c
      subroutine NEBSETUP( IDAT,IWR, IWRNEB,
     $ do_neb, image, neb_step,
     $ nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
     $ Lantikink_neb,Loptall_neb,Lclimb_neb,Lnoninteract,
     $ ndim, dolvatm, rprim,orig,rscale,
     $ natm, imgstat, engyneb, ratm,ratma,ratmb, wksma,
     $ ntmp, tmp, igrouplist, ngroup )
c---------------------------------------------------------------
c Purpose: Initialize and set up a Nudged Elastic Band (NEB)
c          series of calculations.
c
c Written: Peter A. Schultz, 17-January-2002, for v2.54
c
c Revision history:
c  19Nov08-PAS/2.62: mpi-comm rework to accommodate k-parallel
c  16Nov07-PAS/2.60b: bugfix - mprecvi lenmsg argument
c  13Jun07-PAS/2.60: rename MP calls
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  11May05-RMV/2.59: add flags: climb, noninteract
c  10Jul03-PAS/2.56: change master NEB write to IWRNEB
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Identify version/revision status for this NEB
      CHARACTER  version*5,         lastup*11
      DATA       version /'1.10 '/, lastup/'19Nov08-PAS'/
c
      LOGICAL    do_neb
c
c Output: option flags for NEB
      LOGICAL    Lantikink_neb, Loptall_neb
      LOGICAL    Lclimb_neb, Lnoninteract
c
c Input: info about lattice, for internal/external coord converts
      LOGICAL    dolvatm
      DIMENSION  rprim(3,4),orig(3),rscale(3)
c Output: coords for local calculation; initial, final images
      DIMENSION  ratm(3,natm)
      DIMENSION  ratma(3,natm),ratmb(3,natm)
c         status of each image, and computed energy for each
      DIMENSION  imgstat(*), engyneb(*)
c Output: List of first node id in each group of nodes.
c Array dimension is ngroup+1, where igrouplist(ngroup+1) 
c is set to number of nodes for convenience.
      DIMENSION igrouplist(*)
c
c Scratch:
      DIMENSION  wksma(natm,*)
      DIMENSION  tmp(ntmp)
      DIMENSION  itmpbuf(2)
c
c Local declarations:
      LOGICAL    Lrestart_neb
      CHARACTER    label*72, readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
c
c >>>> EXECUTABLE CODE:
c
 9000 format(
     $  / 4x,'Program SEQQUEST: LCAO-PS-DFT electronic structure code'
     $  / 4x,'Written by Peter A. Schultz' / )
 9001 format( 1x,72('-')
     $  / 4x,'***** NEB module for transition state finding *****'
     $  / 4x,'Written by Aidan Thompson, adapted by PAS for SeqQuest'
     $  / 4x,'Based on the NEB method of Hannes Jonnson'
     $ // 4x,'NEB Version: ',a,'        Last revised: ',a
     $  / 4x,'Maximum number of NEB images =',i6
     $  / 1x,72('-') / )
c
c Assign workspace in tmp()
      i_ratm = 1
c
c Find out who is what:
      call MPNODES_G( nodes )
      call MPNODE0_G( neb_master )
      call MPNODE_G( node )
      call MPCOMM_G( icomm_g )
c
      if( node.eq.neb_master )then
c
c       >>> Open NEB output file, and give some vitals:
c
        write(IWR,*) '***** NEB Module invoked *****'
        write(IWR,9001)  version,lastup,nimg_max
c
        write(IWR,*) 'This node is the NEB master'
        write(IWR,*) 'NEB output file= neb_out'
        call FLOPENA( IWRNEB, 'neb_out' )
        if( IWRNEB.ne.IWR )then
          write(IWRNEB,9000)
        endif
        write(IWRNEB,9001)  version,lastup,   nimg_max
c
c       >>> Initialize NEB parameters:
c
        call DAT0NEB( nimg_neb, nimg_max,
     $   ea_neb,eb_neb, spring_neb,
     $   Lantikink_neb,Loptall_neb, 
     $   Lclimb_neb, Lnoninteract, Lrestart_neb,
     $   imgstat, natm, ratm, ratma )
c
c       >>> Input data to define this NEB calculation:
c
        write(IWR,*) '>>>>> reading NEB input data'
c
        call NEBDAT( IWRNEB,IDAT,
     $   neb_step, nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
     $   Lantikink_neb,Loptall_neb,
     $   Lclimb_neb, Lnoninteract, Lrestart_neb,
     $   natm, imgstat, ratma,ratmb,
     $   tmp(i_ratm) )
c-->     ratm_neb(3*natm*nimg_neb)-o
c
        n_ratm = 3*natm*nimg_neb
        i_tmp2 = i_ratm + n_ratm
        ntmp2 = ntmp - n_ratm
c
c       >>> If restart, recover coordinates
c       >>> and recover energy/forces for previously completed images
cpas:       I should consider making NEB restart automatic
c
        if( Lrestart_neb )then
c
c         Load preexisting image geometries from NEB configuration file
c
          write(IWRNEB,*) '>>>>> get restart coords from neb_geom'
          call FLOPENA( iflnebg, 'neb_geom' )
          call NEBRGET( IWRNEB,iflnebg,
     $     neb_step, nimg_neb, ea_neb,eb_neb,
     $     natm, imgstat, ratma,ratmb,
     $     tmp(i_ratm) )
c -->      ratm_neb(3an)-io
          call FLCLOSE( iflnebg )
c
c         Identify and load energies/forces for previously done images
c
          write(IWRNEB,*) '>>>>> restart image E/F from neb_efrc'
          call FLOPENA( iflnebf, 'neb_efrc' )
          call NEBFSTAT( IWRNEB, iflnebf, nimg_neb,
     $     natm, imgstat,
     $     wksma )
c -->      fatmx(3a)-s
          call FLCLOSE( iflnebf )
c
        endif
c
c       Convert input images into internal cartesian bohr coords
        call NEBNTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $   nimg_neb, imgstat,  ratma, ratmb,
     $   tmp(i_ratm),      wksma(1,1),wksma(1,4) )
c -->    ratm_neb(3an)-io  ratm(3a)-s ratm0(3a)-s
c
c       Fill in coordinates for empty images in the band:
        write(IWRNEB,*) '>>>>> fill in empty images in band'
        call NEBFILL( IWRNEB,
     $   natm, nimg_neb, imgstat, ratma,ratmb,
     $   tmp(i_ratm),      wksma )
c-->     ratm_neb(3an)-io  ratmx(6*natm)-s
c
c       Convert all images back to external user coordinates
        call NEBXTRNL( ndim,natm,dolvatm,rprim,orig,rscale,
     $   nimg_neb, imgstat,  ratma, ratmb,
     $   tmp(i_ratm),      wksma(1,1),wksma(1,4) )
c -->    ratm_neb(3an)-io  ratm(3a)-s ratm0(3a)-s
c
        if( .not.  Lrestart_neb )then
c         Save record of all image coordinates to NEB restart file
          call FLOPENA( iflnebg, 'neb_geom' )
          REWIND( unit=iflnebg )
          call NEBRSAV( iflnebg,
     $     neb_step, nimg_neb, ea_neb,eb_neb,
     $     natm, imgstat, ratma,ratmb,
     $     tmp(i_ratm) )
          call FLCLOSE( iflnebg )
        endif
c
c       Write record of all image coordinates to NEB output file
        call NEBRSAV( IWRNEB,
     $   neb_step, nimg_neb, ea_neb,eb_neb,
     $   natm, imgstat, ratma,ratmb,
     $   tmp(i_ratm) )
c
        write(IWRNEB,*) '>>>>> NEBSETUP: Assign images to nodes'
c
        IWRX = IWRNEB
c       End node=neb_master branch
      else
c       Not node=neb_master, receive assignments (image#, coords).
c
c       Assign memory (i_ratm not used, but to be safe):
        n_ratm = 3*natm
        i_tmp2 = i_ratm + n_ratm
        ntmp2 = ntmp - n_ratm
        IWRX = IWR
c
      endif
c
c     Set up partitioning of nodes into image-specific communicators
c      and record in igrouplist()
c
c     There may be a fancier way to assign processors,
c     to get better Scalapack performance.
c
      if( node .eq. neb_master )then
c
        if( nodes .ge. nimg_neb )then
c         More nodes than images: group nodes for each image  
          ngroup = nimg_neb
          nper = nodes / ngroup
          nremain = nodes - nper*ngroup
          inode = 0
          do  igroup=1,ngroup
            igrouplist(igroup) = inode
            if( igroup .le. nremain )then
              inode = inode + nper + 1
            else
              inode = inode + nper
            endif
          enddo
          igrouplist(ngroup+1) = nodes
        else
c         More images than nodes, then make each node into its own group
          ngroup = nodes
          do  igroup=1,ngroup
            igrouplist(igroup) = igroup - 1
          enddo
          igrouplist(ngroup+1) = nodes
        endif
c
c       Tell each node which group it belongs to
c
        do  igroup=1,ngroup
          do  nodei=igrouplist(igroup),igrouplist(igroup+1)-1
            if( nodei .eq. neb_master )then
              mygroup = igroup
            else
              itmpbuf(1) = igroup
              lenmsg = 1
              call MPSENDI( nodei, lenmsg, itmpbuf, icomm_g )
            endif
          enddo
        enddo
c
      else
c
c       Each node receives its group assignment from the neb_master
c
        lenmsg = 1
        call MPRECVI( neb_master, lenmsg, itmpbuf, icomm_g )
        mygroup = itmpbuf(1)
c
      endif
c     
c Split processors up into local communicators:
c
      call MPSPLITGCOMM( mygroup )
      write(IWR,'(1x,a,i4,a)') 
     $     'MPINFO::: NEB Group = ',mygroup  ,'           :::'
c
c At this point, I have coordinates and completion status
c for every image in the band.  Now cycle over groups of 
c nodes and assign images (with coordinates) to groups.
c
cxxx  I should probably put a status/sync step here.
      call NEBASSGN( IWRX, image,neb_step, ierror,
     $ node,neb_master,nodes,
     $ nimg_neb, natm,imgstat, ratm,
     $ tmp(i_ratm), ntmp2, tmp(i_tmp2) ,
     $ igrouplist, ngroup, icomm_g )
c-->   ratm_neb(3an)-io    tmp-s
c
c    That's all Folks!
c
      RETURN
      END
