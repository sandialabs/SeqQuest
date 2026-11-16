c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICINIT
c
c
      subroutine SLICINIT( sliclbl )
c---------------------------------------------------------------
c Purpose: initialize SLIC stuff
c
c Written: Peter A. Schultz, 14-August-2017
c---------------------------------------------------------------
c
      IMPLICIT NONE
      CHARACTER  sliclbl*(*)
      INTEGER  islictyp, ipslic
      REAL*8   pslic
c
      call SLICTSET( sliclbl, islictyp )
      call SLICPSET( sliclbl, pslic,ipslic )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICTSET
c
c
      subroutine SLICTSET( sliclbl, islictyp )
c---------------------------------------------------------------
c Purpose: set up flags for geometry update method
c
c Written: Peter A. Schultz, 30-October-2001, for v2.50 (dmd)
c---------------------------------------------------------------
c
c  Setting igc_type (communicated via common) for a given slic
c Relaxation schemes:
c  igc_fixed = frozen in place
c  igc_free  = free cartesian
c  igc_slic  = slic, 1-fold (spherical) or 2-fold (cylindrical)
c   igt_slic1 = 1-fold atom slic
c   igt_slic2 = 2-fold bridge atom slic
c  igc_frame = frame constraint
c  igc_bond  = bond constraint
c  igc_vgp   = vector-based constraint- in-plane
c  igc_vgd   = vector-based constraint- distance along axis
c  igc_vga   = vector-based constraint- along axis 
c  igc_vgl   = vector-based constraint- along axis, and distance
c  igc_slink = slinky/layer (slab inter-layer) 
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c  Geometry constraint/local-coord (SLICS) types
      COMMON  /SLICTYPS/ igc_fixed, igc_free,
     $ igc_slic,  igt_slic1,igt_slic2,
     $ igc_frame,igc_bond,igc_vgp,igc_vgd,igc_vga,igc_vgl,igc_slink
c
c input declarations:
      CHARACTER*(*)  sliclbl
c local declarations:
      CHARACTER*6    slictyp
c
c >>>> EXECUTABLE CODE:
c
c
      slictyp = sliclbl(1:6)
      if( slictyp(1:4).eq.'init' .or. slictyp(1:4).eq.'INIT' )then
c       Simply initialize the tags
        islictyp = -999
        igc_fixed = 0
        igc_free  = 1
        igc_slic  = 2
        igc_frame = 3
        igc_bond  = 4
        igc_vgp   = 5
        igc_vgd   = 6
        igc_vga   = 7
        igc_vgl   = 8
        igc_slink = 9
        igc_dev   = 10
c
        igt_slic1 = 1
        igt_slic2 = 2
c
      elseif( slictyp.eq.'FIXED ' )then
        islictyp = igc_fixed
      elseif( slictyp.eq.'FREE  ' )then
        islictyp = igc_free
      elseif( slictyp.eq.'SLIC  ' )then
        islictyp = igc_slic
       elseif( slictyp.eq.'SLIC1 ' )then
         islictyp = igt_slic1
       elseif( slictyp.eq.'SLIC2 ' )then
         islictyp = igt_slic2
      elseif( slictyp.eq.'FRAME ' )then
        islictyp = igc_frame
      elseif( slictyp.eq.'BOND  ' )then
        islictyp = igc_bond
      elseif( slictyp.eq.'VGFIXP' )then
        islictyp = igc_vgp
      elseif( slictyp.eq.'VGFIXD' )then
        islictyp = igc_vgd
      elseif( slictyp.eq.'VGFIXA' )then
        islictyp = igc_vga
      elseif( slictyp.eq.'VGFIXL' )then
        islictyp = igc_vgl
      elseif( slictyp.eq.'SLINKY' )then
        islictyp = igc_slink
      elseif( slictyp.eq.'SLDEV ' )then
        islictyp = igc_dev
c
      else
        write(IWR,*) 'SLICTYP not extended to more than initialize'
        call STOPXERR( 'No such slictyp found' )
      endif
c
c    That's all Folks!
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICTGET
c
c
      entry      SLICTGET( sliclbl, islictyp )
c
      if(     islictyp .eq. igc_fixed )then
        slictyp =  'FIX   '
      elseif( islictyp .eq. igc_free  )then
        slictyp =  'FREE  '
      elseif( islictyp .eq. igc_slic  )then
        slictyp =  'SLIC  '
       elseif( islictyp .eq. igt_slic1 )then
         slictyp =  'SL1   '
       elseif( islictyp .eq. igt_slic2 )then
         slictyp =  'SL2   '
      elseif( islictyp .eq. igc_frame )then
        slictyp =  'FRAM  '
      elseif( islictyp .eq. igc_bond  )then
        slictyp =  'BOND  '
      elseif( islictyp .eq. igc_vgp   )then
        slictyp =  'VGP   '
      elseif( islictyp .eq. igc_vgd   )then
        slictyp =  'VGD   '
      elseif( islictyp .eq. igc_vga   )then
        slictyp =  'VGA   '
      elseif( islictyp .eq. igc_vgl   )then
        slictyp =  'VGL   '
      elseif( islictyp .eq. igc_slink )then
        slictyp =  'LYR   '
      elseif( islictyp .eq. igc_dev   )then
        slictyp =  'DEV   '
c
      else
        call STOPXERR( 'Unknown SLIC typ - coding error' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICPSET
c
c
      subroutine SLICPSET( sliclbl, pslic,ipslic )
c---------------------------------------------------------------
c Purpose: set up flags for geometry update method
c
c Written: Peter A. Schultz, 30-October-2001, for v2.50 (dmd)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero = 0.d0, one=1.d0 )
      CHARACTER   sliclbl*(*)
c
      COMMON  / SLICPARS / sl_pars(9), isl_fank
c
      if(    sliclbl(1:4).eq.'INIT' )then
c       Initialize everything to default of zero
        do  i=1,9
          sl_pars(i) = zero
        enddo
        isl_fank = 2
c
c     Record desired parameter into 
      elseif( sliclbl(1:6).eq.'PAR1  ' )then
        sl_pars(1) = pslic
      elseif( sliclbl(1:6).eq.'PAR2  ' )then
        sl_pars(2) = pslic
      elseif( sliclbl(1:6).eq.'PAR3  ' )then
        sl_pars(3) = pslic
      elseif( sliclbl(1:6).eq.'PAR4  ' )then
        sl_pars(4) = pslic
c
      elseif( sliclbl(1:6).eq.'FANK0 ' )then
c       Leave sphere anchor force unchanged
        isl_fank = 0
      elseif( sliclbl(1:6).eq.'FANKF ' )then
c       Propagate full anchor force into anchor force
        isl_fank = 1
      elseif( sliclbl(1:6).eq.'FANKB ' )then
c       Propagate bond-projected force into anchor force
        isl_fank = 2
c
      else
        write(*,*) 'slicpset: sliclbl=',sliclbl
        call MTBUFF
        call STOPXERR( 'Unknown slicpset label' )
      endif

      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICPGET
c
      entry      SLICPGET( sliclbl, pslic,ipslic  )
c
c     Retrieve desired parameter into 
      if(     sliclbl(1:6).eq.'PAR1  ' )then
        pslic = sl_pars(1)
      elseif( sliclbl(1:6).eq.'PAR2  ' )then
        pslic = sl_pars(2)
      elseif( sliclbl(1:6).eq.'PAR3  ' )then
        pslic = sl_pars(3)
      elseif( sliclbl(1:6).eq.'PAR4  ' )then
        pslic = sl_pars(4)
c
      elseif( sliclbl(1:4).eq.'FANK'   )then
        ipslic = isl_fank
c
      else
        call STOPXERR( 'Bad slicpget argument' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLICSUM
c
c
      subroutine SLICSUM( IWR, natm,
     $ nslics,islics,islicats,vslics,    iatframe )
c---------------------------------------------------------------
c Purpose: summarize slic inputs, to output file
c Written: 31-July-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  islics(4,*),islicats(*),vslics(3,*)
c
      DIMENSION  iatframe(3)
c
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      if( iatframe(1).gt.0 ) write(IWR,*) 'DEV/FRAME activated'
c
      if( nslics.eq.0 ) RETURN
c
c  Recap all unified SLIC descriptors
c
      if( iprint.gt.0 .and. nslics.ne.0 )
     $  write(IWR,*) 'SLIC recap: # of slics nslics=',nslics
c
      if( iprint.lt.2 ) RETURN
c
      do  islic=1,nslics
        islictyp = islics(1,islic)
        write(IWR,'(a,i2,a,i3)') 'SLICS #',islic,' typ=',islictyp
c
        if( islictyp.eq.1 .or. islictyp.eq.2 )then
c         SLIC project atom
          iatmproj = islics(2,islic)
c**** DEV should reverse CYL/SPH order to 1-fold/2-fold
          write(IWR,'(2x,a,i3)') '  CYL=1,SPH=2, type=',islictyp
          write(IWR,'(2x,a,i8)') '  SLIC projected atom: ',iatmproj
          write(IWR,'(2x,a,2i8)')
     $     '       anchors: ',islics(3,islic),islics(4,islic)
c
        elseif( islictyp.eq.3 )then
c         FRAME
          write(IWR,*) 'DEV/FRAME activated, not unified into SLICS'
c
        elseif( islictyp.eq.4 )then
c         BOND atoms
          write(IWR,'(2x,a,2i6)') 'BOND length fixed between atoms=',
     $                            islics(2,islic),islics(3,islic)
c
        elseif( islictyp.ge.5 .and. islictyp.le.8 )then
c         VGFIX atom
          islicv = islics(4,islic)
          write(IWR,'(2x,a,2i3,3f12.4)')'  VGFIX: typ, #atoms, vector=',
     $     islictyp, islics(2,islic), (vslics(i,islicv),i=1,3)
c
        elseif( islictyp.eq.9 )then
c         SLINKY/LAYERS atom
          islicv = islics(4,islic)
          write(IWR,'(2x,a,i3,3f12.4)') '  SLINKY: #atoms, vector=',
     $               islics(2,islic), (vslics(i,islicv),i=1,3)
c
        endif
c
      enddo
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
