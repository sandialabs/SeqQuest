c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CELLRUPD
c
c
      subroutine CELLRUPD( do_cell, iucstat, dolvatm,
     $ ndim, rprim, oldprim, orig,rscale,
     $ natm,ratm,ratm0 )
c---------------------------------------------------------------
c Purpose: update atoms, k-points, etc. to new unit cell
c
c Written: Peter A. Schultz, 2-May-2002, v2p53 ( cell opt)
c
c Revision history:
c   8Aug03-PAS/2.56: fix k-vector handling, remove gprim stuff
c  27Aug02-PAS/2.53g: sign error in rt-hand cell to gprim
c   1May02-RPM: adapted from Rick Muller cell code in main program
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0 )
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      LOGICAL    do_cell, dolvatm
      DIMENSION  rprim(3,3), oldprim(3,3)
      DIMENSION  orig(3),rscale(3)
      DIMENSION  ratm(3,natm), ratm0(3,natm)
c Local:
      DIMENSION  ratmlv(3)
      DATA       lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( iucstat.gt.0 ) RETURN
c
      do  id=1,ndim
        if( ABS( rscale(id) - one ) .gt. 1.d-8 )
     $   call STOPXERR( 'scale+cell - incompatible with cell opt' )
      enddo
c
      if( .not. dolvatm )then
C        call STOPXERR( 'cell opt with lattice coords only' )
c
c       Update the atom positions based on the new lattice vectors
c ****  Only necessary if we do not have lattice units, which we
c ****  currently require.  However, the trial code below had
c ****  crashed, and it is blocked for now.
c
        if( lstout.gt.3 )then
          call FLGETIWR( IWR )
          write(IWR,*) 'CELLRUPD: old rprim'
          do  jd=1,3
             write(IWR,'(3f12.6)')  (oldprim(id,jd),id=1,3)
          enddo
          write(IWR,*) 'CELLRUPD: new rprim'
          do  jd=1,3
             write(IWR,'(3f12.6)')  (rprim(id,jd),id=1,3)
          enddo
          write(IWR,*) 'CELLRUPD: ratm positions before update'
          do  iatm=1,natm
            write(IWR,'(3f12.6)')  (ratm(i,iatm),i=1,3)
          enddo
        endif
c
c       Convert atomic coordinates to internal (Cartesian bohr)
        call RINTRNL( ndim,natm, dolvatm, oldprim,orig,rscale,
     $   ratm,ratm0 )
c
c       The following really should be limited by problem dimension,
c       and should also include transform of orig(), and origws().
        do  iatm=1,natm
c         Project on old lattice vectors to get lattice units:
          call LVPROJ( oldprim, ratm0(1,iatm),ratmlv )
c         Expand in new lattice vectors to get get new positions
          call LVXPND( rprim, ratm0(1,iatm),ratmlv )
        enddo
c
c       Convert atomic positions back to user external frame
        call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $   ratm,ratm0 )
c
        if( lstout.gt.3 )then
          write(IWR,*) 'CELLRUPD: Atom positions after update'
          do  iatm=1,natm
            write(IWR,'(3f12.6)')  (ratm(i,iatm),i=1,3)
          enddo
          write(IWR,*) 'CELLRUPD: ratm0 after update'
          do  iatm=1,natm
            write(IWR,'(3f12.6)')  (ratm0(i,iatm),i=1,3)
          enddo
          write(IWR,*) 'CELLRUPD: new rprim'
          do  jd=1,3
            write(IWR,'(3f12.6)')  (rprim(id,jd),id=1,3)
          enddo
        endif
c
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
