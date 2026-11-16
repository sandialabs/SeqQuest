c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> UCATCHK
c
c
      subroutine UCATCHK( IWR, ndim,natm, dolvatm, itypa,znuc,
     $ rprim,orig,rscale,
     $ ratm0,ratm )
c---------------------------------------------------------------
c Purpose: project atom coordinates into the primary unit cell
c
c Written: Peter A. Schultz, 24-September-1997, v2.21
c
c Revision history:
c  29Oct09-PAS/2.62: enforce error when atoms too close; cosmetics
c  23Nov04-PAS/2.58a: fix proximity checks (thx Andy Pineda)
c  20Oct03-PAS?2.58: proximity checks for atoms
c  21Jun01-PAS/2.48: replace STOPs
c   7Jul00-PAS/2.45: unified scaling
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  ratm0() = coordinates of atoms as input
c  ratm()  = coordinates of atoms inside primary unit cell
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input:
      LOGICAL    dolvatm
      DIMENSION  itypa(natm),znuc(*)
      DIMENSION  rprim(3,4),orig(3),rscale(3)
c Input (raw) atomic coordinates
      DIMENSION  ratm0(3,natm)
c Output atomic coordinates (projected into unit cell)
      DIMENSION  ratm(3,natm)
c
c Local declarations:
      DIMENSION  tr(3),ttr(3)
      LOGICAL  outuc,shift
c Set warning at 1/3 less than H-H molecule, kill at less than 1/2 H-H bond
      DATA     vclose, tooclose / 0.4999, 0.2999d0 /
      DATA     lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
      ioutuc = 0
      ijflag = 0
c
      do 100 iatm=1,natm
c       Project atoms into primary unit cell in periodic directions
        call PROJUC( outuc,shift, ndim, ratm0(1,iatm),tr, rprim )
c
        if( outuc )then
c         Atom is outside uc in non-periodic direction - fatal
          ioutuc = ioutuc + 1
c         Get coordinates in original input units/frame of reference
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     ttr, ratm0(1,iatm) )
          write(IWR,*) '***** ERROR: Atom outside unit cell:',iatm
          write(IWR,9010) ttr
 9010     format('>>>>> Coordinates:',3f12.8)
c
        elseif( shift )then
c         Atom is outside uc in periodic direction - safe to shift
          ratm(1,iatm) = tr(1)
          ratm(2,iatm) = tr(2)
          ratm(3,iatm) = tr(3)
c
c         Get coordinates in original input units/frame of reference
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     ttr,tr )
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     tr, ratm0(1,iatm) )
c
          if( lstout.gt.1 )then
            write(IWR,*) '***** NOTE: atom outside unit cell:',iatm
            write(IWR,9011)  tr,ttr
 9011       format(6x,'Old coords:',3f12.8/6x,'Shifted to:',3f12.8)
          endif
c
        else
c         Atom safely within primary unit cell to begin with, copy
          ratm(1,iatm) = ratm0(1,iatm)
          ratm(2,iatm) = ratm0(2,iatm)
          ratm(3,iatm) = ratm0(3,iatm)
        endif
c
c       Check for fatal overlap of atoms
c
  913   format(1x,a,' R(#',i4,'--#',i4,')=',f12.8)
        ityp = itypa(iatm)
        if( znuc(ityp) .ne. zero )then
c
          nda = 0
          ndb = 0
          ndc = 0
          if( ndim.gt.0 ) nda = 1
          if( ndim.gt.1 ) ndb = 1
          if( ndim.gt.2 ) ndc = 1
          do  jatm=1,iatm-1
            jtyp = itypa(jatm)
c
            if( znuc(jtyp) .ne. zero )then
              rijmin = 9999.
c
              do  idc=-ndc,ndc
                do  idb=-ndb,ndb
                  do  ida=-nda,nda
                    do  ir=1,3
                      tr(ir) = ratm(ir,iatm) - ratm(ir,jatm)
     $                      - DBLE(ida)*rprim(ir,1)
     $                      - DBLE(idb)*rprim(ir,2)
     $                      - DBLE(idc)*rprim(ir,3)
                    enddo
                    rij = SQRT( tr(1)**2 + tr(2)**2 + tr(3)**2 )
                    if( rij .lt. rijmin ) rijmin = rij
                  enddo
                enddo
              enddo
              if( rijmin .lt. vclose )then
                if( ijflag .eq. 0 )then
c                 First instance, give generic warning
                  write(IWR,*) '>>>>> WARNING: atoms very close'
                endif
c               Assess seriousness, give specific notice
                if( rijmin .lt. tooclose )then
                  ijflag = 2
                  write(IWR,913) '***** ERROR: ',iatm,jatm, rijmin
                else
                  ijflag = MAX( ijflag, 1 )
                  write(IWR,913) '>>>>> CAUTION: ',iatm,jatm, rijmin
                endif
              endif
c
            endif
c
          enddo
c
        endif
c
c       Close loop over atoms
  100 continue
c
      if( ijflag.gt.1 ) call STOPXERR( 'atoms on top of each other' )
      if( ioutuc.gt.0 ) call STOPXERR( 'uc-atom - atom outside cell' )
c
c    That's all Folks!
c
      RETURN
      END
