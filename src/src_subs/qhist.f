c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QHIST
c
c
      subroutine QHIST0( istatfl, statkey,
     $ ndim0,natm0,ntyp0, itypa, natmnm,typnm,atmnm )
c---------------------------------------------------------------
c Purpose: write out data as running log of code that can be
c          be parsed for post-processing.
c
c Written: Peter A. Schultz, 19-May-2001, for v2.55
c
c Revision history:
c  29Apr09-PAS/2.62: install r_units into fmax (thx ACPineda)
c  21Jan08-PAS/2.61: spin-optimization
c   3Jan08-PAS/2.61: enable skip if file-unit disabled
c  14Jun05-PAS/2.59: File units
c  18Mar04-PAS/2.57d: fix cell energy output
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER  statkey*(*)
c
      DIMENSION  itypa(*)
      CHARACTER*(*)  typnm(*),atmnm(*)
c
      DIMENSION  vdata(3,*)
c
c local declarations:
      DIMENSION  vec3(3)
      CHARACTER  wantlbl*4
      CHARACTER  keychar*1
      DATA  keychar / '@' /
c
      CHARACTER*4  dist
      DATA         dist / 'ANGS' /
      DATA         one     / 1.d0 /
      DATA         au_to_A / 0.5291771d0 /
c
c >>>> EXECUTABLE CODE:
c
 9100 format(a)
 9103 format(5x,3f20.8)
c
c Check that this feature is on ...
      if( istatfl .lt. 0 ) goto 1999
c
      wantlbl = statkey(1:4)
c
      if( wantlbl .eq. 'INIT' )then
c       Initialize system information in file:
c
        REWIND( unit=istatfl )
c
        write(istatfl,9100)  keychar // '=KEYCHAR'
c
        write(istatfl,9100)  keychar // 'DISTANCE UNIT'
        if( dist.eq.'BOHR' )then
          write(istatfl, '(1x, a)' )     'BOHR'
        elseif( dist.eq.'ANGS' )then
          write(istatfl, '(1x, a)' )     'ANGSTROM'
        else
          call FLGETIERR( IERRFL )
          write(IERRFL,*) 'QHIST0: unknown distance unit: ',dist
          call STOPXERR( 'QHIST0 - unknown distance unit' )
        endif
c
        ndim = ndim0
        natm = natm0
        ntyp = ntyp0
        write(istatfl,9100)  keychar // 'DIMENSION'
        write(istatfl, '(1x,i4)' )     ndim
c
        write(istatfl,9100)  keychar // 'NUMBER OF ATOMS'
        write(istatfl, '(1x,i8)' )    natm
        do  iatm0=0,natm-1,20
          iatm1 = iatm0 + 1
          iatms = iatm0 + 20
          if( natm.lt. (iatm0+20) ) iatms = natm
          write(istatfl,'(20i4)')  (itypa(iatm),iatm=iatm1,iatms)
        enddo
c
        write(istatfl,9100)  keychar // 'TYPES'
        write(istatfl, '(1x,i8)' )    ntyp
        do  ityp=1,ntyp
          write(istatfl,'(a)')  typnm(ityp)
        enddo
c
      endif
c
 1999 continue
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QHIST
c
c
      entry QHIST( istatfl, statkey, ndata, rdata, vdata )
c
c Check that this feature is on ...
      if( istatfl .lt. 0 ) goto 2999
c
      if( dist .eq. 'ANGS' )then
        r_units = au_to_A
      elseif( dist .eq. 'BOHR' )then
        r_units = one
      else
        call FLGETIERR( IERRFL )
        write(IERRFL,*) 'QHIST: unknown distance unit: ',dist
        call STOPXERR( 'Unknown distance unit in QHIST' )
      endif
c
      wantlbl = statkey(1:4)
      if( wantlbl .eq. 'IMAG' )then
c       Write out current NEB image number:
        write(istatfl,9100)  keychar // 'IMAGE'
        write(istatfl,'(1x,i8)')  ndata
c
      elseif( wantlbl .eq. 'UCST' )then
c       Write out current cell step number:
        write(istatfl,9100)  keychar // 'UCSTEP'
        write(istatfl,'(1x,i8)')  ndata
c
      elseif( wantlbl .eq. 'CELL' )then
c       Write out current cell vectors:
        write(istatfl,9100)  keychar // 'CELL VECTORS'
        do  idim=1,3
          vec3(1) = vdata(1,idim) * r_units
          vec3(2) = vdata(2,idim) * r_units
          vec3(3) = vdata(3,idim) * r_units
          write(istatfl,9103)  vec3
        enddo
c
      elseif( wantlbl .eq. 'GSTE' )then
c       Write out current geometry step number:
        write(istatfl,9100)  keychar // 'GSTEP'
        write(istatfl,'(1x,i8)')  ndata
c
      elseif( wantlbl .eq. 'COOR' )then
c       Write out current geometry information
        write(istatfl,9100)  keychar // 'COORDINATES'
        natm = ndata
        do  iatm=1,natm
          vec3(1) = vdata(1,iatm) * r_units
          vec3(2) = vdata(2,iatm) * r_units
          vec3(3) = vdata(3,iatm) * r_units
          write(istatfl,9103)  vec3
        enddo
c
      elseif( wantlbl .eq. 'ESPI' )then
c       Write out converged SCF energy
        write(istatfl,9100)  keychar // 'ESPINSCF(Ry): E, e(up), e(dn)'
        write(istatfl,9103)  rdata, vdata(1,1),vdata(2,1)
c
      elseif( wantlbl .eq. 'SP-F' )then
c       Write out new spin polarization
        write(istatfl,9100)  keychar // 'SPIN-POLARIZATION (e)'
        write(istatfl,9103)  rdata
c
      elseif( wantlbl .eq. 'ESCF' )then
c       Write out converged SCF energy
        write(istatfl,9100)  keychar // 'ESCF'
        write(istatfl,9103)  rdata
c
      elseif( wantlbl .eq. 'FORC' )then
c       Write out forces
        write(istatfl,9100)  keychar // 'FORCES'
        natm = ndata
        do  iatm=1,natm
          vec3(1) = vdata(1,iatm) / r_units
          vec3(2) = vdata(2,iatm) / r_units
          vec3(3) = vdata(3,iatm) / r_units
          write(istatfl,9103)  vec3
        enddo
c
      elseif( wantlbl .eq. 'VELO' )then
c       Write out velocities
        write(istatfl,9100)  keychar // 'VELOCITY'
        natm = ndata
        do  iatm=1,natm
          vec3(1) = vdata(1,iatm) * r_units
          vec3(2) = vdata(2,iatm) * r_units
          vec3(3) = vdata(3,iatm) * r_units
          write(istatfl,9103)  vec3
        enddo
c
      elseif( wantlbl .eq. 'FMAX' )then
c       Write out maximum force
        write(istatfl,9100)  keychar // 'FMAX'
        write(istatfl,9103)  rdata / r_units
c
      elseif( wantlbl .eq. 'FRMS' )then
c       Write out rms force
        write(istatfl,9100)  keychar // 'FRMS'
        write(istatfl,9103)  rdata / r_units
c
      elseif( wantlbl .eq. 'ECEL' )then
c       Write out converged geometry energy
        write(istatfl,9100)  keychar // 'ECELL'
        write(istatfl,9103)  rdata
c
      elseif( wantlbl .eq. 'STRE' )then
c       Write out stress tensor
        write(istatfl,9100)  keychar // 'STRESS'
        ndim = 3
        do  idim=1,ndim
          vec3(1) = vdata(1,idim)
          vec3(2) = vdata(2,idim)
          vec3(3) = vdata(3,idim)
          write(istatfl,9103)  (vec3(i),i=1,ndim)
        enddo
c
      elseif( wantlbl .eq. 'MAXS' )then
c       Write out maximum stress
        write(istatfl,9100)  keychar // 'MAXSTRESS'
        write(istatfl,9103)  rdata
c
      elseif( wantlbl .eq. 'RMSS' )then
c       Write out rms stress
        write(istatfl,9100)  keychar // 'RMSSTRESS'
        write(istatfl,9103)  rdata
c
      elseif( wantlbl .eq. 'END ' )then
c       Close the record
        write(istatfl,9100)  keychar // 'END FILE'
c
      else
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '***** This qhist key unknown:',statkey
        call STOPXERR( 'QHIST: unknown keyword' )
c
      endif
c
cpas: machine dependent flush. - I should replace this (utl)
      call FLUSH( istatfl )
c
 2999 continue
c
      RETURN
c
c    That's all Folks!
c
      END
