c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GSTAMP
c
c
      subroutine GSTAMP( IWR, igstep,igeom,igstat,image,
     $ iatmfmt, ntyp,natm, natmnm, itypa, atmnm,typnm, ratm )
c---------------------------------------------------------------
c Purpose: output a complete formatted geometry to a file
c
c Written: Peter A. Schultz, 14-February-2002, for v2.54
c
c Revision history:
c  12Apr15-PAS/2.66: master-pnly write of external file
c   8Nov02-PAS/2.54d: add geom file to output.
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  itypa(*)
      CHARACTER*(*)  typnm, atmnm
      DIMENSION  typnm(*), atmnm(*)
      DIMENSION  ratm(3,natm)
c
c Local:
      CHARACTER*10  imgname
c
c >>>> EXECUTABLE CODE:
c
 9430 format(1x,79('#') )
c
      write(IWR,*)
      write(IWR,9430)
c
      call MPNODE( iproc )
      call MPNODE0( master )
c
      if( image.eq.0 )then
c       Standard SCF geometry relaxation
c
        if( igstat.lt.0 )then
          write(IWR,9431) 'INITIAL',igstep,igeom+1
        elseif( igstat.lt.2 )then
          write(IWR,9431) '  NEW',igstep,igeom+1
        elseif( igstat.ge.2 )then
          write(IWR,9431) 'FINAL',igstep,  igeom
        endif
 9431   format(1x,a,' GEOMETRY for step',i4,', of',i4,' total')
c
        if( iproc.eq.master )then
c         Open geom file with default name:
          call FLOPENA( igeomfl, 'geom' )
          REWIND( unit=igeomfl )
c
          write(igeomfl,'(a,i7)') 'atom, type, position;  step#',igstep
          write(IWR    ,'(a,i7)') 'atom, type, position;  step#',igstep
        endif
c
      else
c       Multiple-image (NEB) calculation.
c
        if( igstat.lt.0 )then
          write(IWR,9432) 'INITIAL',image, igstep
        elseif( igstat.lt.2 )then
          write(IWR,9432) '  NEW',image, igstep
        elseif( igstat.ge.2 )then
          write(IWR,9432) 'FINAL',image, igstep
        endif
 9432   format(1x,a,' GEOMETRY for image',i3,', step',i4 )
c
        if( iproc.eq.master )then
c         Open geom file with image-specific geom file name
          call STRNUM( nstr, imgname, image,image )
          call FLOPENA( igeomfl, 'geom'//imgname(1:nstr) )
          REWIND( unit=igeomfl )
c
          write(igeomfl,'(a,i4,1x,i7)')
     $     'atom, type, position; image#',image
          write(IWR    ,'(a,i4,1x,i7)')
     $     'atom, type, position; image#',image
        endif
c
      endif
c
c Write latest geometry to user output listing file:
      call CONFIGWR( IWR    , iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c
      if( iproc.eq.master )then
c       Write latest geometry to running geometry file:
        call CONFIGWR( igeomfl, iatmfmt, natmnm, ntyp,typnm,
     $   natm, atmnm, itypa, ratm )
c
        call FLCLOSE( igeomfl )
      endif
c
      write(IWR,9430)
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
