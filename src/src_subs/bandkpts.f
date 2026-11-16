c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDKPTS
c
c
      subroutine BANDKPTS( veckn,veck,wtkn,nnscf,nknscfd,nkd,bvl,
     $    nbranch, nbranchd,bsymb,rprim,gprim,dk,nbr,IWR )
c---------------------------------------------------------------
c Purpose: defined bands of k-points for band structure
c
c Written: A.H. Edwards, 100608,
c
c Revision history:
c  21Jan15-PAS/2.66: add-gamma option, add-klist option
c  03May12-ACP: Merged into 2.64
c---------------------------------------------------------------
c
cAHE 100608
c     We start to put substance into this program
c     We start by fixing the number of points on each branch (=15)
c     we will use bvl to determine the appropriate subroutine.
c     This is now almost a dummy program. We will put the first scf veck 
c     into veckn.
c     nnscf will be the number of non-scf k-points
c     We must check at the end of the routine that nnscf <nknscfd
c
c Band structure symbol definitions are taken from:
c   W. Setyawan and S. Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
c
      IMPLICIT REAL*8(a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0  )
      DIMENSION  veckn(3,nknscfd),veck(3,nkd),wtkn(nknscfd)
      DIMENSION  rprim(3,3)
      DIMENSION  gprim(3,3),dk(nknscfd),nbr(nbranchd)
      CHARACTER*3  bvl
      CHARACTER*2  bsymb(2,nbranchd)
c
      CHARACTER  symb1, symb2
c
c >>>> EXECUTABLE CODE:
c
c  We need to calculate the total number of non-SCF points.
c
      write(IWR,*)  'Bravais lattice, bvl = ', bvl
c
      if( bvl .eq. 'gam' )then
c        Add just the gamma point to the wanted k-points
         write(IWR,*) '  Adding gamma point, and just the gamma point'
         nnscf = 1
         veckn(1,1) = zero
         veckn(1,2) = zero
         veckn(1,3) = zero
         wtkn(1) = one
c
      elseif( bvl .eq. 'bcc' )then
c       Body-centered cubic lattice (per Setyawan and Curtarolo)
        call BANDSBCC( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'sc ' )then
c       Simple cubic lattice (standard)
        call BANDSSC( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
      elseif( bvl .eq. 'cub' )then
c       Cubic lattice  (per Setyawan and Curtarolo)(different implementation of sc)
        call BANDSCUB( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'fcc' )then
c       Face-centered cubic lattice (standard)
        call BANDSFCC( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
      elseif( bvl .eq. 'fc2' )then
c       Face-centered cubic  (per Setyawan and Curtarolo) (equivalent to fcc)
        call BANDSFCC2( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'hex' )then
c       Hexagonal lattice
        call BANDSHEX( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'orc' )then
c       Orthorhombic  (per Setyawan and Curtarolo)
        call BANDSORC( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'rhl' )then
c       Rhombohedral  (per Setyawan and Curtarolo)
        call BANDSRHL( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,rprim,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'tet' )then
c       Tetragonal  (per Setyawan and Curtarolo)
        call BANDSTET( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'tri' )then
c       Triclinic  (per Setyawan and Curtarolo)
        call BANDSTRI( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'stc' )then
cxxx264 ??? Simple-triclinic (standard) ???
        call BANDSSTC( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $        nknscfd,gprim,dk,nbr,IWR )
c
      elseif( bvl .eq. 'mcl' )then
c       Simple monoclinic (per Setyawan and Curtarolo)
        call BANDSMCL( nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $       nknscfd,rprim,gprim,dk,nbr,IWR )
      elseif( bvl .eq. 'smc' )then
cxxx264 ???  Simple monoclinic = redundant with bvl='mcl' above?
c       We need to use the real-space vectors to determine:
c       phi, the acute monoclinic angle
c       lambda = c/(b*cos(phi)
c
        write(IWR,*) 'Monoclinic >smc< is a work in progress!'
        call STOPXERR( 'Bravais lattice >smc< not supported' )
c
      else
        write(IWR,*) '***** ERROR: bravais option not supported:', bvl
        call STOPXERR( 'Bravais lattice option not supported' )
      endif
c
      write(IWR,*) 'BANDKPTS: total nnscf,max = ', nnscf, nknscfd
      if( nnscf.gt.nknscfd ) call STOPXERR('nnscf is too large!')
c
c    That's all Folks
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
c End:
