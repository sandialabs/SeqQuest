c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDDATA
c
c
      subroutine BANDDATA( IDAT,IWR, nbranchd,nbranch,bvl,bsymb )
c---------------------------------------------------------------
c Purpose: read band structure data from input deck.
c
c Written: A.H. Edwards, 6-Oct-2008
c
c Revisions:
c  21Jan15/2.66-PAS: add add-gamma option; lnwrite lines
c  03May12-ACP&29Dec12-PAS: Merged into 2.64
c---------------------------------------------------------------
c
c  Arrays/variables for band structure calculations
c
c     bvl = symbol for Bravais Lattice
c            fcc    face centered cubic
c            hex    hexagonal
c            stc    simple triclinic
c            sc     simple cubic
c            ...    there are others, see bandkpts
c           gam     Gamma point, in whatever lattice
c     nbranch = number of branches for band structure
c     bsymb(nbranchd) = symbols for the endpoints of the branches
c         G always represents Gamma
c         Rest from BSW (???) paper
c         per: W. Setyawan and S. Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c 
      CHARACTER*3  bvl
      CHARACTER*2  bsymb
      DIMENSION    bsymb(nbranchd)
c
      PARAMETER  ( nline = 128 )
      CHARACTER  readlbl*6,wantlbl*6,label*(nline)
      EQUIVALENCE( label, readlbl )
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a128)
 9100 format(1x,'# For definitions of symbols, use and cite (G=gamma):'
     $ /1x,'# W. Setyawan and S. Curtarolo, ',
     $       'Comp. Mat. Sci. 49, 299 (2010).' )
c
      read(IDAT,8000)  label
c
      wantlbl = 'bravai'
      if( readlbl.eq.wantlbl )then
        call LNWRIT( IWR, label, nline )
c       Read bravais lattice type
        read(IDAT,*)  bvl
        write(IWR,*)  bvl
c
        read(IDAT,8000)  label
        if( readlbl.eq. 'nbranc' )then
          call LNWRIT( IWR, label, nline )
          read(IDAT,*)  nbranch
          write(IWR,*)  nbranch,' = number of branches'
        else
          write(IWR,*) readlbl,' not >nbranch<  ... bye!'
          call STOPXERR( 'b.s. keyword >nbranch< not found' )
        endif
c
        read(IDAT,8000)  label
        if( readlbl.eq. 'symbol' )then
          call LNWRIT( IWR, label, nline )
          write(IWR,9100)
          do  i=1,nbranch
            read(IDAT,*)  bsymb(i)
c           Put in error step
            write(IWR,*)  bsymb(i)
          enddo
        else
          write(IWR,*)  readlbl,'  not >symbol<  bye!'
          call STOPXERR( 'b.s. keyword >symbol< not found' )
        endif
c
      elseif( readlbl.eq. 'gamma ' .or. readlbl.eq. 'Gamma ' )then
c       Just add gamma point to end of list
        call LNWRIT( IWR, label, nline )
        bvl = 'gam'
        nbranch = 1
        bsymb(1) = 'GG'
      else
        write(IWR,*)  readlbl,' not >bravais<  ... bye!'
        call STOPXERR( 'b.s. keyword >bravais< not found' )
      endif
c
      read(IDAT,8000)  label
      if(  readlbl.ne. 'end ba' )then
        write(IWR,*)  readlbl, ' not >end bands< ...  bye!'
        call STOPXERR( 'final b.s. keyword >end bands< not found' )
      else
        call LNWRIT( IWR, label, nline )
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
c End:
