c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POPCLOSE
c
c
      subroutine POPCLOSE( icloseocc,dokdefct,ikdefct,
     $ efermi,egap, nstbulk,
     $ etemp,edegen, elnumbr,esumsp, occfull,
     $ nstate,nk, wtk, eigval,eigpop,npop,
     $ wlvl )
c---------------------------------------------------------------
c Purpose: set level populations using closed shell, and defect
c          sampling (gamma pt defect level)
c
c Written: Peter A. Schultz, 16-February-2005, for v2.59 code
c
c Revision history:
c  17Oct08-PAS/2.62: fix esumsp, installing proper wtk()
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.d0 )
c
c Input:
      LOGICAL    dokdefct
      DIMENSION  wtk(nk), eigval(nstate,nk)
c Output:
      DIMENSION  eigpop(nstate,nk)
c Scratch:
      DIMENSION  wlvl(*)
c
c Local:
      LOGICAL    closeshl
c
c >>>> EXECUTABLE CODE:
c
c Clear out occupancy levels:
      do  ik=1,nk
        call MKZERO( nstate, eigpop(1,ik) )
      enddo
c
c If we have no electrons, skip everything
      if( elnumbr .le. zero )then
        npop = 0
        esumsp = 0
        gap = -one
        efermi = eigval(1,1)
        goto 900
      endif
c
c Now go over each k-point and populate with
c a fixed number of electrons in each:
c
c First, check if doing defect calculation, in
c which case will populate only bulk (non-defect)
c states in the following:
      elbulk = elnumbr
      if( dokdefct ) elbulk = occfull*nstbulk
c
c Second, set closed shell flag if to do 0K occupations
      closeshl = .false.
      if( icloseocc .eq. 2 ) closeshl = .true.
      npop = 0
      do 1000 ik=1,nk
c       Use wt of one - to get closed-shell e-count
        eltot1k = elbulk
        wtocc = one
c
        call POP1K( closeshl,
     $   fermilvl,gap1k, ehomo1k,elumo1k,
     $   etemp,edegen, eltot1k,occfull,
     $   nstate,npop1k, wtocc,
     $   eigval(1,ik),eigpop(1,ik), wlvl )
c
        npop = MAX( npop, npop1k )
        if( ik.eq.1 )then
          efermi = fermilvl
          ehomo = ehomo1k
          elumo = elumo1k
        else
          efermi = MAX( efermi, fermilvl )
          ehomo = MAX( ehomo, ehomo1k )
          elumo = MIN( elumo, elumo1k )
        endif
 1000 continue
      egap = elumo - ehomo
c
c If we are doing a defect-sampling calculation,
c need to go back and put the defect electrons into
c the defect states at the defect k-point (presumably
c the gamma point):
      if( dokdefct )then
        if( npop .ne. nstbulk )then
cxxx:start
cxxx      New feature - dump added diagnostic if fails.
          call FLGETIWR( IWR )
          write(IWR,*) 'ERROR/POPCL: npop, nstbulk=',npop,nstbulk
          write(IWR,*) 'POPCL: levels(-3/+5) and occups near Ef:'
  913     format( 1x,i4,5(1x,f10.4,f7.3) )
          do  i=nstbulk-3,nstbulk+5
            write(IWR,913)  i, (eigval(i,ik),eigpop(i,ik), ik=1,nk)
          enddo
cxxx:end
          call STOPXERR( 'defect state is in the VB' )
        endif
c
        ikdef = IABS( ikdefct )
        if( ikdefct .lt. 0 )then
c         k-point outside of special point sample, blank bulk occs
          do  ilvl=1,nstbulk
            eigpop(ilvl,ikdef) = zero
          enddo
        endif
c
c       Occupy the defect state(s):
c        ... nstumo=number of available states above bulk states
        nstumo = nstate - nstbulk
c        ... eldefct=number of defect electrons
        eldefct = elnumbr - occfull*DBLE( nstbulk )
c        ... rescale "full" to get total electron count right
        if(  wtk(ikdef).le.zero ) call STOPXERR( 'wt(kdefct)=0' )
        eldef = eldefct / wtk(ikdef)
        fulldef = occfull / wtk(ikdef)
        wtdef = one
        ilvldef = nstbulk + 1
        call POP1K( closeshl,
     $   fermidef,gapdef, ehomodef,elumodef,
     $   etemp,edegen, eldef  ,fulldef,
     $   nstumo,npopdef, wtdef,
     $   eigval(ilvldef,ikdef),eigpop(ilvldef,ikdef), wlvl )
c
        npop = npop + npopdef
        efermi = fermidef
        egap = elumodef - ehomo
      endif
c
c Sum occupation-weighted eigenvalues to get "band energy" sum:
c
      esumsp = zero
      do  ik=1,nk
        esum1k = zero
        wtki = wtk(ik)
        do  ilvl=1,npop
          esum1k= esum1k + wtki*eigpop(ilvl,ik)*eigval(ilvl,ik)
        enddo
        esumsp = esumsp + esum1k
      enddo
c
  900 continue
c
      RETURN
c
c    That's all Folks!
c
      END
