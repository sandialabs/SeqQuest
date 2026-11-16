c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WKMEM
c
c
      subroutine PWKMEM( IWR, lvlout, maxwkd,
     $ do_gga,do_spin,do_blas3,idosolv,
     $ mat,nmat,nptr, ncplx,nk,nlat, norb,nstate, noad,natm, nhiste,
     $ mxang, ndbox,nbox,noboxs,
     $ near2c,irhoopt, idwf,  memwf,memwfw, iwf,iwfw,
     $ icoskr,isinkr, is1at,
     $ ioint,ioofk, ieigval,ieigpop,
     $ imbroye, membroye,
     $ i00, i01s, i01,i02,i03,i04,
     $ i05,i06,i07,i08,i09,i10,i11,i12 )
c---------------------------------------------------------------
c Purpose: to allot/assign memory for the whole code, with the
c          big workspace wk().  Check for memory sufficiency.
c          Make decision about on-the-fly v. stored grid orbs
c
c Written: Peter A. Schultz, 7-June-2001, for 2.47
c
c Revision history:
c   9Feb07-PAS/2.60: memory for in-core broyden; harmonize tp-serial;
c                    patch 24apr07, 4jun07
c  30Jun06-PAS/2.60: tighten memory usage in i03/i04/i05
c  11Sep05-PAS/2.59: static mem for linear dependence code (rpm)
c  18Aug05-PAS/2.59-tp21: reworking solver mem
c   5Nov04-PAS/2.59-tp: temporary hack to get needed mem to pschroed
c   3Oct02-PAS/2.54: reduce mem needed in setup
c  22Jan02-PAS/2.52: int-overflow max reset to (2**31-1)
c  25Oct01-PAS/2.50: int-overflow fix
c  20Jul01-PAS/2.49: spin-polarized dft (2*eigs, iwfw, purge fpop)
c  21JUn01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c **************************************************************
c *** THIS IS THE MOST DANGEROUS AND CRITICAL ROUTINE IN THE ***
c *** ENTIRE CODE. IT MANAGES AND ASSIGNS MEMORY EVERYWHERE. ***
c **************************************************************
c
c This routine partitions all big memory, everywhere in the code.
c Hence any changes in memory usage and requirements in Quest
c should be checked and validated against this routine, and
c appropriate updates made.
c
c We will be assigning "pointers" for all "irregular" N**2+ arrays
c that need to carry data through an scf cycle ...
c    ... needed in all major routines [mat0/grid0/scf/frc]:
c   icoskr  -> coskr(nk,nlat)
c   isinkr  -> sinkr(nk,nlat)
c    ... could be left until scf, but now needed in [mat0->scf]:
c   is1at   -> s1atom(noad,noad,nk,natm)
c    ... used only in scf, generated at end of set [set->scf]:
c   ioint   -> orbint(2,norb,nk)      [iff near2c=2, else skipped]
c   ioofk   -> orbofk(2,natm,norb,nk) [iff near2c=2, else skipped]
c    ... used only during scf [scf]:
c   ieigval -> eigval(nstate,nk)
c   ieigpop -> eigpop(nstate,nk)
c These should be staged carefully to be "allocated" only when
c necessary.
c
c And, the biggies, we assign "pointers" for the big workspaces
c that carry orbital matrices and grid fields:
c   i01s, i01,i02, ... ,i12
c The i00 denotes the first available (unused) space within the
c array, and should be simply "1", but just to be safe ...
c
c We will compute how much memory is required to keep grid bf, both
c fully in memory and on-the-fly, and flag whether there is enough
c memory to keep grid bf:
c in memory, or to compute them on-the-fly:
c   idwf = 0 ... compute grid orbitals on the fly
c        = 1 ... grid orbitals stored in memory
c
c We assign a pointer for the grid bf. The location of grid bf's
c depends on whether the bf's are saved (in which case this routine
c must make sure that other workspace does not step on them), or
c computed on-the-fly, in which case we can use any idle space:
c   iwf
c
c We assign space for grid workspace.  The amount of memory needed
c is dependent on the scheme chosen to evaluate the density
c (DM=density-matrix v. EF=eigenfcn), and whether we use blas3 or
c not. I only examine DM v. EF, and ignore blas3 (picking biggest):
c   iwfw
c
c Note that the grid bf and grid workspace is all single precision,
c while (for right now in 2.47) everything else is double.
c
c I need to know what the largest possible integer is in order to fend
c off integer overflow problems in some key parts of the following.
c I will not to do a full check of everything because it is too painful
c to code.  I will only do an overflow check for those operations that
c involve *optional* memory.  If required space goes over this value,
c the system is too big to do anyway, and I may as well let code bomb.
c 2**31 = 2 147 483 648
c For a 4-byte integer, (2**31-1) is the (safe) maximum:
      PARAMETER  ( maxint = 2 147 483 647 )
c
c input arrays:
c      ... flags that affect memory needs:
      LOGICAL    do_gga, do_spin, do_blas3
c      ... information about box breakup, and occupancies:
      DIMENSION  ndbox(2,2,3),noboxs(nbox)
c
c >>>> EXECUTABLE CODE:
c
      lstout = lvlout
c
      mem0 = i00 - 1
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Compute workspace needed for the "irregular" named spaces:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      ... needed in all major routines [mat0/grid0/scf/frc]:
      mcoskr  = nk*nlat
      msinkr  = nk*nlat
c      ... could be left until scf, but now needed in [mat0->scf]
      ms1at   = noad*noad*nk*natm
c      ... used only in scf, generated at end of set [set->scf]:
      moint   = 2*norb*nk
      moofk   = 2*natm*norb*nk
c      ... used only during scf [scf]:
      meigval = nstate*nk
      meigpop = nstate*nk
      if( do_spin )then
        meigval = 2*meigval
        meigpop = 2*meigpop
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Allocate space needed for "irregular" and regular spaces:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icoskr  = i00
      isinkr  = icoskr + mcoskr
      is1at   = isinkr + msinkr
c     coskr/sinkr/s1at is core allocation in setup
      i1set   = is1at + ms1at
      mem0set = i1set - 1
      inext   = is1at + ms1at
      if( near2c.eq.2 )then
c       Allocate space for two-center (orbint and orbofk) defect terms:
        ioint = inext
        ioofk = ioint + moint
        inext = ioofk + moofk
      else
c       Skip two-center contributions to atom density defects:
        ioint = is1at
        ioofk = is1at
      endif
      ieigval = inext
      ieigpop = ieigval + meigval
c Up to here is core memory in scf phase:
      i1scf   = ieigpop + meigpop
      mem0scf = i1scf - 1
c
c Bigger of Hamiltonian matrix, grid field (for first four slots):
      nwk = MAX( nmat, nptr )
c
      i01s = i1set
c
      i01 = i1scf
      i02 = i01 + nwk
      i03 = i02 + nwk
      i04 = i03 + nwk
      i05 = i04 + nwk
c GRIDRHO: i03i04 needs norb^2 (+nptr) or ncplx*norb*nstate before i05
c VSLOMAT: i03(i04) needs potentially 1.5*norb^2 before i05
c Hence, following "cmem" compaction does not work
cmem      n0304 = MAX( nwk+nptr, ncplx*norb*nstate )
cmemc      longer = nwk + nwk - n0304
cmemc      if( longer.lt.0 ) longer = 0
cmem      i05 = i03 + n0304
      i06 = i05 + nptr
      i07 = i06 + nptr
      i08 = i07 + nptr
      i09 = i08 + nptr
      i10 = i09 + nptr
      i11 = i10 + nptr
      i12 = i11 + nptr
      memixx = 4*nwk + 8*nptr
cmem      memixx = 2*nwk + n0304 + 8*nptr
      if( (mem0scf+memixx+1) .ne. (i12+nptr) ) call STOPXERR('wkmemi12')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Compute memory needed for grid basis functions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c NB: Grid orbitals are r*4, but memory in wk() is r*8.
c
c Survey boxes to get grid bf size: total, and biggest in single box
      nrwfbmx = 0
      nrwfall = 0
      nrmax = 0
      do  ibox=1,nbox
        nobox = noboxs(ibox)
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
        if( nrbox .gt. nrmax ) nrmax = nrbox
        nrwf = ncplx*nk*nobox*nrbox
        if( nrwf.gt.nrwfbmx ) nrwfbmx = nrwf
c
c with int overflow check:
        if( nrwfall .le. (maxint-nrwf) .and. nrwfall .ge. 0 )then
          nrwfall = nrwfall + nrwf
        else
          nrwfall = -1
        endif
      enddo
c
      if( nrwfall .gt. 0 )then
        if( lstout.gt.1 ) write(IWR,9010)  nrwfbmx,nrwfall
      else
        if( lstout.gt.1 ) write(IWR,9011)  nrwfbmx
      endif
 9010 format(/1x,'Need',i12,'-length r*4 array to pack grid wfcns/box',
     $       /1x,'Need',i12,'-length r*4 array to pack all grid wfcns')
 9011 format(/1x,'Need',i12,'-length r*4 array to pack grid wfcns/box',
     $ /1x,'Would need too much length array to pack all grid wfcns')
c
cwf4: fit r*4 orbs into r*8 space (all and biggest box):
      memwfa = ( nrwfall + 1 ) / 2
      memwfb = ( nrwfbmx + 1 ) / 2
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c            Compute workspace needed for grid routines
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c nc=ncplx,nk=nk,no=nobox,nsh=biggest(2l+1),nr=nrbox,np=npop=nstate
c            with blas3      without
c gridrho/dm  nc*no*nr     nc*nk*nsh*nr
c gridrho/ef  nc*np*nr     nc*nk*np *nr
c vslomat     nc*no*nr     nc*nk*nsh*nr
c vslofrc     nc*no*nr*nk  nc*nk*   *nr
c
      memdm = memwfb
      memef1 = memwfb
      memef2 = (ncplx*nk*nstate*nrmax+1)/2
      memef = MAX( memef1,memef2 )
      if( irhoopt .eq. 1 )then
        memwfw = memdm
      elseif( irhoopt .eq. 2 )then
        memwfw = memef
      else
        memwfw = MAX( memdm, memef )
      endif
      if( lstout.gt.1 ) write(IWR,9012)  2*memwfw
 9012 format( 1x,'Need',i12,'-length r*4 array for grid wfcns scratch')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Compute memory requirements in major code bottlenecks
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c >>>>> setup phase
c
c setup/"mat0set"
      memmat0 = MAX( 4*nmat,  (3*nmat+mxang) )
c setup/"grid0set"
      if( do_gga )then
        memgr0 = 14*nptr
      else
        memgr0 = 5*nptr
      endif
c setup/
      memset = mem0set + MAX( memmat0, memgr0 )
c
c >>>>> scf phase
c
c scf/eigensolvers
      if( ncplx .eq. 1 )then
c       Real eigensolver/erlsolv
CTP        memeig = 4*mat
ctp:start (rpm)
c       How many processors on this image?
        call MPNODES( nprocs )
        if( nprocs .eq. 1 )then
c         Serial solver, just needs ...
          memeig = 4*mat
        else
c         Parallel solver, needs extra for now ...
          memeig = 5*mat
        endif
ctp:end
      elseif( ncplx .eq. 2 )then
c       Complex eigensolver/eigsolv
        memeig = mat + norb*(4*norb+2*nstate) + 7*norb
        if( idosolv.eq.3 )then
c         Linear dependence treatment in SCHROED needs more
          memeig = mat + norb*(6*norb) + 7*norb
        endif
ctp:start
c       Allot extra memory for complex parallel solver pschroed:
        call MPNODES( nprocs )
        if( nk.eq.1 .and. nprocs.gt.1 ) memeig = memeig + 3*norb*norb
ctp:end
      else
        call STOPXERR( 'wkmem-nc/ncplx is disallowed value' )
      endif
c scf/doengy,blend
      if( do_spin )then
        memengy = 3*nmat + 2*mat
        nvecbl = 2*nmat
      else
        memengy = nmat + 2*mat
        nvecbl = nmat
      endif
      call MPNODES( nprocs )
      memblnd = 2*nvecbl + 2*nvecbl
c scf/blendhistory
c     Determine how much memory needed to in-core blend history
      call MPNODES( nprocs )
      nveclocl = ( nvecbl + nprocs - 1 ) / nprocs
      memblnd = 2*nvecbl + 2*nveclocl
      memhiste = 2*nveclocl + 2*nhiste*nveclocl
c scf/*before wf*
      memscf = mem0scf + MAX( memixx, memeig, memengy, memblnd )
c
c >> grid orbitals
c
c Allocate space for grid wf temp work space(location is critical!!!):
c Must see memory usage around gridrho/vslomat/vslofrc for this to work.
      iwfw = i07
c Compute space available for grid wf ...
c  ... idiot-check that temp wf space does not spill over end of ixx,
c      to ensure that stored grid orbitals are safe
      memscfw = iwfw + memwfw
      if( memscf .gt. memscfw ) memscfw = memscf
c  ... memory available (and used) for in-core grid wf:
      iwfall = memscfw + 1
      maxwfa = maxwkd - memscfw
      memall = memscfw + memwfa
c  ... and for on-the-fly grid wf:
      iwffly = iwfw + memwfw
      memfly = iwffly + memwfb - 1
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Now we put it all together and make decisions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      minmem = MAX( memset, memscf, memfly )
c
      write(IWR,9020)  minmem
 9020 format(1x,'Minimum maxwkd for scf to complete=',i12)
c     Note: vslofrc uses extra memory that may exceed this minimum
c
      if( minmem .gt. maxwkd )then
        write(IWR,*) '>>>>> need maxwkd=',minmem
        write(IWR,*) '>>>>> have maxwkd=',maxwkd
        call STOPXERR( 'maxwkd  - insufficient memory for scf calc' )
      endif
c
      if( nrwfall .gt. 0 .and. memall .le. maxwkd )then
c       All wave fcns is not too many, and all wvfcns fit
        write(IWR,9021)  memwfa, maxwfa
 9021   format(1x,'All grid wfcns stored in',i10,' of',i10,' available')
        idwf = 1
        iwf = iwfall
        memwf = memwfa
        memtopscf = memall
      elseif( memfly .le. maxwkd )then
c       Have to do on-the-fly grid orbitals
        write(IWR,9022)  memwfb, maxwfa
 9022   format(1x,'Box-wise grid wfcns use ',i10,';  ',i10,' available')
        if( maxwfa .lt. memwfb ) write(IWR,*) ' for in-core grid wfcns'
        idwf = 0
        iwf = iwffly
        memwf = memwfb
c       Top mem usage is either on-the-fly scratch, or scf sandbox:
        memtopscf = MAX( memfly, memscf )
      else
        call STOPXERR( 'wkmem   - logic failure in WKMEM routine' )
      endif
c
c See if we can squeeze in-core broyden into remaining space::
c
      memfree = maxwkd - memtopscf
      if( memhiste .le. memfree )then
c       Have room for full in-core broyden
        imbroye = memtopscf + 1
        membroye = memhiste
      else
c       Not room for full in-core broyden - set it to zero
        imbroye = 1
        membroye = 0
      endif
c
c    That's all Folks!
c
      RETURN
      END
