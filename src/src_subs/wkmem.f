c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WKMEM
c
c
      subroutine WKMEM( IWR, lvlout, maxwkd,
     $ do_gga,do_spin,do_blas3, Lijkmat,
     $ idosolv, do_kppsolve,do_psolv, do_eigvecs,Lnonscf,
     $ mat,nmat,nmatkp,  nptr, nhiste,
     $ ncplx,nk,nlat, norb,nstate, noad,natm,ntyp,itypa,norba,
     $ mxang, ndbox,nbox,noboxs,mxrb,
     $ near2c,irhoopt, idwf,  memwf,memwfw, iwf,iwfw,
     $ icoskr,isinkr, is1at,
     $ ioint,ioofk, ieigval,ieigpop,
     $ idbr, imblnde, memblnde,
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
c  14Dec14-PAS/PAS/2.66: adapt (minor!) for ASD e-blend
c  17Jan13-PAS+ACP/2.64: k-parallelparallel memory reductions, and flags
c  26Feb13-PAS/2.63b: bugfix - parallel stripe v. square(ijk) consistency
c   9Dec08-PAS/2.62: reduced mem p-schroed
c  15Oct08-PAS/2.62: stripe-parallel mat0set/frc
c  22Sep09-PAS/2.62: split dmat spin file; reduce doengy mem; err retrace
c  23Sep08-PAS/2.61h: error in m4 memory arithmetic
c  24May08-PAS/2.61d: reduce blend-assigned mem to get large spin-calcs
c  17May08-PAS/2.61c: fix the vslomat memory allocation
c  12Feb08-PAS/2.61: pack i03/i04 tighter through vslomat/gridrho
c   4Jan08-PAS/2.61: reduced (dist) memory for grid wvfcns
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
c (DM=density-matrix v. EF=eigenfcn), and whether we use blas3 or not.
c I only examine DM v. EF, and ignore blas3 flag (picking biggest):
c   iwfw
c
c Note that the grid bf and grid workspace is all single precision,
c while (for right now in 2.47) everything else is double precision.
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
c Input arrays:
c      ... flags that affect memory needs:
      LOGICAL    do_gga, do_spin, do_blas3
      LOGICAL    Lijkmat
      LOGICAL    do_kppsolve, do_psolv
c
      LOGICAL    do_eigvecs, Lnonscf
c       do_eigvecs = flag determines if eigenvectors (and memory) are needed
c       Lnonscf    = flag to turn off contributions to memory not-used by non-SCF
c                    parts of the code. (e.g. Broyden)
c
c      ... information about atom types, and norb/atom:
      DIMENSION  itypa(natm),norba(ntyp)
c      ... information about box breakup, and occupancies:
      DIMENSION  ndbox(2,2,3),noboxs(nbox)
c
c Local declarations:
      DATA  lvldbg / 0 /
      LOGICAL  flagdum
      LOGICAL  USE_KPSTORE
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM  / .false. /
c
c >>>> EXECUTABLE CODE:
c
      USE_KPSTORE = .true.
c     For nk=1, old scheme is smaller by nmatpr:
      if( nk .eq. 1 ) USE_KPSTORE = .false.
cxxx: declare dev-version of wkmem:
      write(IWR,*) 'DEV264h/wkmem: perlsolv dist-mem'
      write(IWR,*) 'DEV/wkmem: KPSTORE=',USE_KPSTORE
      write(IWR,*) 'DEV/wkmem: psolv version, input do_psolv=',do_psolv
      lstout = MAX( lvlout, lvldbg )
c
  100 continue
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
      nspin = 1
      if( do_spin )then
        nspin = 2
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
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Determine allocation sizes needed for stripe-parallel setup routines:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
c     
c Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
      nmatstr = ncplx*nk*matstr
      if( lstout .gt. 2 )then
        write(IWR,*) 'WKMEM: full striped matrix size, all-k storage=',
     $                       matstr,nmatstr
      endif
c Need maximum value of norba for memory allocation.
      maxnorba = norba(1)
      do  ityp=2,ntyp
        maxnorba = MAX( maxnorba, norba(ityp) )
      enddo
c           
c Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
      nmatpr = mymatpr*ncplx*nk
      mxnmatpr = maxmatpr*ncplx*nk
      if( nprocs.gt.1 .and. lstout.gt.2 )then
        write(IWR,*)'WKMEM: striped-mat; nprocs,active=',nprocs,nprocstr
        write(IWR,*)'WKMEM: striped/proc: my, max=',nmatpr,mxnmatpr
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Determine allocation sizes for other parts of the code:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Bigger of Hamiltonian matrix, grid field (for first couple slots):
      nwk = MAX( nmatkp, nptr )
c
c Find the largest number of orbitals in a box
      noboxmax = 1
      do  ibox=1,nbox
        nobox = noboxs(ibox)
        if( nobox .gt. noboxmax ) noboxmax = nobox
      enddo
c Get local MP info (here, the k-parallel partitions)
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
c
      call KPMINE( nk, nkloc, nk0 )
c
c Biggest headache is budgeting around vslomat (gridrho less so)
c How much scratch space do we need (beginning at i03, i05 reserved)?
c noboxmax is max number of orbitals overlapping a given box
c
c vslomat requirements
c  vijk:     nspin*nk*norb^2 for Hams, in i01
c  scratch: in i01+nspin*nmat (all these r4)
c   without blas3:   nk*norb^2
c   with blas3   :   nk*norb^2 + nc*nk*noboxmax^2 (i.e. add last term)
c   with parallel:   nk*norb^2 + nk*norb^2 (i.e. add last term)
c                    (last term is target in matrix reduce)
c gridrho requirements
c  dmatrd: norb^2[r8], buffer to read in dmat from disk, k-by-k,
c    and transfer in into [r4] dmat(nc*nk*norb^2)=nmat space.
c  DM-blas3: nc*nk*norb^2 + nc*nk*noboxmax^2
c  DM-w/o  : nc*nk*norb^2
c  EF-blas3: nc*nk*norb*nstate + nc*nk*noboxmax*nstate
c  EF-w/o  : nc*nk*norb*nstate
c
c For vslomat, compute scratch spave needed, counting from i03,
c     Beginning with space devoted to running ham vij [r8] storage
      mvslomat = -nwk
      if( do_spin ) mvslomat = mvslomat + nwk
c     Now work out amount of r4 workspace needed
c
c     This to compute r4 dense ham and r4 full ham (and r4 reduce also)
      m4vsl1 = nmatkp
cxxxcpas: acp had converted "ncplx" to "2" also ... is this needed? Apparantly, no.
      if( do_blas3 ) m4vsl1 = m4vsl1 + ncplx*nkloc*noboxmax**2
c
      m4vsl2 = nmatkp
c     If nprocs .gt. 1, need another nmatkp as buffer for MPI_Reduce
C      call MPNODES( nprocs ) ! switch to k-parallel
      call MPNODES_K( nprocs )
      if( nprocs .gt. 1 ) m4vsl2 = nmatkp + nmatkp
c
      m4vsl = MAX( m4vsl1, m4vsl2 )
      mvslomat = mvslomat + (m4vsl+1)/2
c
c For gridrho, use larger blas3-active requirement, too.
c  However, know that vslomat memory requirements trump this.
c  Lots of space to play with here, dmat into i01, dmatrd(mat) into i02.
c  And then enough room for everything in i02 after dmat reads done.
      mgridrho = 0
c  Compute memory needs offset from i03, assuming i02(nwk)
      mgridrho = -nwk
c  dmatrd space for rhdmat/rheigvec
c  Assume i01 (dmat), need r8(mat)=r4(2mat) in i02 buffer to read dmatrd:
      m4grh1 = 2*mat
c  dmatrd space for grhodm*/grhoef*
c  And then r4 memory needed as scratch after dmat obtained
      m4grh2 = 0
C      if( do_blas3 ) m4gr2 = ncplx*nkloc*noboxmax**2
      if( do_blas3 )then
        if( irhoopt.eq.1 )then
c         irhoopt.eq.1 (DM-density matrix)
          m4grh2 = ncplx*nkloc*noboxmax**2
        else
c         irhoopt.eq.2 (EF-eigenfunction)
          m4grh2 = ncplx*nkloc*noboxmax*nstate
        endif
c         Alternately since noboxmax<=norb and nstate<=norb, universal:
c         m4grh2 = ncplx*nkloc*noboxmax*norb
      endif
c dmatrd space
      m4grh = MAX( m4grh1, m4grh2 )
c Convert to r8 units
      mgridrho = mgridrho + (m4grh+1)/2
c
c For grdovlp, it can poke above ceiling, too, in k-parallel (Code patch by ACP)
c   NB: currently, arrays are declared and used as nk rather than kparallel nkloc
c   NB: this is cause to push k-parallel into various overlap arrays used by grdovlp <sigh>
      meikr = 2*nk*mxrb
c grdovlp requirement relative to i03 (hence -2nwk)
      mgrdovlp = meikr + MAX( 2*norb*nk, noad*noad*nk*natm ) - 2*nwk
c
      m0304 = MAX( 2*nptr, mvslomat, mgridrho, mgrdovlp )
      if( lstout .gt. 2 )
     $ write(IWR,9130) 'WKMEM: max orbs in box noboxmax=',noboxmax
      if( lstout .gt. 3 )
     $ write(IWR,9130) 'DEV/WKMEM: m0304,2nwk,mgrdo,mvsl,2nptr=',
     $                     m0304,2*nwk,mgrdovlp,mvslomat,2*nptr
c
      i01s = i1set
c
      i01 = i1scf
      i02 = i01 + nwk
      i03 = i02 + nwk
C      i04 = i03 + nwk
C      i05 = i04 + nwk
      i04 = i03 + nwk
      i05 = i03 + m0304
      if( (i04+nptr) .gt. i05 ) i04 = i05 - nptr
      i06 = i05 + nptr
      i07 = i06 + nptr
      i08 = i07 + nptr
      i09 = i08 + nptr
      i10 = i09 + nptr
      i11 = i10 + nptr
      i12 = i11 + nptr
cmem      memixx = 4*nwk + 8*nptr
      memixx = 2*nwk + m0304 + 8*nptr
c Idiot check that total-mem matches incrementally added mem amount
      if( (mem0scf+memixx+1) .ne. (i12+nptr) ) call STOPXERR('wkmemi12')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Compute memory needed for grid basis functions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c NB: Grid orbitals are stored r*4, but space in wk() is declared r*8.
c     Hence, we have 2X declared space to store wvfcns
c
c Survey boxes to get grid bf size: total, and biggest in single box
      nrwfbmx = 0
      nrwfall = 0
      nrmax = 0
      do 400  ibox=1,nbox
c
c       Each processor will do every nprocs'th box in this loop
        if( MOD( ibox, nprocs_k ) .ne. iproc_k ) goto 400
c
        nobox = noboxs(ibox)
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
        if( nrbox .gt. nrmax ) nrmax = nrbox
cdev        nrwf = ncplx*nk*nobox*nrbox
        nrwf = ncplx*nkloc*nobox*nrbox
        if( nrwf .gt. nrwfbmx ) nrwfbmx = nrwf
c
c       Add to full mem requirement, with int overflow check:
        if( nrwfall .le. (maxint-nrwf) .and. nrwfall .ge. 0 )then
          nrwfall = nrwfall + nrwf
        else
          nrwfall = -1
        endif
c
  400 continue
c
      if( nrwfall .gt. 0 )then
        if( lstout .gt. 1 ) write(IWR,9040)  nrwfbmx,nrwfall
      else
        if( lstout .gt. 1 ) write(IWR,9041)  nrwfbmx
      endif
 9040 format(/1x,'Need',i12,'-length r*4 array to pack grid wfcns/box',
     $       /1x,'Need',i12,'-length r*4 array to pack all grid wfcns')
 9041 format(/1x,'Need',i12,'-length r*4 array to pack grid wfcns/box',
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
c nc=ncplx,nk=nkloc,no=nobox,nsh=biggest(2l+1),nr=nrbox,np=npop=nstate
c            with blas3      without
c gridrho/dm  nc*no*nr     nc*nk*nsh*nr
c gridrho/ef  nc*np*nr     nc*nk*np *nr
c vslomat     nc*no*nr     nc*nk*nsh*nr
c vslofrc     nc*no*nr*nk  nc*nk*   *nr
c
      memdm = memwfb
      memef1 = memwfb
cdev      memef2 = (ncplx*nk*nstate*nrmax+1)/2
      memef2 = (ncplx*nkloc*nstate*nrmax+1)/2
      memef = MAX( memef1,memef2 )
      if( irhoopt .eq. 1 )then
        memwfw = memdm
      elseif( irhoopt .eq. 2 )then
        memwfw = memef
      else
        memwfw = MAX( memdm, memef )
      endif
      if( lstout .gt. 1 ) write(IWR,9042)  2*memwfw
 9042 format( 1x,'Need',i12,'-length r*4 array for grid wfcns scratch')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Compute memory requirements in major code bottlenecks
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c >>>>> setup phase
c
c Default to the (usually much) smaller distributed stripe mat0set
      Lijkmat = .false.
c
c setup/"mat0set"/ijk-tp
      memijkv = 3*nmat + mat
      memijka = 2*nmat + mxang
      if( nprocs .gt. 1 )then
c       Count up tp-reduce space as well:
        memijkv = 4*nmat
        memijka = 2*nmat + mxang
      endif
      memijk  = MAX( memijkv, memijka )
c
c setup/"mat0set"/stripes ...
c   ... Merge all-k striped matrices
c NB: use mxnmatpr, NOT local nmatpr, to ensure consistency across procs
      memstr1v = 2*mxnmatpr + nmatstr + mat
C      memseta = 2*nmatpr + nmatstr + MAX( mxang, mat )
      memstr1a = 2*mxnmatpr + mxang
      memstr1 = MAX( memstr1v, memstr1a )
c   ... Merge k-by-k striped matrices reduced memory footprint store
      memstr2v = 3*mxnmatpr + MAX( matstr*ncplx, 3*norb*maxnorba ) + mat
      memstr2a = 2*mxnmatpr + mxang
      memstr2 = MAX( memstr2v, memstr2a )
c
      if( DEBUG_MEM .or. lvlout .gt. 3 )then
      write(IWR,*) 'DEV/WKMEM: memstr2v,memstr2a=',memstr2v,memstr2a
      write(IWR,*) 'DEV/WKMEM: nmatpr,matstr,mat=',nmatpr,matstr,mat
      write(IWR,*) 'DEV/WKMEM: norb,maxnorba,mxang=',norb,maxnorba,mxang
      endif
c
c     Record stripe memory requirements appropriate to the branch chosen
      if( .not. USE_KPSTORE )then
        memstr = memstr1
      else
        memstr = memstr2
      endif
c
      if( DEBUG_MEM .or. lstout .gt. 3 )then
        write(IWR,*) 'MEM/WKMEM: memstr(old,new)=',memstr1,memstr2
        write(IWR,*) 'MEM/WKMEM: USE_KPSTORE=',USE_KPSTORE,memstr
      endif
c
c Decide whether to revert back to ijkmat for setup:
      msetstr = mem0set + memstr
      msetijk = mem0set + memijk
c
      if( DEBUG_MEM .or. lstout .gt. 3 )
     $ write(IWR,*) 'DEV/MEM: msetstr,msetijk=',msetstr,msetijk
      if( msetstr.gt.maxwkd .and. msetijk.le.maxwkd )then
c       The stripe-parallel does not fit, but the ijk-mat does
        Lijkmat = .true.
      endif
      if( ncplx.eq.2 .and. natm.le.4 .and. msetijk.lt.msetstr )then
c       With very small systems, the ijk-mat will be better
        Lijkmat = .true.
      endif
cxxx264c MAT0SET1 is broken (PAS: really?, I've not seen this)
cxxx264      Lijkmat = .false.
cxxx263b: Should have reduce bcast here to ensure procs have same info
c
c
      if( DEBUG_MEM .or. lstout .gt. 3 )then
        write(IWR,*) 'DEV/MEM: nmatpr,nmatstr=',nmatpr,nmatstr
        write(IWR,*) 'DEV/MEM: nmat,mat,mxang=',nmat,mat,mxang
        memold = MAX( 4*nmat,  (3*nmat+mxang) )
        write(iWR,*) 'DEV/MEM: memstr,ijk,old=',memstr,memijk,memold
      endif
c
      memmat0 = memstr
      if( Lijkmat )then
c       We will use the ijkmat branch in the setup
        memmat0 = memijk
        if( lstout .gt. 0 )
     $   write(IWR,*) 'WKMEM: using ijk+Reduce setup integrals'
      else
        if( lstout .gt. 0 )
     $   write(IWR,*) 'WKMEM: using stripe-parallel+Merge  setup'
      endif
c
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
c
c       How many processors on this image?
        call MPNODES( nprocs )
        if( nprocs .eq. 1 )then
c         Serial lapack solver, just needs ...
          memeig = 4*mat
        else
cxxx264   Memory needs for distributed parallel real solver needs doing
c          memeig = 5*mat
          call EIGPROCS( nprocs, nprow,npcol, npsolver )
          call EIGMDIST( matdist, norb, npsolver )
          memeig = 4*matdist + mat
        endif
c
      elseif( ncplx .eq. 2 )then
c       Complex eigensolver/eigsolv
c
c       Serial lapack complex eigensolver memory needs (as baseline):
        memeig = mat + norb*(4*norb+2*nstate) + 7*norb
        if( idosolv.eq.3 )then
c         Linear dependence treatment in SCHROED needs more
          memeig = mat + norb*(6*norb) + 7*norb
        endif
c
c       Allot distributed memory for complex parallel solver pschroed:
c       ***** MUST be compatible with branching in PEIGSOLV *****
c
        call MPNODES( nprocs )
c
c       Determine if we MUST avoid serial-solves due to memory limits
        memscfmx = maxwkd - i01 + 1
        memeigmx = 2*(memscfmx/3)
        if( nprocs.ge.4 .and. memeig.gt.memeigmx )then
c         We want some slack above memeig for the blend, etc.
c         In future design, we might use an "npsolv" rather than do_psolv
          do_psolv = .true.
          if( DEBUG_MEM .or. lstout .gt. 0 )then
          write(IWR,*) 'DEV/wkmem: turned psolv to do_psolv=',do_psolv
          write(IWR,*) '>>>>> mem needed for seq solves=',memeig
          write(IWR,*) '>>>>> limited to ',memeigmx,' of free ',memscfmx
          endif
        endif
c
cxxx264: the solver (memory) branches and their criteria need to be cleaned up
        if( do_kppsolve
C  The following is now determined in main (outside wkmem)
C     $     .and. ( norb.ge.200*(nprocs/nk+1) )
     $                  )then
c         New k-parallel-parallel solver, invokes ScaLAPACK solvers within k-point groups
c         Sizing routine in pschroedkp.f
          if( lstout .gt. 0 )
     $    write(IWR,*) 'WKMEM: Use k-parallel-parallel solver'
          flagdum = .false.
ctp
CTPc kppmemx call is mpi-/blacs-dependent, needs be removed for serial
CTPc kppmem0 is mpi-/blacs-independent, and is good enough for parallel, too
CTP         if( DEBUG_MEM )then
CTP            call KPPMEMX( flagdum, do_eigvecs,
CTP     $                    norb,nstate, memeig )
CTP            write(IWR,*) 'DEV/WKMEM: kppmemx memeig=',memeig
CTP          endif
          call KPPMEM0( flagdum, do_eigvecs,
     $                  norb,nstate, memeig )
          if( DEBUG_MEM .or. lstout .gt. 2 )
     $    write(IWR,*) 'DEV/WKMEM: kppmem0 memeig=',memeig
c
        elseif( nprocs.eq.1 .or.
     $         (  .not. do_psolv .and.
     $           ( norb.lt.800 .or. nk.gt.4 .or. nprocs.lt.(4*nk) )  )
     $                           )then
c         Use serial lapack solver
c          IF 1 proc (already serial)
c          OR if not compelled to be parallel-solves AND:
c              a small problem (little benefit from parallel)
c           OR more than 4 k-points (do k-distrbuted serial solves)
c           OR too few procs for a parallel-parallel ...
c               (this last is likely moot with the new k-parallel-parallel)
          if( DEBUG_MEM .or. lstout .gt. 0 )
     $    write(IWR,*) 'WKMEM: Use serial complex solver in PSCHROED2'
          if( DEBUG_MEM .or. lstout .gt. 2 )
     $    write(IWR,*) 'DEV/WKMEM: serial memeig=', memeig
c
        else
c         Parallel-solver pschroed
c         Get number of active procs
          call EIGPROCS( nprocs, nprow,npcol, npsolver )
c         Get size of distributed matrix
C          call EIGMDIST( matdist, norb, npsolver )
          call PSCHWKSZ0(            npsolver, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
c         Memory bottlenecks in pschroed have following requirements:
          nwpreph = 2*matdist + 3*mat
          nwpreps = 4*matdist + 2*mat
          nwprep = MAX( nwpreph, nwpreps )
          nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
C          nwclose = 6*matdist + 2*norb*nstate ! this is more precise?
          nwclose = 6*matdist + 2*norb*norb
          nwpsch = MAX( nwprep, nwsolv, nwclose )
          memeig = nwpsch
c
          if( DEBUG_MEM .or. lstout .gt. 0 )
     $    write(IWR,*) 'WKMEM: Use parallel complex solver in PSCHROED'
          if( DEBUG_MEM .or. lstout .gt. 2 )
     $    write(IWR,*) 'DEV/WKMEM: parallel memeig=',memeig
        endif
c
      else
        call STOPXERR( 'wkmem-nc/ncplx is disallowed value' )
      endif
c
c scf/doengy
c  Need to track all nspin hams, and also dmat and h0mat for single k
      memengy = nspin*nmatkp + 2*mat
c
c scf/blend&blendhistory
c  Length of full k-parallel blend vector
      nvecbl = nspin*nmatkp
c
c     Determine size of data-distributed blend Hams
      call MPNODES_K( nprocbl )
      nveclocl = ( nvecbl + nprocbl - 1 ) / nprocbl
c
c     Now minimum memory we need to run blend ...
c      ... upon input (and in gesser) use two full sets of Hams:
      memblndi = 2*nvecbl
c      ... and, if k-parallel, further need read buffer(mat) in h0readkp:
      if( nkloc .lt. nk ) memblndi = memblndi + mat
c      ... in blend, pass a full set plus three local scratch spaces:
      memblnds = nvecbl + 3*nveclocl
c      ... and the workspace we need is the larger of these two:
      memblnd = MAX( memblndi, memblnds )
      if( Lnonscf ) memblnd = 0
c
c     Determine how much memory would need to keep blend history in-core:
c     (Note: ASD has nhiste=0 need)
      memhiste = 2*nveclocl + 2*nhiste*nveclocl
      if( Lnonscf ) memhiste = 0
c
c scf/*before wf*
c     memscf is largest extent of scratch memory used in scf.
c     Beyond memscf can be used for storage, wil survive from on iter into next
      memscf = mem0scf + MAX( memixx, memeig, memengy, memblnd )
c
      if( DEBUG_MEM .or. lstout .gt. 3 )then
        write(IWR,9130) 'WKMEM: mem0scf=',mem0scf
        write(IWR,9130) 'WKMEM: memixx/eig/engy/blnd=',
     $   memixx, memeig, memengy, memblnd
      endif
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
      if( lstout .gt. 0 ) write(IWR,9020)  minmem
 9020 format(1x,'Minimum maxwkd for scf to complete=',i12)
c     Note: vslofrc uses extra memory that may exceed this minimum
c
      if( lstout .gt. 2 )then
        write(IWR,9130)  'WKMEM: from memset,memscf,memfly=',
     $   memset,memscf,memfly
      endif
 9130 format(1x, a, 6(1x,i12) )
c
      if( minmem .gt. maxwkd )then
        write(IWR,*) '>>>>> need maxwkd=',minmem
        write(IWR,*) '>>>>> have maxwkd=',maxwkd
        if( lstout .lt. 5 )then
c         Go back and give more detail ...
          write(IWR,*) ' retracing to get more detail on memory needs'
          lstout = 5
          goto 100
        endif
        call STOPXERR( 'maxwkd  - insufficient memory for scf calc' )
      endif
c
      if( nrwfall .gt. 0 .and. memall .le. maxwkd )then
c       All wave fcns is not too many, and all wvfcns fit
        write(IWR,9021) 'All grid wfcns stored in', memwfa,
     $                  ' of', maxwfa,' available'
 9021   format(1x, a,i10,a,i10,a )
        idwf = 1
        iwf = iwfall
        memwf = memwfa
        memtopscf = memall
      elseif( memfly .le. maxwkd )then
c       Have to do on-the-fly grid orbitals
        write(IWR,9021)  'Box-wise grid wfcns use ',memwfb,';  ',
     $                    maxwfa,' available'
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
c See if we can squeeze in-core broyden into remaining space:
c
      memfree = maxwkd - memtopscf
      if( memhiste .le. memfree )then
c       Have room for full in-core broyden - set memblnde location, size
        write(IWR,9021) 'WKMEM: in-core Broyden! uses',memhiste,
     $                  ' out of', memfree,' available'
        imblnde = memtopscf + 1
        memblnde = memhiste
        idbr = 1
      else
c       Not room for full in-core broyden - set memblnde to 0
        imblnde = 1
        write(IWR,9021) 'WKMEM: on-disk Broyden! for in-core',memhiste,
     $                  ' needed, but only', memfree,' available'
        memblnde = 0
        idbr = 0
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MEMCHK
c
c
c###############################################################
c  Purpose: check for integrity of in-core stored arrays,
c           e.g., grid orbitals or blending history,
c           to check if memory is stepped on by memory overrun in job
c           This is solely debugging, a user should never see this
c
c Written: Peter A. Schultz, 11-Jabuary-2012, for 2.63
c          based on a suggestion from Andy Pineda
c
c Revisions:
c  none
c###############################################################
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MEMCHKIWF
c
c
      subroutine MEMCHKIWF( idwf, itypchk, labelwf, wvfcns, wfiwf )
c---------------------------------------------------------------
c Purpose: check the grid orbitals
c  idwf = flag if grid orbitals are indeed stored
c---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      REAL    wvfcns, wfiwf, wfiwf1
      DIMENSION  wvfcns(*)
      CHARACTER*(*) labelwf
c
c Local declarations:
      DATA  lvldbg / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( idwf .le. 0 ) goto 999
c
      call FLGETIWR( IWR )
c
 1800 format(a,a,1p4d10.3)
      if( itypchk .eq. 0 )then
c       Initialize
        wfiwf = wvfcns(1)
        if( lvldbg .gt. 0 )
     $  write(IWR,1800) 'IWFINIT/wf(1): ',labelwf, wfiwf
      elseif( itypchk .eq. 1 )then
c       Check against recorded value
        wfiwf1 = wvfcns(1)
        if( lvldbg .gt. 0 )
     $  write(IWR,1800) 'IWFCHK/wf(1): ',labelwf, wfiwf1
        if( wfiwf1 .ne. wfiwf )then
          write(IWR,1800) '***** IWFBUG *****',labelwf,wfiwf,wfiwf1
          call STOPXERR( 'Memory error - grid wave fcns overwritten' )
        endif
      elseif( itypchk .eq. 2 )then
c       Just output the value without check
        wfiwf1 = wvfcns(1)
        write(IWR,1800) 'IWFCHK/wf(1): ',labelwf, wfiwf1
      endif
c
  999 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MEMCHKIBR
c
c
      subroutine MEMCHKIBR( idbr, itypchk, labelbr, wkbroy, savbr )
c---------------------------------------------------------------
c Purpose: check the blending data
c  Complicated by the need to track IF blend data exists yet
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  wkbroy(*)
      CHARACTER*(*) labelbr
c
c Local declarations:
      DATA  initbr  / 0 /
c
      DATA  lvldbg / 0 /
c
c >>>> EXECUTABLE CODE:
c
c  If do not have in-memory broyden, bounce out ...
      if( idbr .eq. 0 ) goto 999
      if( idbr .eq. -1 )then
c       Initialize
        initbr = 0
        goto 999
      endif
c
      thisbr = wkbroy(idbr)
c
      call FLGETIWR( IWR )
c
 1800 format(a,a,1p4d10.3)
      if( itypchk .eq. 0 )then
c       Initialize
        savbr = thisbr
        initbr = 1
        if( lvldbg .gt. 0 )
     $  write(IWR,1800) 'IBR-INIT/br(1): ',labelbr, savbr
      elseif( initbr.ne.0 )then
c       There has been an initialization with history data
        if( itypchk .eq. 1 )then
c         Check against recorded value
          if( lvldbg .gt. 0 )
     $    write(IWR,1800) 'IBRCHK/br(1): ',labelbr, thisbr
          if( thisbr .ne. savbr )then
            write(IWR,1800) '***** IBRBUG *****',labelbr,savbr,thisbr
            call STOPXERR( 'Memory error - blend history overwritten' )
          endif
        elseif( itypchk .eq. 2 )then
c         Just output the value without a check
          write(IWR,1800) 'IBRCHK/br(1): ',labelbr, thisbr
        endif
      endif
c
  999 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPPMEM0
c
c
      subroutine KPPMEM0( per_process, do_eigvecs,
     $ norb,nstate,  memeig )
c---------------------------------------------------------------
c Purpose: compute memory requirements (size of wmat argument) of 
c          PSCHROEDKP for wkmem using Scalapack calls.
c
c Written:  Andrew C. Pineda,  9-February-2010 for v2.62
c Revision history:
c  4Jan13-PAS/2.64: clean-up
c---------------------------------------------------------------
c
c Return argument:
c    memeig = process value (per_process = .true.)
c           = max value over all processes (per_process = .false.)
c
      INTEGER     NB
      PARAMETER  (NB=32)
C      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
      INTEGER  norb
      LOGICAL  per_process
      LOGICAL  do_eigvecs
c
c Output arguments
      INTEGER  memeig
c
c Local variables
      integer  memeiglocal
      logical  per_process_local
c
      integer  mat, matdist
      integer  nwpreph, nwpreps, nwprep, nwsolv
      integer  nwclose, nwpsch, nwleft
      integer  iwHmat, iwOvlp, iwVecs
      integer  iwZwork, iwRwork, iwIwork, iwiwk1, iwlast, iwEnd
      integer  lwzwork,lzwork, lrwork, liwork
      integer  iwMat1, iwMat2, iwVeco
c
      integer  CONTEXT_K
      integer  nprocs_k, icomm
      integer  nprow_k, npcol_k, npsolver_k
c
      integer  ierr
c
      LOGICAL  USE_SCALAPACK_FORMULAS
      DATA    USE_SCALAPACK_FORMULAS / .true. /
C      LOGICAL  OPTIMAL_EIGENVECTORS
C      DATA    OPTIMAL_EIGENVECTORS   / .true. /
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/KPPMEM0: do_eigvecs=',do_eigvecs
      endif
      mat = norb*norb
c  Get Local MP info
c     Image level info
      call MPCOMM( icomm )
c     K-point level info
      call MPNODES_K( nprocs_k )
c
      memeiglocal = 0
c
      if( nprocs_k .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
        lwzwork = MAX( 4*norb, mat )
        lrwork  = 7*norb
        liwork  = 5*norb
c
c       Close serial solver memory needs
      else
c       USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
c
        call EIGPROCS( nprocs_k, nprow_k, npcol_k, npsolver_k )
c
c       Determine pzhegvx workspace requirements
c
C        if( USE_SCALAPACK_FORMULAS )then
Cc         Get EXACT mpi/blacs-aware work space using this context:
CC          CONTEXT_K = -1
C          call MPCOMM_K( icomm_k )
C          CONTEXT_K = icomm_k
C          call BLACS_GRIDINIT( CONTEXT_K, 'R', nprow_k, npcol_k )
C          call PSCHWKSZX( CONTEXT_K, npsolver_k, do_eigvecs,
C     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
C          lwork = lzwork
C          call BLACS_GRIDEXIT( CONTEXT_K )
C        else
c         Use default mpi/blacs-avoiding naive work space:
          call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
C        endif
      endif
c
c     Initialize memory allocations within icomm_k.
c     MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-Zwork, d-R8-Rwork)
c
cxxx264: list solver memory allocs as implemented, but not used here
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwiwk1  + norb
      iwEnd   = iwlast  - 1
c
      iwMat1  = iwOvlp
      iwMat2  = MAX( iwMat1 + mat, iwVecs )
      iwVeco  = iwZwork
c
c     The solver memory allocations have the following requirements:
      nwpreph = 2*matdist + 3*mat
      nwpreps = 4*matdist + 2*mat
      nwprep = MAX( nwpreph, nwpreps )
      nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
c     nwclose = 6*matdist + 2*norb*nstate
c     because all vectors are brought back into full norb*norb ...
cxx264: nwclose (eigeenec store) needs to be dealt with (norb*nstate(
      nwclose = 6*matdist + 2*norb*norb
      nwpsch = MAX( nwprep, nwsolv, nwclose )
c
      nwleft = nwpsch - nwsolv
      if( DEBUG_MEM )then
        write(IWR,*) 'DEV/KPPMEM0: spare solver mem=',nwleft
      endif
c
cxxx264: is this mem-req for LAPACK solve cognizant of do_eigvecs?
      memeiglocal = MAX( iwEnd, iwMat2+2*mat, iwVeco+2*mat )
c
c >>>>> Put it all together and report to final result
c
      per_process_local = per_process
c
C      if( per_process_local ) then
c        Report the local memeig value
         memeig = memeiglocal
C      else
Cc       Report the maximum value of memeig for any process in icomm
Cc       Need barrier to prevent a crash
C        call MPBARRIER( icomm )
C        icount = 1
C        call MPI_ALLREDUCE( memeiglocal, memeig, icount,
C     $       MPI_INTEGER, MPI_MAX, icomm,
C     $       ierr )
C      endif
      if( DEBUG_MEM )then
        write(IWR,'(a,4i12)') 'DEV/KPPMEM0: works,Z,R,I=',
     $                         lwzwork, lrwork, liwork
        write(IWR,'(a,3i12)') 'DEV/KPPMEM0: nwsch,memeig=',
     $                         nwpsch, memeiglocal, memeig
        write(IWR,*) 'DEV/KPPMEM0: icomm,nprocs_k=',icomm,nprocs_k
        call FLUSH( IWR )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHWKSZ0
c
c
      subroutine PSCHWKSZ0(          nprocs, do_eigvecs,
     $ norb, matdist, lwzwork, lzwork, lrwork, liwork )
c---------------------------------------------------------------
c Purpose: compute exact memory requirements for PSCHROEDKP solvers,
c          using Scalapack (p)zhegvx calls.
c      
c Written:  Peter A. Schultz, 10-January-2013, for 2.64
c           Derived from PSCHROEDKPSZ by Andrew C. Pineda (9Feb10)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Work array space for
c    (1)  serial solver zhegvx
c    (2)  parallel solver pzhegvx
c
      INTEGER     NB
      PARAMETER  (NB=32)
C      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
C      INTEGER  icontxt
      INTEGER  nprocs
      LOGICAL  do_eigvecs
c
      INTEGER  norb
c
c Output arguments
      INTEGER  matdist
      INTEGER  lwzwork, lzwork, lrwork, liwork
c
c Local variables
      INTEGER  mat
      INTEGER  nprow, npcol, npsolver
c
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/WKSZ-0: norb,nprocs=',norb,nprocs
      endif
cxx:end

      mat = norb*norb
c
      if( nprocs .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
cxxx264        lwzwork = 2*MAX( 2*norb, mat ) ! I am sure we don't need 2mat -pas
        lwzwork =   MAX( 4*norb, mat )
        lzwork  = lwzwork / 2
        lrwork  = 7*norb
        liwork  = 5*norb
cxxx264 - need to scope a do_eigvecs-aware memory for single-proc
C        memeig = 2*mat + 2*mat + 2*norb*nstate + mat + 7*norb
c
c       Close serial solver zhegvx memory needs
      else
c       Parallel, USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
        call EIGPROCS( nprocs, nprow, npcol, npsolver )
c       NOTE: following is a crude approximation to matdist, NOT exact
        call EIGMDIST( matdist, norb, npsolver )
c
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Simplified appxximate work memory requirements
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      Simplified naive values
        lwzwork = 2*matdist
        lzwork = lwzwork / 2
        lrwork = 2*matdist
        liwork = 6*norb
        if( .not. do_eigvecs )then
c         No eigenvectors, this becomes much less:
          lzwork = (NB+1) * (norb+1)
          lrwork =9*norb
        endif
c
c       Close pzhegvx parallel solver work array memory needs
      endif
c
      if( DEBUG_MEM )then
        write(IWR,'(a,3i12)') 'DEV/WKSZ-0: works, Z,R,I=',
     $                         lwzwork,lrwork,liwork
        call FLUSH( IWR )
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
