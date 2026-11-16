c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DOENGY
c
c
      subroutine DOENGY( engytotl, IWR,
     $ idmatfl,idmatsfl, iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $ nspin, nearopt,near2c,
     $ atm0engy,esnsgr,vesgrxc,excengy,esnsengy, ves1ctr,ves1c0,
     $ do_field, engydfld,engy0fld, engypol,
     $ xc0val,exctot,excslo,exc1c, corslow,
     $ norb, nk, natm,ntyp,nshld,
     $ itypa,numshl,lshel, znuc,occij, wtk,
     $ rhoatom, defatom,defatdn,
     $ vslnear,vesnear,excnear, norbd,
     $ vmat,          dmat, h0mat )
c      nmatloc*nspin  mat   mat   ; nmatloc=mat*nkloc
c---------------------------------------------------------------
c Purpose: perform energy analysis for this iteration
c
c Written: Peter A. Schultz, 7-June-1994 for v2.10
c
c Revision history:
c  13Feb12-PAS&ACP/2.63: k-parallel
c  22Sep08-PAS/2.62: split e/vxc0, h0t/nl files, dmat/spin file
c  20Jun07-PAS/2.60: merge serial and tp
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  30Jul05-PAS/2.59-tp: master-only matrix files
c  17Sep03-PAS,JTK/2.57: install E-field terms; bulk polarization
c   3Oct02-PAS/2.54: rework storage of ham0 pieces
c  21Jul01-PAS/2.49: compact dmat; spin-polarized dft
c  21May01-PAS/2.47: extracted traces for (non-stardard) 1c es
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  16Sep98-PAS/2.27: tweak the core nearby diagnostic
c  14Jul98-PAS/2.23: change energy to rid delh trace; split basis
c  26Oct97-PAS/2.21: clean up exc calc with output levels
c   3Jun97-PAS/2.17: naming conventions / do0trace used for es
c                    another exc correction
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      LOGICAL    do_field
      DIMENSION  xc0val(*)
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  znuc(ntyp), occij(nshld,nshld,ntyp), wtk(nk)
      DIMENSION  rhoatom(natm,*)
      DIMENSION  defatom(natm,*), defatdn(natm,*)
      DIMENSION  vslnear(norbd,nspin),vesnear(natm)
      DIMENSION  excnear(norbd,nspin)
c
C      DIMENSION  vmat(norb*norb,nkloc,nspin)
      DIMENSION  vmat(norb*norb,*)
c
c Scratch arrays:
      DIMENSION  dmat(norb*norb), h0mat(norb*norb)
c
c Local declarations:
      DIMENSION  vtrk(2)
      DATA  zero,half,one / 0.d0,0.5d0,1.d0 /
      DATA  lstout / 1 /
c
c Local/task-parallel data
      DIMENSION  traces(8)
c
c >>>> EXECUTABLE CODE:
c
 9010 format( 4x, a, (1x,3f20.10) )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  First compute the contributions from the nearby corrections
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      spinfac = one
      if( nspin.eq.2 ) spinfac = half
      vesnr0 = zero
      vesnr1 = zero
      vesnr2 = zero
      excnr1 = zero
      excnr2 = zero
      excnrc = zero
      slonr1 = zero
      slonr2 = zero
c
      do 101 iatm=1,natm
        defatm0 = znuc(itypa(iatm)) - rhoatom(iatm,1)
        vesnr0 = vesnr0 + vesnear(iatm)*defatm0
c
        vesnr1 = vesnr1 + vesnear(iatm)*defatom(iatm,1)
        vesnr2 = vesnr2 + vesnear(iatm)*defatom(iatm,2)
        excnr1 = excnr1 + excnear(iatm,1)*defatom(iatm,1)
        excnr2 = excnr2 + excnear(iatm,1)*defatom(iatm,2)
        slonr1 = slonr1 + vslnear(iatm,1)*defatom(iatm,1)
        slonr2 = slonr2 + vslnear(iatm,1)*defatom(iatm,2)
c
        defcore = rhoatom(iatm,2) - rhoatom(iatm,3)
        excnrc = excnrc + excnear(iatm,1)*defcore*spinfac
  101 continue
c
      if( nspin .eq. 2 )then
        do 102 iatm=1,natm
          vesnr1 = vesnr1 + vesnear(iatm)*defatdn(iatm,1)
          vesnr2 = vesnr2 + vesnear(iatm)*defatdn(iatm,2)
          excnr1 = excnr1 + excnear(iatm,2)*defatdn(iatm,1)
          excnr2 = excnr2 + excnear(iatm,2)*defatdn(iatm,2)
          slonr1 = slonr1 + vslnear(iatm,2)*defatdn(iatm,1)
          slonr2 = slonr2 + vslnear(iatm,2)*defatdn(iatm,2)
c
          defcore = rhoatom(iatm,2) - rhoatom(iatm,3)
          excnrc = excnrc + excnear(iatm,2)*defcore*spinfac
  102   continue
      endif
c
      if( near2c.eq.0 )then
        slonr = zero
        excnr = zero
        vesnr = zero
      elseif( near2c.eq.1 )then
        slonr = slonr1
        excnr = excnr1
        vesnr = vesnr1
      elseif( near2c.eq.2 )then
        slonr = slonr1 + slonr2
        excnr = excnr1 + excnr2
        vesnr = vesnr1 + vesnr2
      endif
      esnsnr = half*(vesnr-vesnr0)
c
      if( lstout.gt.2 )then
        write(IWR,*) 'DOENGY: nearby construction energy analysis'
        write(IWR,910)
  910   format(1x,'E-nearby',9x,'total',10x,'1-ctr',10x,'2-ctr')
  911   format(2x,a6,2x,4f15.8)
        write(IWR,911) 'v-near',slonr,slonr1,slonr2
        write(IWR,911) 'exc-nr',excnr,excnr1,excnr2,excnrc
        write(IWR,911) 'ves-nr',vesnr,vesnr1,vesnr2
        write(IWR,911) 'ves0nr',vesnr0,vesnr0,zero,esnsnr
        write(IWR,*)
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Compute contributions from variety of traces
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      matsiz = norb*norb
      vtrc    = zero
      vtrcdn  = zero
      pseutrc = zero
      ekintrc = zero
      vxc0trc = zero
      exc0trc = zero
      h0trc   = zero
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c     And need nkloc to compute array offset for ispin!=1
      call KPMINE ( nk, nkloc, nk0 )
c
c  Only masters read files, do traces, and run the big calcs
      if( iproc .ne. master .and. iproc_k .ne. master_k ) goto 300
c     Skip if no local k-points, EXCEPT the image master who must read
      if( nkloc .eq. 0 .and. iproc .ne. master ) goto 300
c
      if( iproc .eq. master )then
        REWIND( unit=idmatfl )
        if( nspin .eq. 2 )then
c         We have spin-dependent potentials and density matrices,
c          position the spin-dn dmats
          REWIND( unit=idmatsfl )
        endif
        REWIND( unit=iham0fl )
        REWIND( unit=ivxc0fl )
        REWIND( unit=iexc0fl )
        REWIND( unit=ih0tfl )
        REWIND( unit=ih0nlfl )
      endif
c
      ikloc = 0
      do 250 ik=1,nk
c
        ispin = 1
c       Inform the k-parallel:
        lenmsg = matsiz
        call KPMASTER( nk, nprocs, ik, kmasterpr )
c       NB: kmasterpr = rank of master of ik, in image comm, DIFFERENT from
c           master_k = rank k-point master of this proc, in comm_k
c
c       Retrieve dmat(k;up) and send to its k-master
c       (This might need to be managed via tickets)
c
        if( iproc .eq. master )then
c         Image master reads in density matrix for this k-vector
          call READBIG( idmatfl, matsiz, dmat )
          if( kmasterpr .ne. master )then
c           ... and sends this dmat off to its k-master
            call MPSENDR8( kmasterpr, lenmsg, dmat , icomm )
            goto 220
          endif
        elseif( iproc .eq. kmasterpr )then
c         The k-master receives its dmat from master
          call MPRECVR8( master, lenmsg, dmat , icomm )
        else
c         This processor not involved with this k-point, skip
          goto 250
        endif
c
        ikloc = ikloc + 1
c
c       Trace with the delta-Ham
        vtrk(1) = zero
c       call DOTRACE( vtrk(1), norb, dmat, vmat(1,ikloc,ispin=1) )
        call DOTRACE( vtrk(1), norb, dmat, vmat(1,ikloc) )
c       The next line only matters on the image master
        vtrc = vtrc + vtrk(1)
c
  220   continue
c
        if( nspin .eq. 2 )then
c         The above was up-spin(in dmat), now must do down-spin(in h0mat):
          ispin = 2
          lenmsg = matsiz
          if( iproc .eq. master )then
c           Image master reads in dn-spin density matrix for this k-vector
            call READBIG( idmatsfl, matsiz, h0mat )
            if( kmasterpr .ne. master )then
c             ... and sends this dmat off to its k-master
              call MPSENDR8( kmasterpr, lenmsg, h0mat, icomm )
              goto 240
            endif
          elseif( iproc .eq. kmasterpr )then
c           The k-master receives its dmat from master
            call MPRECVR8( master, lenmsg, h0mat, icomm )
          else
c           This processor not involved in this k-point, skip
            goto 250
          endif
c
c         Trace with the dn-spin delta-Ham
          vtrk(2) = zero
c         call DOTRACE( vtrk(2), norb, h0mat, vmat(1,ikloc,ispin=2) )
          call DOTRACE( vtrk(2), norb, h0mat, vmat(1,ikloc+nkloc) )
c
c         The next line only matters on the image master
          vtrcdn = vtrcdn + vtrk(2)
        endif
c
  240   continue
c
c       At this point, traces with del-Hams are done,
c       the up-dmat is in dmat(), and the dn-dmat (if spin) is in h0mat.
c       Now need to do traces with iter-independent Ham components.
c       These are done only by image master, which has Ham0s on disk.
c
        lentrc = nspin
        if( iproc .eq. master )then
c
c         Combine up- and dn-spin density matrices:
          if( nspin.eq.2 ) call DAXPY( matsiz, one, h0mat,1, dmat,1 )
c
c         Full spherical atom trace term from matrix on ihamsfl:
          call READBIG( iham0fl, matsiz, h0mat )
          call DOTRACE( h0trc  , norb, dmat, h0mat )
c
c         Cutoff spherical-atom Vxc term from matrix on ihamsfl:
          call READBIG( ivxc0fl, matsiz, h0mat )
          call DOTRACE( vxc0trc, norb, dmat, h0mat )
c
c         Cutoff spherical-atom Exc term from matrix on ihamsfl:
          call READBIG( iexc0fl, matsiz, h0mat )
          call DOTRACE( exc0trc, norb, dmat, h0mat )
c
c         Kinetic energy from matrix on ihamsfl:
          call READBIG( ih0tfl , matsiz, h0mat )
          call DOTRACE( ekintrc, norb, dmat, h0mat )
c
c         Non-local pseudopotential energy from matrix on ihamsfl:
          call READBIG( ih0nlfl , matsiz, h0mat )
          call DOTRACE( pseutrc, norb, dmat, h0mat )
c
          if( master .ne. kmasterpr )then
c           Get del-Ham traces from k-master (if it is not master)
            call MPRECVR8( kmasterpr, lentrc, vtrk, icomm )
            vtrc = vtrc + vtrk(1)
            if( nspin .eq. 2 ) vtrcdn = vtrcdn + vtrk(2)
          endif
c
c         End master-only component traces
        elseif( iproc .eq. kmasterpr )then
c         The k-master sends its del-Ham traces to image master
c
c         The following if-statement is redundant with this if-black
          if( master .ne. kmasterpr )then
            call MPSENDR8( master, lentrc, vtrk, icomm )
          endif
c
        endif
c
c       End of k-point loop acquiring dmat(k) and doing traces
  250 continue
c
      vtrc = vtrc + vtrcdn
c
      traces(1) = vtrc
      traces(2) = h0trc
      traces(3) = vxc0trc
      traces(4) = exc0trc
      traces(5) = ekintrc
      traces(6) = pseutrc
c
c  Non-master processors return here ...
  300 continue
c
c     Broadcast result to everyone from master (synchronization step)
c
      lentrc = 6
      if( nprocs.gt.1 ) call MPBCAST8( master, lentrc, traces, icomm )
      vtrc    = traces(1)
      h0trc   = traces(2)
      vxc0trc = traces(3)
      exc0trc = traces(4)
      ekintrc = traces(5)
      pseutrc = traces(6)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Put all the terms together and output analysis
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c non-standard 1c es term:
      esns1c = half*(ves1ctr-ves1c0)
      vxcengy = vxc0trc + vtrc - (vesgrxc + ves1ctr + vesnr)
      esnsengy = esnsgr + esns1c + esnsnr + engydfld
      engnoxc = h0trc - vxc0trc + atm0engy + esnsengy
c
c Break down details of xc energy evaluation
      exc0ptnr = xc0val(3)
      exc0wtnr = xc0val(4)
      cor0fix = -xc0val(5) + xc0val(6)
      corfast =  xc0val(7)
      exc0fix = -xc0val(1) + xc0val(2) + cor0fix
      exc0slo =  xc0val(8) - xc0val(9)
c
      if( nearopt.eq.1 )then
        exc0nr = exc0ptnr
      elseif( nearopt.eq.2 )then
        exc0nr = exc0wtnr
      else
        exc0nr = zero
      endif
      excnr0 = excnr - exc0nr
      eslofix = exc0slo - exc0nr
c
c Compute total xc energy using different schemes:
      excraw = exctot
      excpas = exc0trc + excslo + exc1c + excnr + corfast + cor0fix
      excnew = excpas + exc0slo - exc0nr
c
c Select xc energy to use, and construct total energy:
      excengy = excnew
      engytotl = engnoxc + excengy
c
      ecoulo   = engnoxc - ekintrc - pseutrc
c
c More detailed breakdown of xc-energy as diagnostics:
      if( lstout.gt.3 )then
        write(IWR,9010) 'VXCtrc0,vtrc   =',vxc0trc,vtrc
        write(IWR,9010) 'VESves1c,1c0,es=',ves1ctr,ves1c0,esns1c
        write(IWR,9010) 'VXCvesgr,1c,nr =',vesgrxc,ves1ctr,vesnr
        write(IWR,9010) ' esns(sl,1c,nr)=',esnsgr,esns1c,esnsnr
        write(IWR,*)'$$$$$$$$$$$$$$$$$$$ EXC ANALYSIS $$$$$$$$$$$$$$$$$'
        write(IWR,9010) 'EXC0(gr,at,fix)= ',xc0val(1),xc0val(2),exc0fix
        write(IWR,9010) 'CORE(fst,sl,fix)=',corfast,corslow,cor0fix
      endif
      etotraw = engnoxc + excraw
      etotpas = engnoxc + excpas
      etotnew = engnoxc + excnew
      if( lstout.gt.3 )then
        if( nearopt.eq.1 )then
          write(IWR,9010) 'EXCnr(0,it,dif)=',exc0ptnr,excnr,excnr0
        elseif( nearopt.eq.2 )then
          write(IWR,9010) 'EXCwt(0,it,dif)=',exc0wtnr,excnr,excnr0
        else
          excnr0 = zero
        endif
        write(IWR,9010) 'EXC0slo,nr,dif=',exc0slo,exc0nr,eslofix
        write(IWR,9010) 'EXCgrd,1c,nr0=', exctot, exc1c, excnr0
        write(IWR,9010) 'EXCtrc,slo,nr=', exc0trc,excslo,excnr
        write(IWR,9010)' **** total energy without EXC =       ',engnoxc
        write(IWR,9010) 'EXC(raw,pas,new) =',excraw, excpas, excnew
        write(IWR,9010) 'Etot with EXC =   ',etotraw,etotpas,etotnew
        write(IWR,*) ' '
      endif
c     Note that "corslow"=INT(delexc*rhocor) contained in excslo
c
c Output full energy analysis:
c
      write(IWR,*)
     $'######################  ENERGY ANALYSIS  #######################'
      if( lstout.gt.0 )then
        write(IWR,9010)'Spherical-atom es ENERGY correction =', atm0engy
        write(IWR,9010)'Full spherical atom trace ENERGY =   ', h0trc
        write(IWR,9010)'Spherical-atom-Vxc trace ENERGY =    ',-vxc0trc
        write(IWR,9010)'non-spherical es ENERGY =            ', esnsengy
        write(IWR,9010)'Total Exc ENERGY =                   ', excengy
        write(IWR,*)
        write(IWR,9010)'total kinetic trace ENERGY =   ', ekintrc
        write(IWR,9010)'non-local ps-pot trace ENERGY =', pseutrc
        write(IWR,9010)'total electrostatic ENERGY =   ', ecoulo
        write(IWR,9010)'Total Vxc ENERGY =             ', vxcengy
        if( do_field )then
        es_nofield = esnsengy - engydfld
        write(IWR,9010)'non-spherical ENERGY w/o field =', es_nofield
        write(IWR,9010)'Total electric-field ENERGY =   ', engydfld
        endif
      endif
      write(IWR,*)
      write(IWR,9010)'TOTAL ENERGY (Rydberg) =             ', engytotl
      if( ABS( engypol ) .gt. 1.d-6 )then
        write(IWR,9010) 'Bulk polarization ENERGY(Ry) =',-ABS( engypol )
        write(IWR,9010) 'TOTAL ENERGY w/bulk polarization=',
     $    engytotl - ABS( engypol )
      endif
      write(IWR,*)
      write(IWR,*)
     $'################################################################'
c
      RETURN
c
c    That's all Folks!
c
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
