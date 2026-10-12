c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     drive the MiMiC MM client
c
      subroutine mimic_loop
      use domdec, only: rank
      use iounit, only: iout
      use mimic
      use neigh, only: ineigup
      use mpi
      implicit none
#ifdef MIMIC
      integer :: request
      logical :: done
c
c     timing of the client loop, reported once at exit: the time spent
c     blocked on the next request (waiting for the server) and the time
c     spent handling each category of request
c
      integer, parameter :: tc_positions = 1, tc_field = 2,
     &   tc_energy = 3, tc_forces = 4, tc_induced = 5,
     &   tc_sendforces = 6, tc_sendenergy = 7, tc_senddipoles = 8,
     &   tc_sendmultipoles = 9, tc_other = 10, ntcat = 10
      character(len=36), parameter :: tclabel(ntcat) = (/
     &   'receive positions + re-decomposition',
     &   'receive electric field              ',
     &   'compute energy                      ',
     &   'compute forces                      ',
     &   'compute induced multipoles          ',
     &   'send forces                         ',
     &   'send energies (total, induction)    ',
     &   'send induced dipoles                ',
     &   'send multipoles                     ',
     &   'other requests (setup, info, exit)  ' /)
      real*8 :: tcat(ntcat), twait, twall, tstart, t0, t1, t2
      integer :: ncount(ntcat), nreq, icat
c
      if (rank == 0) write(iout,10)
   10 format (/,' Entering MiMiC QM/MM client loop')
c
c     reject the polarization solvers that ignore the external field
c
      call mimic_check_polar
c
c     the client does not know the server's time step, so rebuild the
c     neighbor lists for every new geometry instead of every ineigup
c     steps of a placeholder dt
c
      ineigup = 1
c
c     report the client loop only: mdinit already ran a full gradient,
c     which filled the per-term and PME timers
c
      call clear_timers
      tcat = 0.0d0
      ncount = 0
      twait = 0.0d0
      nreq = 0
      tstart = mpi_wtime()
c
      done = .false.
      do while (.not.done)
         t0 = mpi_wtime()
         call mimic_get_request(request)
         t1 = mpi_wtime()
         twait = twait + t1 - t0
         nreq = nreq + 1
c
         if (request == MCL_SEND_CLIENT_ID) then
            call mimic_send_client_id
         else if (request == MCL_SEND_CLIENT_NAME) then
            call mimic_send_client_name
         else if (request == MCL_SEND_CLIENT_RUNTYPE) then
            call mimic_send_client_runtype
         else if (request == MCL_SEND_CLIENT_APIVER) then
            call mimic_send_client_apiver
         else if (request == MCL_SEND_BOX_VECTORS) then
            call mimic_send_box
         else if (request == MCL_SEND_NUM_PARTICLES) then
            call mimic_send_natoms
         else if (request == MCL_SEND_NUM_PARTICLE_SPECIES) then
            call mimic_send_num_species
         else if (request == MCL_SEND_PARTICLE_SPECIES_IDS) then
            call mimic_send_species_ids
         else if (request == MCL_SEND_SPECIES_ELEMENTS) then
            call mimic_send_species_elements
         else if (request == MCL_SEND_SPECIES_MASSES) then
            call mimic_send_species_masses
         else if (request == MCL_SEND_SPECIES_LABELS) then
            call mimic_send_species_labels
         else if (request == MCL_SEND_NUM_FRAGMENTS) then
            call mimic_send_num_fragments
         else if (request == MCL_SEND_NUM_PARTICLES_IN_FRAGMENTS) then
            call mimic_send_fragment_natoms
         else if (request == MCL_SEND_PARTICLE_IDS_IN_FRAGMENTS) then
            call mimic_send_fragment_ids
         else if (request == MCL_SEND_PARTICLE_MULTIPOLE_ORDER) then
            call mimic_send_multipole_order
         else if (request == MCL_SEND_PARTICLE_MULTIPOLES) then
            call mimic_send_multipoles
         else if (request == MCL_SEND_UPDATED_PARTICLE_MULTIPOLES) then
            call mimic_send_multipoles
         else if (request == MCL_SEND_NUM_CONSTR_BONDS) then
            call mimic_send_num_constr_bonds
         else if (request == MCL_SEND_CONSTR_BOND_PARTICLE_IDS) then
            call mimic_send_constr_bond_ids
         else if (request == MCL_SEND_CONSTR_BOND_LENGTHS) then
            call mimic_send_constr_bond_lengths
         else if (request == MCL_SEND_NUM_CONSTR_ANGLES) then
            call mimic_send_num_constr_angles
         else if (request == MCL_SEND_PARTICLE_POSITIONS) then
            call mimic_send_positions
         else if (request == MCL_SEND_NUM_BOUNDARY_BONDS) then
            call mimic_send_num_boundary_bonds
         else if (request == MCL_SEND_BOUNDARY_BOND_ATOMS) then
            call mimic_send_boundary_bonds
         else if (request == MCL_RECV_PARTICLE_POSITIONS) then
            call mimic_recv_positions
            call mimic_prepare
         else if (request == MCL_RECV_EL_FIELD_ON_PARTICLE) then
            call mimic_recv_field
         else if (request == MCL_COMPUTE_ENERGY) then
            call mimic_energy
         else if (request == MCL_COMPUTE_FORCES) then
            call mimic_forces
         else if (request == MCL_COMPUTE_INDUCED_MULTIPOLES) then
            call mimic_induced_multipoles
         else if (request == MCL_SEND_PARTICLE_FORCES) then
            call mimic_send_forces
         else if (request == MCL_SEND_ENERGY) then
            call mimic_send_energy
         else if (request == MCL_SEND_INDUCTION_ENERGY) then
            call mimic_send_induction_energy
         else if (request == MCL_SEND_PARTICLE_IND_DIPOLES_AVG) then
            call mimic_send_induced_dipoles('averaged')
         else if (request == MCL_SEND_PARTICLE_IND_DIPOLES) then
            call mimic_send_induced_dipoles('d')
         else if (request == MCL_EXIT) then
            done = .true.
         else
            if (rank == 0) write(iout,20) request
   20       format (/,' MIMIC_LOOP  --  unhandled request',i8, !< cAA add MCL translate
     &                ', leaving client loop')
            done = .true.
         end if
c
         t2 = mpi_wtime()
         icat = timing_category(request)
         tcat(icat) = tcat(icat) + t2 - t1
         ncount(icat) = ncount(icat) + 1
      end do
      twall = mpi_wtime() - tstart
c
      call mimic_print_timings
c
      call mimic_finalize
      if (rank == 0) write(iout,30)
   30 format (/,' Leaving MiMiC QM/MM client loop')
      return
c
      contains
c
c     map a request tag to its timing category; the MCL tags are
c     variables, not constants, hence the if chain
c
      integer function timing_category(req)
      integer :: req
      if (req == MCL_RECV_PARTICLE_POSITIONS) then
         timing_category = tc_positions
      else if (req == MCL_RECV_EL_FIELD_ON_PARTICLE) then
         timing_category = tc_field
      else if (req == MCL_COMPUTE_ENERGY) then
         timing_category = tc_energy
      else if (req == MCL_COMPUTE_FORCES) then
         timing_category = tc_forces
      else if (req == MCL_COMPUTE_INDUCED_MULTIPOLES) then
         timing_category = tc_induced
      else if (req == MCL_SEND_PARTICLE_FORCES) then
         timing_category = tc_sendforces
      else if (req == MCL_SEND_ENERGY .or.
     &         req == MCL_SEND_INDUCTION_ENERGY) then
         timing_category = tc_sendenergy
      else if (req == MCL_SEND_PARTICLE_IND_DIPOLES_AVG) then
         timing_category = tc_senddipoles
      else if (req == MCL_SEND_PARTICLE_MULTIPOLES .or.
     &         req == MCL_SEND_UPDATED_PARTICLE_MULTIPOLES) then
         timing_category = tc_sendmultipoles
      else
         timing_category = tc_other
      end if
      end function
c
c     whole-run timing report, averaged over the ranks like mdstattime.
c     The MiMiC table is printed with the VERBOSE keyword. Tinker's
c     standard tables follow when forces were computed, selected as
c     usual by the VERBOSE, VERBOSEFORCESTIME and VERBOSEPMETIME
c     keywords: they are averaged per MCL_COMPUTE_FORCES request, and
c     their percentages refer to the wall time of the client loop.
c     Collective; verbose comes from the key file, so every rank takes
c     the same branch
c
      subroutine mimic_print_timings
      use domdec, only: nproc, COMM_TINKER
      use inform, only: iprint, verbose, verboseforcestime,
     &                  verbosepmetime
      use timestat, only: timestep
      real*8 :: buf(ntcat+2), pct
      integer :: i, ierr
c
      if (verbose) then
         buf(1:ntcat) = tcat
         buf(ntcat+1) = twait
         buf(ntcat+2) = twall
         if (rank == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, buf, ntcat+2, MPI_REAL8,
     &                      MPI_SUM, 0, COMM_TINKER, ierr)
         else
            call MPI_REDUCE(buf, buf, ntcat+2, MPI_REAL8, MPI_SUM, 0,
     &                      COMM_TINKER, ierr)
         end if
c
         if (rank == 0) then
            buf = buf / dble(nproc)
            write(iout,110) nproc
  110       format (/,' ',72('='),/,' MiMiC client timings over the',
     &                ' whole run, averaged over',i5,' rank(s)',
     &              /,' ',72('-'),/,' request',30x,'count',
     &                '     total (s)   average (s)  % wall')
            call print_row('waiting for server requests         ',
     &                     nreq, buf(ntcat+1), buf(ntcat+2))
            do i = 1, ntcat
               call print_row(tclabel(i), ncount(i), buf(i),
     &                        buf(ntcat+2))
            end do
            pct = 0.0d0
            if (buf(ntcat+2) > 0.0d0)  pct = 100.0d0
     &         * (buf(ntcat+2) - buf(ntcat+1)) / buf(ntcat+2)
            write(iout,120) buf(ntcat+2), buf(ntcat+2) - buf(ntcat+1),
     &                      pct
  120       format (' ',72('-'),
     &            /,' wall time of the client loop',15x,f14.6,
     &            /,' busy (wall - waiting)',22x,f14.6,14x,f8.2,
     &            /,' ',72('='))
         end if
      end if
c
c     Tinker's standard tables, one block over all force requests
c
      if (ncount(tc_forces) > 0) then
         if (rank == 0 .and. (verbose .or. verboseforcestime .or.
     &       verbosepmetime))  write(iout,130) ncount(tc_forces)
  130    format (/,' MiMiC: the averages below are per',
     &             ' MCL_COMPUTE_FORCES request (',i0,' requests)',
     &           /,' and the percentages refer to the wall time of',
     &             ' the client loop')
         timestep = twall
         iprint = ncount(tc_forces)
         call mdstattime(ncount(tc_forces), 0.0d0)
      end if
      end subroutine
c
      subroutine clear_timers
      use timestat
      timeinte = 0.0d0
      timereneig = 0.0d0
      timecommpos = 0.0d0
      timeparam = 0.0d0
      timenl = 0.0d0
      timegrad = 0.0d0
      timered = 0.0d0
      timetp = 0.0d0
      timecommforces = 0.0d0
      timestep = 0.0d0
      timebonded = 0.0d0
      timevdw = 0.0d0
      timeelec = 0.0d0
      timepolar = 0.0d0
      timecleargrad = 0.0d0
      timereal = 0.0d0
      timerec = 0.0d0
      timecommforcesrec = 0.0d0
      timecommforcesreal = 0.0d0
      timegrid = 0.0d0
      timefft = 0.0d0
      timescalar = 0.0d0
      timegrid2 = 0.0d0
      timerecreccomm = 0.0d0
      end subroutine
c
      subroutine print_row(label, cnt, ttot, tw)
      character(len=*) :: label
      integer :: cnt
      real*8 :: ttot, tw, tavg, pct
      tavg = 0.0d0
      if (cnt > 0)  tavg = ttot / dble(cnt)
      pct = 0.0d0
      if (tw > 0.0d0)  pct = 100.0d0 * ttot / tw
      write(iout,140) label, cnt, ttot, tavg, pct
  140 format (' ',a36,i8,f14.6,f14.6,f8.2)
      end subroutine
#else
      write(iout,40)
   40 format (/,' MIMIC_LOOP  --  this binary was not built with',
     &          ' MiMiC support; reconfigure with --enable-mimic')
      call fatal
#endif
      end
c
c     ################################################################
c     ##  subroutine mimic_prepare  --  re-decompose for the newly  ##
c     ##  installed (server-supplied) coordinates                   ##
c     ################################################################
c
c     Runs the same standalone routines a normal integrator step uses to
c     redistribute atoms over domains and rebuild the neighbor lists,
c     without any integration. All of it depends on the coordinates
c     only, so it runs once per MCL_RECV_PARTICLE_POSITIONS, not per
c     compute request; mdinit has already prepared the initial ones.
c     The phases are timed into the same timestat buckets as in verlet.
c
      subroutine mimic_prepare
      use cutoff, only: use_list
      use timestat
      use mpi
      implicit none
      integer, save :: step = 0
      real*8 :: time0, time1
c
      step = step + 1
c
      time0 = mpi_wtime()
      call reassign
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
      time0 = mpi_wtime()
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
c
      time0 = mpi_wtime()
      call reinitnl(step)
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
      time0 = mpi_wtime()
      call mechanic_up_para(step)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
      time0 = mpi_wtime()
      if (use_list)  call nblist(step)
      time1 = mpi_wtime()
      timenl = timenl + time1 - time0
      end
c
c     ################################################################
c     ##  subroutine mimic_energy  --  MM energy and induced        ##
c     ##  dipoles at the current coordinates, no forces             ##
c     ################################################################
c
c     Answers MCL_COMPUTE_ENERGY. The energy stays rank-local; it is
c     reduced by whichever send request asks for it.
c
      subroutine mimic_energy
      use mimic, only: mimic_stash_energy
      implicit none
      real*8 energy
      real*8 :: epot
c
      epot = energy ()
      call mimic_stash_energy(epot)
      end
c
c     ################################################################
c     ##  subroutine mimic_induced_multipoles  --  induced dipoles  ##
c     ##  only, at the current coordinates and external field       ##
c     ################################################################
c
c     Answers MCL_COMPUTE_INDUCED_MULTIPOLES. Solves the polarization
c     equations alone, without the energy or the gradient. The steps
c     energy() would run before the solver are repeated here: empole
c     normally rebuilds rpole, which epolar skips under use_mpole. The
c     external field and the MiMiC-group zeroing are applied inside the
c     solvers, which start from the previous dipoles (polgsf = 1).
c
      subroutine mimic_induced_multipoles
      use ewald, only: aewald, apewald
      use mimic, only: mimic_stash_dipoles
      use mpole, only: npole
      use potent, only: use_polar, use_pmecore, use_chgflx
      implicit none
c
      if (use_polar .and. npole > 0) then
         if (use_chgflx)  call alterchg
         aewald = apewald
         call chkpole(.false.)
         call rotpole
c
c     the solvers rejected by mimic_check_polar never reach this point
c
         if (use_pmecore) then
            call newinduce_pme
         else
            call newinduce_pme2
         end if
      end if
c
      call mimic_stash_dipoles
      end
c
c     ################################################################
c     ##  subroutine mimic_forces  --  MM energy, induced dipoles   ##
c     ##  and gradient at the current coordinates                   ##
c     ################################################################
c
c     Answers MCL_COMPUTE_FORCES. Both the energy and the gradient stay
c     rank-local; the gradient is scattered into a sparse global buffer
c     and reduced only when the forces are requested. The stash holds
c     dE/dx; the sign flip to a physical force is done when sending.
c
      subroutine mimic_forces
      use deriv, only: resetForcesRec
      use domdec, only: nbloc
      use mimic, only: mimic_stash_energy, mimic_stash_forces
      use timestat
      use mpi
      implicit none
      real*8, allocatable :: derivs(:,:)
      real*8 :: epot, time0, time1
c
      allocate(derivs(3,nbloc))
      derivs = 0.0d0
c
c     gradient() does not clear the reciprocal-space accumulators, and
c     allocstep only runs when new positions arrive, so a repeated
c     request at the same geometry would otherwise add to the last one
c
      time0 = mpi_wtime()
      call resetForcesRec
c
c     evaluate the MM energy and gradient
c
      call gradient(epot, derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
c
      time0 = mpi_wtime()
      call commforces(derivs)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1 - time0
c
      call mimic_stash_forces(derivs)
      call mimic_stash_energy(epot)
c
      deallocate(derivs)
      end
