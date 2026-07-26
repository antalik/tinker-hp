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
      implicit none
#ifdef MIMIC
      integer :: request
      logical :: done
c
      if (rank == 0) write(iout,10)
   10 format (/,' Entering MiMiC QM/MM client loop')
c
c     reject the polarization solvers that ignore the external field
c
      call mimic_check_polar
c
      done = .false.
      do while (.not.done)
         call mimic_get_request(request)
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
         else if (request == MCL_RECV_EL_FIELD_ON_PARTICLE) then
            call mimic_recv_field
         else if (request == MCL_COMPUTE_ENERGY) then
            call mimic_energy
         else if (request == MCL_COMPUTE_FORCES) then
            call mimic_forces
         else if (request == MCL_SEND_PARTICLE_FORCES) then
            call mimic_send_forces
         else if (request == MCL_SEND_ENERGY) then
            call mimic_send_energy
         else if (request == MCL_SEND_INDUCTION_ENERGY) then
            call mimic_send_induction_energy
         else if (request == MCL_SEND_PARTICLE_INDUCED_DIPOLES) then
            call mimic_send_induced_dipoles
         else if (request == MCL_EXIT) then
            done = .true.
         else
            if (rank == 0) write(iout,20) request
   20       format (/,' MIMIC_LOOP  --  unhandled request',i8,
     &                ', leaving client loop')
            done = .true.
         end if
      end do
c
      call mimic_finalize
      if (rank == 0) write(iout,30)
   30 format (/,' Leaving MiMiC QM/MM client loop')
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
c     without any integration.
c
      subroutine mimic_prepare
      use cutoff, only: use_list
      implicit none
      integer, save :: step = 0
c
      step = step + 1
c
      call reassign
      call reassignpme(.false.)
      call commpos
      call commposrec
c
      call reinitnl(step)
      call mechanic_up_para(step)
      call allocstep
      if (use_list)  call nblist(step)
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
      call mimic_prepare
      epot = energy ()
      call mimic_stash_energy(epot)
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
      use domdec, only: nbloc
      use mimic, only: mimic_stash_energy, mimic_stash_forces
      implicit none
      real*8, allocatable :: derivs(:,:)
      real*8 :: epot
c
      call mimic_prepare
c
      allocate(derivs(3,nbloc))
      derivs = 0.0d0
c
c     evaluate the MM energy and gradient
c
      call gradient(epot, derivs)
      call commforces(derivs)
c
      call mimic_stash_forces(derivs)
      call mimic_stash_energy(epot)
c
      deallocate(derivs)
      end
