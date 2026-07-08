c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     drive the MiMiC MM client
c
      subroutine mimic_loop
      use atoms, only: n
      use domdec, only: rank
      use iounit, only: iout
      use mimic
      implicit none
#ifdef MIMIC
      integer :: request
      logical :: done
      real*8 :: energy
      real*8, allocatable :: force(:,:)
c
      if (rank == 0) write(iout,10)
   10 format (/,' Entering MiMiC QM/MM client loop')
c
      allocate(force(3,n))
      force = 0.0d0
      energy = 0.0d0
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
         else if (request == MCL_COMPUTE_FORCES) then
            call mimic_forces(energy, force)
         else if (request == MCL_SEND_PARTICLE_FORCES) then
            call mimic_send_forces(force)
         else if (request == MCL_SEND_ENERGY) then
            call mimic_send_energy(energy)
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
      deallocate(force)
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
c     ##  subroutine mimic_forces  --  MM energy+forces at current  ##
c     ##  (server-supplied) coordinates, gathered on rank 0         ##
c     ################################################################
c
c     Re-decomposes for the freshly installed coordinates, evaluates the
c     MM gradient, and reduces the global potential energy and the full
c     (3,n) gradient onto rank 0. It orchestrates the same standalone
c     routines a normal integrator step uses (reassign, gradient,
c     commforces, ...) without any integration. The returned "force"
c     array is the gradient dE/dx; the sign flip to a physical force is
c     done when sending (mimic_send_forces).
c
      subroutine mimic_forces(energy, force)
      use atoms, only: n
      use cutoff, only: use_list
      use domdec, only: nloc, glob, nbloc, COMM_TINKER
      use energi, only: esum
      use mpi
      implicit none
      real*8 :: energy
      real*8 :: force(3,n)
      real*8, allocatable :: derivs(:,:), fglob(:,:)
      real*8 :: epot
      integer :: i, iglob, ierr
      integer, save :: step = 0
c
      step = step + 1
c
c     redistribute atoms to domains for the new coordinates
c
      call reassign
      call reassignpme(.false.)
      call commpos
      call commposrec
c
      allocate(derivs(3,nbloc))
      derivs = 0.0d0
      call reinitnl(step)
      call mechanic_up_para(step)
      call allocstep
      if (use_list)  call nblist(step)
c
c     evaluate the MM energy and gradient
c
      call gradient(epot, derivs)
      call reduceen(epot)
      call commforces(derivs)
c
c     assemble the global gradient on rank 0 by summing each rank's
c     owned-atom contributions into a sparse (3,n) buffer
c
      allocate(fglob(3,n))
      fglob = 0.0d0
      do i = 1, nloc
         iglob = glob(i)
         fglob(1,iglob) = derivs(1,i)
         fglob(2,iglob) = derivs(2,i)
         fglob(3,iglob) = derivs(3,i)
      end do
      call MPI_REDUCE(fglob, force, 3*n, MPI_REAL8, MPI_SUM, 0,
     &                COMM_TINKER, ierr)
      energy = epot
c
      deallocate(fglob)
      deallocate(derivs)
      end
