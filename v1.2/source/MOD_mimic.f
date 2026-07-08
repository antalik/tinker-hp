c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mimic  --  MiMiC interface (MCL client) wrap         ##
c     ##                                                              ##
c     ##################################################################
c
c     Thin wrapper around the MCL for the Tinker-HP MM client.
c
c     use_mimic   .true. when running as a MiMiC MM client
c
c     mimic_init hands COMM_WORLD to MCL right after MPI is set up. MCL
c     overwrites the communicator in place with this client's sub-
c     communicator, so COMM_TINKER and the domain decomposition are then
c     built on the MM subset. Compiled with the real MCL calls only when
c     configured with --enable-mimic (-DMIMIC); otherwise the routines
c     are stubs. mcl_initialize always runs in a MiMiC-enabled build;
c     use_mimic (client behaviour) is switched on only when the MIMIC
c     integrator is selected (keyword INTEGRATOR MIMIC).
c
      module mimic
#ifdef MIMIC
      use mcl
#endif
      implicit none
      private
      public :: use_mimic
      public :: mimic_init
      public :: mimic_get_request
      public :: mimic_finalize
      public :: mimic_abort
      public :: mimic_send_client_id
      public :: mimic_send_client_name
      public :: mimic_send_client_runtype
      public :: mimic_send_client_apiver
      public :: mimic_send_box
      public :: mimic_send_natoms
      public :: mimic_send_num_species
      public :: mimic_send_species_ids
      public :: mimic_send_species_elements
      public :: mimic_send_species_masses
      public :: mimic_send_species_labels
      public :: mimic_send_num_fragments
      public :: mimic_send_fragment_natoms
      public :: mimic_send_fragment_ids
      public :: mimic_send_multipole_order
      public :: mimic_send_multipoles
      public :: mimic_send_num_constr_bonds
      public :: mimic_send_constr_bond_ids
      public :: mimic_send_constr_bond_lengths
      public :: mimic_send_num_constr_angles
      public :: mimic_send_positions
      public :: mimic_send_num_boundary_bonds
      public :: mimic_send_boundary_bonds
      public :: mimic_recv_positions
      public :: mimic_send_forces
      public :: mimic_send_energy
#ifdef MIMIC
c     re-exported request tags so the client loop can dispatch on them
      public :: MCL_EXIT
      public :: MCL_SEND_CLIENT_ID
      public :: MCL_SEND_CLIENT_NAME
      public :: MCL_SEND_CLIENT_RUNTYPE
      public :: MCL_SEND_CLIENT_APIVER
      public :: MCL_SEND_BOX_VECTORS
      public :: MCL_SEND_NUM_PARTICLES
      public :: MCL_SEND_NUM_PARTICLE_SPECIES
      public :: MCL_SEND_PARTICLE_SPECIES_IDS
      public :: MCL_SEND_SPECIES_ELEMENTS
      public :: MCL_SEND_SPECIES_MASSES
      public :: MCL_SEND_SPECIES_LABELS
      public :: MCL_SEND_NUM_FRAGMENTS
      public :: MCL_SEND_NUM_PARTICLES_IN_FRAGMENTS
      public :: MCL_SEND_PARTICLE_IDS_IN_FRAGMENTS
      public :: MCL_SEND_PARTICLE_MULTIPOLE_ORDER
      public :: MCL_SEND_PARTICLE_MULTIPOLES
      public :: MCL_SEND_NUM_CONSTR_BONDS
      public :: MCL_SEND_CONSTR_BOND_PARTICLE_IDS
      public :: MCL_SEND_CONSTR_BOND_LENGTHS
      public :: MCL_SEND_NUM_CONSTR_ANGLES
      public :: MCL_SEND_PARTICLE_POSITIONS
      public :: MCL_SEND_NUM_BOUNDARY_BONDS
      public :: MCL_SEND_BOUNDARY_BOND_ATOMS
      public :: MCL_RECV_PARTICLE_POSITIONS
      public :: MCL_COMPUTE_FORCES
      public :: MCL_SEND_PARTICLE_FORCES
      public :: MCL_SEND_ENERGY
#endif
c
      logical :: use_mimic = .false.
c     MCL program id of the MiMiC server (always 0)
      integer, parameter :: MIMIC_SERVER = 0
c     cached mapping of atoms to dense force-field-type "species"
      integer :: mimic_nspec
      logical :: mimic_species_built = .false.
      integer, allocatable :: mimic_spec_of_atom(:)
      integer, allocatable :: mimic_spec_type(:)
      save
c
      contains
c
c     split COMM_WORLD via MCL when running as a MiMiC MM client
c
      subroutine mimic_init(comm)
      implicit none
      integer, intent(inout) :: comm
#ifdef MIMIC
c     always join the MCL communicator in a MiMiC-enabled build; whether
c     we actually act as a client is decided later from the selected
c     integrator (use_mimic, set when the MIMIC integrator is chosen)
      call mcl_initialize(comm)
#endif
      end subroutine
c
c     receive the next request from the server and share it with all ranks
c
      subroutine mimic_get_request(request)
      use domdec, only: rank, COMM_TINKER
      use mpi
      implicit none
      integer :: request
#ifdef MIMIC
      integer :: ierr
      if (rank == 0) then
         call mcl_receive(request, 1, MCL_REQUEST, MIMIC_SERVER)
      end if
      call MPI_BCAST(request, 1, MPI_INT, 0, COMM_TINKER, ierr)
#else
      request = 0
      call mimic_disabled
#endif
      end subroutine
c
c     disconnect from the MiMiC server
c
      subroutine mimic_finalize
      implicit none
#ifdef MIMIC
      call mcl_finalize()
#endif
      end subroutine
c
c     abort the whole MiMiC job through MCL
c
      subroutine mimic_abort(code)
      implicit none
      integer :: code
#ifdef MIMIC
      logical :: initialized
      call mcl_is_initialized(initialized)
      if (initialized)  call mcl_abort(code)
#endif
      end subroutine
c
c     send our MCL program id
c
      subroutine mimic_send_client_id
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: id
      if (rank /= 0) return
      call mcl_get_program_id(id)
      call mcl_send(id, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the client program name
c
      subroutine mimic_send_client_name
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      character(len=*), parameter :: cname = 'TINKER-HP'
      integer :: ln
      if (rank /= 0) return
      ln = len(cname)
      call mcl_send(ln, 1, MCL_LENGTH, MIMIC_SERVER)
      call mcl_send(cname, ln, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the calculation tag
c
      subroutine mimic_send_client_runtype
      use domdec, only: rank
      use iounit, only: iout
      use potent, only: use_charge, use_mpole, use_polar
      implicit none
#ifdef MIMIC
      integer :: runtype
      if (use_mpole .and. use_polar) then
         runtype = MCL_RUNTYPE_MM_POL_AMOEBA
      else if (use_charge) then
         runtype = MCL_RUNTYPE_MM
      else
         if (rank == 0) write(iout,10)
   10    format (/,' MIMIC_SEND_CLIENT_RUNTYPE  --  unsupported force',
     &             ' field: only fixed-charge and polarizable AMOEBA',
     &           /,' (use_mpole .and. use_polar) are implemented')
         call fatal
      end if
      if (rank == 0) call mcl_send(runtype,1,MCL_DATA,MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send our MCL API version
c
      subroutine mimic_send_client_apiver
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: ver(3)
      if (rank /= 0) return
      call mcl_get_api_version(ver)
      call mcl_send(ver, 3, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the periodic box vectors (rows a,b,c) in bohr
c
      subroutine mimic_send_box
      use boxes, only: lvec
      use domdec, only: rank
      use units, only: bohr
      implicit none
#ifdef MIMIC
      real*8 :: buf(9)
      integer :: i, j, k
      if (rank /= 0) return
      k = 0
      do i = 1, 3
         do j = 1, 3
            k = k + 1
            buf(k) = lvec(i,j) / bohr
         end do
      end do
      call mcl_send(buf, 9, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of particles
c
      subroutine mimic_send_natoms
      use atoms, only: n
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      if (rank /= 0) return
      call mcl_send(n, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     build the atom -> dense-species map once (species = force-field type)
c
      subroutine mimic_build_species
      use atoms, only: n, type
      use sizes, only: maxtyp
      implicit none
#ifdef MIMIC
      integer :: i, t, ns
      integer :: tmap(maxtyp)
      if (mimic_species_built) return
      if (allocated(mimic_spec_of_atom)) deallocate(mimic_spec_of_atom)
      allocate(mimic_spec_of_atom(n))
      tmap = 0
      ns = 0
      do i = 1, n
         t = type(i)
         if (tmap(t) == 0) then
            ns = ns + 1
            tmap(t) = ns
         end if
         mimic_spec_of_atom(i) = tmap(t)
      end do
      mimic_nspec = ns
      if (allocated(mimic_spec_type)) deallocate(mimic_spec_type)
      allocate(mimic_spec_type(ns))
      do t = 1, maxtyp
         if (tmap(t) /= 0) mimic_spec_type(tmap(t)) = t
      end do
      mimic_species_built = .true.
#endif
      end subroutine
c
c     send the number of distinct species
c
      subroutine mimic_send_num_species
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      call mimic_build_species
      if (rank /= 0) return
      call mcl_send(mimic_nspec, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the per-atom species id
c
      subroutine mimic_send_species_ids
      use atoms, only: n
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      call mimic_build_species
      if (rank /= 0) return
      call mcl_send(mimic_spec_of_atom, n, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the element (atomic number) of each species
c
      subroutine mimic_send_species_elements
      use domdec, only: rank
      use katoms, only: atmnum
      implicit none
#ifdef MIMIC
      integer, allocatable :: buf(:)
      integer :: s
      call mimic_build_species
      if (rank /= 0) return
      allocate(buf(mimic_nspec))
      do s = 1, mimic_nspec
         buf(s) = atmnum(mimic_spec_type(s))
      end do
      call mcl_send(buf, mimic_nspec, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the mass (amu) of each species
c
      subroutine mimic_send_species_masses
      use domdec, only: rank
      use katoms, only: weight
      implicit none
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      integer :: s
      call mimic_build_species
      if (rank /= 0) return
      allocate(buf(mimic_nspec))
      do s = 1, mimic_nspec
         buf(s) = weight(mimic_spec_type(s))
      end do
      call mcl_send(buf, mimic_nspec, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the comma-delimited element symbols of the species
c
      subroutine mimic_send_species_labels
      use domdec, only: rank
      use katoms, only: symbol
      use sizes, only: maxtyp
      implicit none
#ifdef MIMIC
      character(len=4*maxtyp) :: labels
      character(len=3) :: sym
      integer :: s, p, ls
      call mimic_build_species
      if (rank /= 0) return
      p = 0
      do s = 1, mimic_nspec
         if (s > 1) then
            p = p + 1
            labels(p:p) = ','
         end if
         sym = symbol(mimic_spec_type(s))
         ls = len_trim(sym)
         if (ls == 0) ls = 1
         labels(p+1:p+ls) = sym(1:ls)
         p = p + ls
      end do
      call mcl_send(p, 1, MCL_LENGTH, MIMIC_SERVER)
      call mcl_send(labels, p, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of fragments (molecules)
c
      subroutine mimic_send_num_fragments
      use domdec, only: rank
      use molcul, only: nmol
      implicit none
#ifdef MIMIC
      if (rank /= 0) return
      call mcl_send(nmol, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of atoms in each fragment
c
      subroutine mimic_send_fragment_natoms
      use domdec, only: rank
      use molcul, only: nmol, imol
      implicit none
#ifdef MIMIC
      integer, allocatable :: buf(:)
      integer :: m
      if (rank /= 0) return
      allocate(buf(nmol))
      do m = 1, nmol
         buf(m) = imol(2,m) - imol(1,m) + 1
      end do
      call mcl_send(buf, nmol, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the atom ids grouped by fragment (molecule order, 1-based)
c
      subroutine mimic_send_fragment_ids
      use atoms, only: n
      use domdec, only: rank
      use molcul, only: kmol
      implicit none
#ifdef MIMIC
      if (rank /= 0) return
      call mcl_send(kmol, n, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the multipole order (0 = point charges)
c
      subroutine mimic_send_multipole_order
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: order
      if (rank /= 0) return
      order = 0
      call mcl_send(order, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the per-atom point charges (fixed-charge FF; mpole deferred)
c
      subroutine mimic_send_multipoles
      use atoms, only: n
      use charge, only: pchg, chglist
      use domdec, only: rank
      use potent, only: use_charge
      implicit none
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      integer :: i, k
      if (rank /= 0) return
      allocate(buf(n))
      buf = 0.0d0
      if (use_charge) then
         do i = 1, n
            k = chglist(i)
            if (k /= 0) buf(i) = pchg(k)
         end do
      end if
      call mcl_send(buf, n, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of constraint (rigid) bonds
c
      subroutine mimic_send_num_constr_bonds
      use domdec, only: rank
      use freeze, only: use_rattle, nrat
      implicit none
#ifdef MIMIC
      integer :: nc
      if (rank /= 0) return
      nc = 0
      if (use_rattle) nc = nrat
      call mcl_send(nc, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the atom pairs of each constraint bond (2*nrat, 1-based)
c
      subroutine mimic_send_constr_bond_ids
      use domdec, only: rank
      use freeze, only: use_rattle, nrat, irat
      implicit none
#ifdef MIMIC
      integer, allocatable :: buf(:)
      integer :: k, nc
      if (rank /= 0) return
      nc = 0
      if (use_rattle) nc = nrat
      allocate(buf(max(1,2*nc)))
      do k = 1, nc
         buf(2*k-1) = irat(1,k)
         buf(2*k)   = irat(2,k)
      end do
      call mcl_send(buf, 2*nc, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the target length of each constraint bond in bohr
c
      subroutine mimic_send_constr_bond_lengths
      use domdec, only: rank
      use freeze, only: use_rattle, nrat, krat
      use units, only: bohr
      implicit none
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      integer :: k, nc
      if (rank /= 0) return
      nc = 0
      if (use_rattle) nc = nrat
      allocate(buf(max(1,nc)))
      do k = 1, nc
         buf(k) = krat(k) / bohr
      end do
      call mcl_send(buf, nc, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of constraint angles (none)
c
      subroutine mimic_send_num_constr_angles
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: na
      if (rank /= 0) return
      na = 0
      call mcl_send(na, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the current particle positions in bohr
c
      subroutine mimic_send_positions
      use atoms, only: n, x, y, z
      use domdec, only: rank
      use units, only: bohr
      implicit none
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      integer :: i
      if (rank /= 0) return
      allocate(buf(3*n))
      do i = 1, n
         buf(3*i-2) = x(i) / bohr
         buf(3*i-1) = y(i) / bohr
         buf(3*i)   = z(i) / bohr
      end do
      call mcl_send(buf, 3*n, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the number of QM/MM boundary bonds (exclusions deferred)
c
      subroutine mimic_send_num_boundary_bonds
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: nb
      if (rank /= 0) return
      nb = 0
      call mcl_send(nb, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the QM/MM boundary bond atoms (empty; exclusions deferred)
c
      subroutine mimic_send_boundary_bonds
      use domdec, only: rank
      implicit none
#ifdef MIMIC
      integer :: buf(1)
      if (rank /= 0) return
      call mcl_send(buf, 0, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     receive updated positions from the server and install them on all
c     ranks (bohr -> Angstrom, then broadcast)
c
      subroutine mimic_recv_positions
      use atoms, only: n, x, y, z
      use domdec, only: rank, COMM_TINKER
      use units, only: bohr
      use mpi
      implicit none
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      integer :: i, ierr
      if (rank == 0) then
         allocate(buf(3*n))
         call mcl_receive(buf, 3*n, MCL_DATA, MIMIC_SERVER)
         do i = 1, n
            x(i) = buf(3*i-2) * bohr
            y(i) = buf(3*i-1) * bohr
            z(i) = buf(3*i)   * bohr
         end do
         deallocate(buf)
      end if
      call MPI_BCAST(x, n, MPI_REAL8, 0, COMM_TINKER, ierr)
      call MPI_BCAST(y, n, MPI_REAL8, 0, COMM_TINKER, ierr)
      call MPI_BCAST(z, n, MPI_REAL8, 0, COMM_TINKER, ierr)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the forces (Hartree/bohr) gathered on rank 0; the argument is
c     the global gradient dE/dx (kcal/mol/Angstrom), so force = -grad
c
      subroutine mimic_send_forces(grad)
      use atoms, only: n
      use domdec, only: rank
      use units, only: bohr, hartree
      implicit none
      real*8 :: grad(3,*)
#ifdef MIMIC
      real*8, allocatable :: buf(:)
      real*8 :: fac
      integer :: i
      if (rank /= 0) return
      fac = bohr / hartree
      allocate(buf(3*n))
      do i = 1, n
         buf(3*i-2) = -grad(1,i) * fac
         buf(3*i-1) = -grad(2,i) * fac
         buf(3*i)   = -grad(3,i) * fac
      end do
      call mcl_send(buf, 3*n, MCL_DATA, MIMIC_SERVER)
      deallocate(buf)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the potential energy in Hartree
c
      subroutine mimic_send_energy(energy)
      use domdec, only: rank
      use units, only: hartree
      implicit none
      real*8 :: energy
#ifdef MIMIC
      real*8 :: e
      if (rank /= 0) return
      e = energy / hartree
      call mcl_send(e, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
#ifndef MIMIC
c
c     reached only if a MiMiC routine is called in a non-MiMiC build
c
      subroutine mimic_disabled
      use iounit, only: iout
      implicit none
      write(iout,10)
   10 format (/,' MIMIC  --  this binary was not built with MiMiC',
     &          ' support; reconfigure with --enable-mimic')
      call fatal
      end subroutine
#endif
      end module
