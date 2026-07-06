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
c     are stubs. For bring-up the client is activated at runtime by
c     setting the environment variable TINKER_MIMIC.
c
      module mimic
#ifdef MIMIC
      use mcl
#endif
      use domdec, only: rank
      use iounit, only: iout
      use potent, only: use_charge, use_mpole, use_polar
      implicit none
      private
      public :: use_mimic
      public :: mimic_init
      public :: mimic_get_request
      public :: mimic_send_client_id
      public :: mimic_send_client_name
      public :: mimic_send_client_runtype
      public :: mimic_send_client_apiver
      public :: mimic_finalize
      public :: mimic_abort
c     re-exported so the client loop only needs use mimic
      public :: rank
      public :: iout
#ifdef MIMIC
      public :: MCL_EXIT
      public :: MCL_SEND_CLIENT_ID
      public :: MCL_SEND_CLIENT_NAME
      public :: MCL_SEND_CLIENT_RUNTYPE
      public :: MCL_SEND_CLIENT_APIVER
#endif
c
      logical :: use_mimic = .false.
c     MCL program id of the MiMiC server (always 0)
      integer, parameter :: MIMIC_SERVER = 0
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
      character(len=32) envval
      integer ln, st
      call get_environment_variable('TINKER_MIMIC',envval,ln,st)
      if (st.eq.0 .and. ln.gt.0 .and. envval(1:1).ne.'0') then
         use_mimic = .true.
         call mcl_initialize(comm)
      end if
#endif
      end subroutine
c
c     receive the next request from the server
c
      subroutine mimic_get_request(request)
      implicit none
      integer request
#ifdef MIMIC
      call mcl_receive(request, 1, MCL_REQUEST, MIMIC_SERVER)
#else
      request = 0
      call mimic_disabled
#endif
      end subroutine
c
c     send our MCL program id
c
      subroutine mimic_send_client_id
      implicit none
#ifdef MIMIC
      integer id
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
      implicit none
#ifdef MIMIC
      character(len=*), parameter :: name = 'TINKER-HP'
      integer ln
      ln = len(name)
      call mcl_send(ln, 1, MCL_LENGTH, MIMIC_SERVER)
      call mcl_send(name, ln, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send the calculation tag
c
      subroutine mimic_send_client_runtype
      implicit none
#ifdef MIMIC
      integer runtype
      if (use_mpole .and. use_polar) then
         runtype = MCL_RUNTYPE_MM_POL_AMOEBA
      else if (use_charge) then
         runtype = MCL_RUNTYPE_MM
      else
         if (rank.eq.0) write(iout,10)
   10    format (/,' MIMIC_SEND_CLIENT_RUNTYPE  --  unsupported force',
     &             ' field: only fixed-charge and polarizable AMOEBA',
     &           /,' (use_mpole .and. use_polar) are implemented')
         call fatal
      end if
      call mcl_send(runtype, 1, MCL_DATA, MIMIC_SERVER)
#else
      call mimic_disabled
#endif
      end subroutine
c
c     send our MCL API version
c
      subroutine mimic_send_client_apiver
      implicit none
#ifdef MIMIC
      integer ver(3)
      call mcl_get_api_version(ver)
      call mcl_send(ver, 3, MCL_DATA, MIMIC_SERVER)
#else
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
      integer code
#ifdef MIMIC
      call mcl_abort(code)
#endif
      end subroutine
c
#ifndef MIMIC
c
c     reached only if a MiMiC routine is called in a non-MiMiC build
c
      subroutine mimic_disabled
      implicit none
      write(iout,10)
   10 format (/,' MIMIC  --  this binary was not built with MiMiC',
     &          ' support; reconfigure with --enable-mimic')
      call fatal
      end subroutine
#endif
      end module
