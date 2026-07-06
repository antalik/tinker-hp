c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     drive the MiMiC MM client
c
      subroutine mimic_loop
      use mimic
      implicit none
#ifdef MIMIC
      integer request
      logical done
c
      if (rank.eq.0) write(iout,10)
   10 format (/,' Entering MiMiC QM/MM client loop')
c
      done = .false.
      do while (.not.done)
         call mimic_get_request(request)
         if (request.eq.MCL_SEND_CLIENT_ID) then
            call mimic_send_client_id
         else if (request.eq.MCL_SEND_CLIENT_NAME) then
            call mimic_send_client_name
         else if (request.eq.MCL_SEND_CLIENT_RUNTYPE) then
            call mimic_send_client_runtype
         else if (request.eq.MCL_SEND_CLIENT_APIVER) then
            call mimic_send_client_apiver
         else if (request.eq.MCL_EXIT) then
            done = .true.
         else
            if (rank.eq.0) write(iout,20) request
   20       format (/,' MIMIC_LOOP  --  unhandled request',i8,
     &                ', leaving client loop')
            done = .true.
         end if
      end do
c
      call mimic_finalize
      if (rank.eq.0) write(iout,30)
   30 format (/,' Leaving MiMiC QM/MM client loop')
#else
      write(iout,40)
   40 format (/,' MIMIC_LOOP  --  this binary was not built with',
     &          ' MiMiC support; reconfigure with --enable-mimic')
      call fatal
#endif
      end
