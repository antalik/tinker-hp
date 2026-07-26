c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mimicgroup  --  MiMiC atom subset and MM exclusions  ##
c     ##                                                              ##
c     ##################################################################
c
c     The MIMIC_GROUP keyword marks a subset of atoms handed to the MiMiC
c     server (e.g. the QM region in a QM/MM run). For these atoms, the
c     following exclusions are made:
c       - charges are zeroed (kcharge)
c       - permanent multipoles are zeroed (kmpole)
c       - polarizabilities are zeroed (kpolar), so no induced dipoles
c       - group-group vdW pairs are skipped (ehal/elj), group-MM kept
c       - bonded terms dropped iff ALL of its atoms are in the group
c
c     The multipole sites themselves are kept in the pole list with their
c     local frames intact, so that atom and pole indices are unchanged.
c
c     nmimic       number of MiMiC-group atoms
c     mimic_group  per-atom logical flag, .true. for MiMiC-group atoms
c
      module mimicgroup
      implicit none
      integer :: nmimic = 0
      logical, allocatable :: mimic_group(:)
      save
      private :: all_in_mimic
c
      contains
c
c     read the MIMIC_GROUP keyword(s) into the mimic_group mask, expanding
c     the negative-number ranges as the ACTIVE/INACTIVE keywords do
c
      subroutine mimic_read_group
      use atoms, only: n
      use keys, only: nkey, keyline
      implicit none
      integer :: i, k, j, next
      integer, allocatable :: lst(:)
      character(len=20) :: keyword
      character(len=240) :: record, string
c
      if (.not. allocated(mimic_group))  allocate(mimic_group(n))
      mimic_group = .false.
      allocate(lst(n+1))
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext(record, keyword, next)
         call upcase(keyword)
         string = record(next:240)
         if (keyword(1:12) == 'MIMIC_GROUP ' .or.
     &       keyword(1:12) == 'MIMIC-GROUP ') then
            do i = 1, n+1
               lst(i) = 0
            end do
            read(string,*,err=10,end=10) (lst(i),i=1,n)
   10       continue
            i = 1
            do while (lst(i) /= 0)
               if (lst(i) > 0) then
                  k = lst(i)
                  if (k >= 1 .and. k <= n)  mimic_group(k) = .true.
                  i = i + 1
               else
                  do k = abs(lst(i)), abs(lst(i+1))
                     if (k >= 1 .and. k <= n)  mimic_group(k) = .true.
                  end do
                  i = i + 2
               end if
            end do
         end if
      end do
c
      deallocate(lst)
      nmimic = count(mimic_group)
      end subroutine
c
c     .true. if all atoms in the list are in the MiMiC group (valid indices)
c
      logical function all_in_mimic(list)
      implicit none
      integer :: list(:)
      integer :: j
      all_in_mimic = .true.
      do j = 1, size(list)
         if (list(j) < 1) then
            all_in_mimic = .false.
            return
         end if
         if (.not. mimic_group(list(j))) then
            all_in_mimic = .false.
            return
         end if
      end do
      end function
c
c     zero the induced dipoles of the MiMiC group after a solve
c
c     Setting polarity to zero makes the initial guess mu = polarity*ef
c     vanish, but it does not make the converged dipole vanish exactly:
c     the solvers floor the diagonal of the polarizability at tinypol
c     (1e-5), so a zero-polarity site relaxes to mu ~ tinypol*ef instead
c     of zero. That is negligible in the energy but would put nonsense in
c     the induced dipoles handed to the MiMiC server, so clear them here.
c
      subroutine mimic_zero_induced
      use atoms, only: n
      use mpole, only: pollist
      use polar, only: uind, uinp
      implicit none
      integer :: i, ii, j
c
      if (nmimic == 0)  return
      if (.not. allocated(mimic_group))  return
c
      do i = 1, n
         if (.not. mimic_group(i))  cycle
         ii = pollist(i)
         if (ii == 0)  cycle
         do j = 1, 3
            uind(j,ii) = 0.0d0
            uinp(j,ii) = 0.0d0
         end do
      end do
      end subroutine
c
c     drop every bonded term whose atoms are all in the MiMiC group by
c     zeroing its (globally indexed, shared-memory) force constants; terms
c     that cross the boundary keep at least one non-group atom and are left
c     intact. tortor (grid-based) is not handled here yet.
c
      subroutine mimic_exclude_bonded
      use domdec, only: hostrank, hostcomm
      use bond, only: nbond, ibnd, bk
      use angle, only: nangle, iang, ak
      use tors, only: ntors, itors,
     &                tors1, tors2, tors3, tors4, tors5, tors6
      use urey, only: nurey, iury, uk
      use strbnd, only: nstrbnd, isb, sbk
      use opbend, only: nopbend, iopb, opbk
      use opdist, only: nopdist, iopd, opdk
      use imptor, only: nitors, iitors, itors1, itors2, itors3
      use pitors, only: npitors, ipit, kpit
      use strtor, only: nstrtor, ist, kst
      use angtor, only: nangtor, iat, kant
      use angang, only: nangang, iaa, kaa
      use mpi
      implicit none
      integer :: i, j, k, ierr
c
      if (nmimic == 0)  return
c
c     shared-memory arrays are written on the node master only
c
      if (hostrank == 0) then
c
c        bond stretch (2 atoms)
c
         do i = 1, nbond
            if (all_in_mimic(ibnd(1:2,i)))  bk(i) = 0.0d0
         end do
c
c        angle bend (3 atoms)
c
         do i = 1, nangle
            if (all_in_mimic(iang(1:3,i)))  ak(i) = 0.0d0
         end do
c
c        Urey-Bradley (1-3 pair; test the 3 angle atoms)
c
         do i = 1, nurey
            if (all_in_mimic(iury(1:3,i)))  uk(i) = 0.0d0
         end do
c
c        torsion (4 atoms)
c
         do i = 1, ntors
            if (all_in_mimic(itors(1:4,i))) then
               tors1(:,i) = 0.0d0
               tors2(:,i) = 0.0d0
               tors3(:,i) = 0.0d0
               tors4(:,i) = 0.0d0
               tors5(:,i) = 0.0d0
               tors6(:,i) = 0.0d0
            end if
         end do
c
c        improper torsion (4 atoms)
c
         do i = 1, nitors
            if (all_in_mimic(iitors(1:4,i))) then
               itors1(:,i) = 0.0d0
               itors2(:,i) = 0.0d0
               itors3(:,i) = 0.0d0
            end if
         end do
c
c        out-of-plane distance (4 atoms)
c
         do i = 1, nopdist
            if (all_in_mimic(iopd(1:4,i)))  opdk(i) = 0.0d0
         end do
c
c        pi-orbital torsion (6 atoms)
c
         do i = 1, npitors
            if (all_in_mimic(ipit(1:6,i)))  kpit(i) = 0.0d0
         end do
c
c        stretch-bend (references an angle: isb(1,.))
c
         do i = 1, nstrbnd
            j = isb(1,i)
            if (all_in_mimic(iang(1:3,j)))  sbk(:,i) = 0.0d0
         end do
c
c        out-of-plane bend (references an angle: iopb(.))
c
         do i = 1, nopbend
            j = iopb(i)
            if (all_in_mimic(iang(1:3,j)))  opbk(i) = 0.0d0
         end do
c
c        stretch-torsion (references a torsion: ist(1,.))
c
         do i = 1, nstrtor
            j = ist(1,i)
            if (all_in_mimic(itors(1:4,j)))  kst(:,i) = 0.0d0
         end do
c
c        angle-torsion (references a torsion: iat(1,.))
c
         do i = 1, nangtor
            j = iat(1,i)
            if (all_in_mimic(itors(1:4,j)))  kant(:,i) = 0.0d0
         end do
c
c        angle-angle (couples two angles: iaa(1,.), iaa(2,.))
c
         do i = 1, nangang
            j = iaa(1,i)
            k = iaa(2,i)
            if (all_in_mimic(iang(1:3,j)) .and.
     &          all_in_mimic(iang(1:3,k)))  kaa(i) = 0.0d0
         end do
      end if
c
      call MPI_BARRIER(hostcomm, ierr)
      end subroutine
c
      end module
