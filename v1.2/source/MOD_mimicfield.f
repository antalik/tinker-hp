c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c     Ecole Polytechnique Federale de Lausanne
c
c     ###############################################################
c     ##                                                           ##
c     ##  module mimicfield  --  external electric field on atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c     Storage for the electric field the MiMiC server evaluates at the
c     atomic sites (e.g. the field of the QM density in a polarizable
c     QM/MM run). The field enters the AMOEBA induced-dipole equations
c     as a constant contribution to the right hand side, added to both
c     the d and p fields.
c
c     This module deliberately does not depend on the MCL library: it
c     is used by the induced-dipole solvers, which must still compile
c     in a build configured without MiMiC support.
c
c     mimic_efield     external field at each atom, in Tinker units of
c                      e/Angstrom**2, replicated on every rank
c     mimic_field_set  .true. once a field has been received
c
      module mimicfield
      implicit none
      real*8, allocatable :: mimic_efield(:,:)
      logical :: mimic_field_set = .false.
      save
c
      contains
c
c     allocate the field buffer on first use and zero it
c
      subroutine mimic_alloc_field
      use atoms, only: n
      implicit none
      if (.not. allocated(mimic_efield)) then
         allocate(mimic_efield(3,n))
         mimic_efield = 0.0d0
      end if
      end subroutine
c
c     add the external field to the right hand side of the induced
c     dipole equations; ef is indexed by the local pole index and holds
c     nrhs (= 2) fields, the d field and the p field
c
      subroutine mimic_add_field(nrhs, npl, ef)
      use atmlst, only: poleglob
      use mpole, only: ipole
      implicit none
      integer :: nrhs, npl
      real*8 :: ef(3,nrhs,*)
      integer :: i, j, k, iipole, iglob
c
      if (.not. mimic_field_set)  return
      if (.not. allocated(mimic_efield))  return
c
      do i = 1, npl
         iipole = poleglob(i)
         iglob = ipole(iipole)
         do k = 1, nrhs
            do j = 1, 3
               ef(j,k,i) = ef(j,k,i) + mimic_efield(j,iglob)
            end do
         end do
      end do
      end subroutine
c
      end module
