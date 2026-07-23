!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module particlemod
  use kindmod

  implicit none

!-- working struct for one particle in the transport kernel; storage is
!-- struct-of-arrays (below) for SIMD-friendly access (OVERHAUL Phase 5)
  type packet
     real(dp) :: x, y, z
     real(dp) :: mu, om
     real(dp) :: e, e0, wl
  end type packet
!
!-- secondary particle properties
  type packet2
     real(dp) :: dist          !particle travel distance
     character(4) :: stat    !particle status: live, cens, flux, dead
     integer :: ix, iy, iz   !positional cell indices
     integer :: ic, ig       !index into compressed domain arrays, group index
     integer :: itype        !IMC or DDMC type
     integer :: ipart, istep !particle number and transport step number
     integer :: idist        !transport distance identifier
  end type packet2
!
!-- particle storage: struct of arrays (prt_npartmax)
  real(dp),allocatable :: prt_x(:)
  real(dp),allocatable :: prt_mu(:)
  real(dp),allocatable :: prt_om(:)
  real(dp),allocatable :: prt_e(:)
  real(dp),allocatable :: prt_e0(:)
  real(dp),allocatable :: prt_wl(:)
  logical,allocatable :: prt_isvacant(:)  !(prt_npartmax)
!
!-- transverse coordinates are constant in 1-D spherical geometry: not
!-- stored per particle; set once at allocation from the grid edges
  real(dp) :: prt_y0 = 0d0
  real(dp) :: prt_z0 = 0d0
!
  integer :: prt_npartmax

  save

  contains

  subroutine particle_alloc(ltalk,y0,z0)
!--------------------------------------------------
    implicit none
    logical,intent(in) :: ltalk
    real(dp),intent(in) :: y0,z0

    integer :: n

!-- transverse coordinates (1-D)
    prt_y0 = y0
    prt_z0 = z0

!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_x(prt_npartmax),prt_mu(prt_npartmax), &
      prt_om(prt_npartmax),prt_e(prt_npartmax),prt_e0(prt_npartmax), &
      prt_wl(prt_npartmax),prt_isvacant(prt_npartmax))
    prt_isvacant = .true.
!
!-- print size only on master rank
    if(ltalk) then
      n = int(6_i8*(storage_size(prt_x,kind=i8)/8_i8)* &
        size(prt_x,kind=i8)/1024_i8) !kB
      write(6,*) 'ALLOC particles:',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif
!
!-- output
    if(ltalk) then
       write(6,*)
       write(6,*) 'particle array:'
       write(6,*) '===================='
       write(6,*) 'npart :',prt_npartmax,nint(prt_npartmax/1000d0),'k'
       write(6,*)
    endif

  end subroutine particle_alloc


  pure subroutine prt_gather(i,ptcl)
!-- load one particle from the SoA storage into a working packet
    integer,intent(in) :: i
    type(packet),intent(out) :: ptcl
    ptcl%x = prt_x(i)
    ptcl%y = prt_y0
    ptcl%z = prt_z0
    ptcl%mu = prt_mu(i)
    ptcl%om = prt_om(i)
    ptcl%e = prt_e(i)
    ptcl%e0 = prt_e0(i)
    ptcl%wl = prt_wl(i)
  end subroutine prt_gather


  subroutine prt_scatter(i,ptcl)
!-- store a working packet back into the SoA storage (y,z not stored);
!-- not pure: assigns module storage
    integer,intent(in) :: i
    type(packet),intent(in) :: ptcl
    prt_x(i) = ptcl%x
    prt_mu(i) = ptcl%mu
    prt_om(i) = ptcl%om
    prt_e(i) = ptcl%e
    prt_e0(i) = ptcl%e0
    prt_wl(i) = ptcl%wl
  end subroutine prt_scatter

end module particlemod
! vim: fdm=marker
