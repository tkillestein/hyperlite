!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module particlemod
  use kindmod

  implicit none

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
  type(packet),allocatable,target :: prt_particles(:)  !(prt_npartmax)
  logical,allocatable :: prt_isvacant(:)  !(prt_npartmax)
!
  integer :: prt_npartmax

  save

  contains

  subroutine particle_alloc(ltalk)
!--------------------------------------------------
    implicit none
    logical,intent(in) :: ltalk

    integer :: n

!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax),prt_isvacant(prt_npartmax))
    prt_isvacant = .true.
!
!-- print size only on master rank
    if(ltalk) then
      n = int((storage_size(prt_particles,kind=i8)/8_i8)* &
        size(prt_particles,kind=i8)/1024_i8) !kB
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


end module particlemod
! vim: fdm=marker
