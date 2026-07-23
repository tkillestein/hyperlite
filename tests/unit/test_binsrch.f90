!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
program test_binsrch
!***********************************************************************
! Unit tests for binsrch: compare against a straightforward linear scan
! on a log-spaced grid (the wavelength-grid use case), plus edge and
! out-of-range behavior in both widerange modes.
!***********************************************************************
  use kindmod
  use randommod
  implicit none
!
  interface
    pure function binsrch(x,arr,n,widerange)
      use kindmod
      integer :: binsrch
      integer,intent(in) :: n
      real(dp),intent(in) :: x
      real(dp),intent(in) :: arr(n)
      logical,intent(in) :: widerange
    end function binsrch
  end interface
!
  integer,parameter :: n = 101
  real(dp) :: arr(n),x
  integer :: i,nfail
  type(rnd_t) :: st
!
  nfail = 0
!-- log-spaced grid, like the wavelength groups
  do i=1,n
    arr(i) = 1d3*10d0**(1.3d0*(i-1)/(n-1))
  enddo
!
!-- random interior points against the linear reference
  call rnd_init(1,0)
  call rnd_seed_particle(st,1)
  do i=1,10000
    call rnd_r(x,st)
    x = arr(1) + x*(arr(n)-arr(1))
    call check(binsrch(x,arr,n,.false.)==linref(x),'interior random')
    call check(binsrch(x,arr,n,.true.)==linref(x),'interior random wide')
  enddo
!
!-- every left edge maps to its own interval; right end to n-1
  do i=1,n-1
    call check(binsrch(arr(i),arr,n,.false.)==i,'left edge')
  enddo
  call check(binsrch(arr(n),arr,n,.false.)==n-1,'right endpoint clamped')
!
!-- out-of-range: clamped without widerange, 0/n with it
  call check(binsrch(arr(1)*0.5d0,arr,n,.false.)==1,'below clamped')
  call check(binsrch(arr(n)*2d0,arr,n,.false.)==n-1,'above clamped')
  call check(binsrch(arr(1)*0.5d0,arr,n,.true.)==0,'below wide')
  call check(binsrch(arr(n)*2d0,arr,n,.true.)==n,'above wide')
!
!-- n=2 quick-return branch
  call check(binsrch(1.5d0,[1d0,2d0],2,.false.)==1,'n=2 inside')
  call check(binsrch(0.5d0,[1d0,2d0],2,.true.)==0,'n=2 below wide')
  call check(binsrch(2.5d0,[1d0,2d0],2,.true.)==2,'n=2 above wide')
  call check(binsrch(0.5d0,[1d0,2d0],2,.false.)==1,'n=2 below clamped')
!
  if(nfail > 0) then
    write(6,*) 'test_binsrch: FAILED, nfail=', nfail
    stop 1
  endif
  write(6,*) 'test_binsrch: all tests passed'
!
contains
!
  subroutine check(cond,label)
    logical,intent(in) :: cond
    character(*),intent(in) :: label
    if(.not.cond) then
      nfail = nfail + 1
      write(6,*) 'FAIL: ', label, x
    endif
  end subroutine check
!
  pure function linref(xx) result(j)
!-- linear-scan reference: interval j with arr(j) <= xx < arr(j+1)
    real(dp),intent(in) :: xx
    integer :: j
    do j=1,n-2
      if(xx < arr(j+1)) return
    enddo
    j = n-1
  end function linref
!
end program test_binsrch
! vim: fdm=marker
