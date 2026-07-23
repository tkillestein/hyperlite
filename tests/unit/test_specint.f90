!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
program test_specint
!***********************************************************************
! Unit tests for specint (composite-Simpson integral of x^m/(e^x-1)):
! - summed over a group grid it must reproduce the analytic Planck
!   integral  int_0^inf x^3/(e^x-1) dx = pi^4/15
! - m=2 analog: int_0^inf x^2/(e^x-1) dx = 2 zeta(3)
! - additivity over a split interval
!***********************************************************************
  use kindmod
  implicit none
!
  interface
    elemental function specint(x1,x2,m,nn) result(ss)
      use kindmod
      real(dp) :: ss
      integer,intent(in) :: m
      real(dp),intent(in) :: x1,x2
      integer,intent(in),optional :: nn
    end function specint
  end interface
!
  real(dp),parameter :: pi = 3.14159265358979324d0
  real(dp),parameter :: planck3 = pi**4/15d0        !int x^3/(e^x-1)
  real(dp),parameter :: planck2 = 2.40411380631919d0 !2*zeta(3)
  integer,parameter :: nseg = 200
  real(dp) :: x1,x2,s,a,b,c
  integer :: i,nfail
!
  nfail = 0
!
!-- m=3: sum over log-spaced segments spanning (1e-4, 60)
  s = 0d0
  do i=1,nseg
    x1 = 1d-4*(6d5)**(dble(i-1)/nseg)
    x2 = 1d-4*(6d5)**(dble(i)/nseg)
    s = s + specint(x1,x2,3)
  enddo
  call check(abs(s/planck3-1d0)<1d-6,'Planck integral m=3',s)
!
!-- m=2
  s = 0d0
  do i=1,nseg
    x1 = 1d-4*(6d5)**(dble(i-1)/nseg)
    x2 = 1d-4*(6d5)**(dble(i)/nseg)
    s = s + specint(x1,x2,2)
  enddo
  call check(abs(s/planck2-1d0)<1d-6,'Planck integral m=2',s)
!
!-- additivity: [a,c] = [a,b] + [b,c] at fine n
  a = 0.5d0
  b = 2.7d0
  c = 9d0
  s = specint(a,b,3,1000) + specint(b,c,3,1000)
  call check(abs(s/specint(a,c,3,2000)-1d0)<1d-9,'additivity',s)
!
!-- x2=0 quick return
  call check(specint(1d0,0d0,3)==0d0,'x2=0 returns 0',0d0)
!
  if(nfail > 0) then
    write(6,*) 'test_specint: FAILED, nfail=', nfail
    stop 1
  endif
  write(6,*) 'test_specint: all tests passed'
!
contains
!
  subroutine check(cond,label,val)
    logical,intent(in) :: cond
    character(*),intent(in) :: label
    real(dp),intent(in) :: val
    if(.not.cond) then
      nfail = nfail + 1
      write(6,*) 'FAIL: ', label, val
    endif
  end subroutine check
!
end program test_specint
! vim: fdm=marker
