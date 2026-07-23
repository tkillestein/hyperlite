!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
program test_philox
!***********************************************************************
! Unit tests for randommod (Philox-4x32-10, OVERHAUL Phase 5b):
! - Random123 known-answer vectors for the bare philox4x32 block
! - draws lie strictly in (0,1)
! - per-particle keying is deterministic and independent of prior state
! - uniformity sanity check on the sample mean
!***********************************************************************
  use kindmod
  use randommod
  implicit none
!
  integer :: nfail = 0
!
  call kat_vectors
  call range_and_determinism
  call uniformity
!
  if(nfail > 0) then
    write(6,*) 'test_philox: FAILED, nfail=', nfail
    stop 1
  endif
  write(6,*) 'test_philox: all tests passed'
!
contains
!
  subroutine check(cond,label)
    logical,intent(in) :: cond
    character(*),intent(in) :: label
    if(.not.cond) then
      nfail = nfail + 1
      write(6,*) 'FAIL: ', label
    endif
  end subroutine check
!
  subroutine kat_vectors
!-- Random123 v1.14.0 known-answer vectors (tests/kat_vectors), philox4x32 r=10
    integer(i8) :: ctr(4),key(2),out(4),expct(4)
!-- all-zero counter and key
    ctr = 0
    key = 0
    expct = [int(z'6627e8d5',i8),int(z'e169c58d',i8), &
             int(z'bc57ac4c',i8),int(z'9b00dbd8',i8)]
    call philox4x32(ctr,key,out)
    call check(all(out==expct),'KAT zeros')
!-- all-ones counter and key
    ctr = int(z'ffffffff',i8)
    key = int(z'ffffffff',i8)
    expct = [int(z'408f276d',i8),int(z'41c83b0e',i8), &
             int(z'a20bc7c6',i8),int(z'6d5451fd',i8)]
    call philox4x32(ctr,key,out)
    call check(all(out==expct),'KAT ones')
!-- digits of pi
    ctr = [int(z'243f6a88',i8),int(z'85a308d3',i8), &
           int(z'13198a2e',i8),int(z'03707344',i8)]
    key = [int(z'a4093822',i8),int(z'299f31d0',i8)]
    expct = [int(z'd16cfe09',i8),int(z'94fdcceb',i8), &
             int(z'5001e420',i8),int(z'24126ea1',i8)]
    call philox4x32(ctr,key,out)
    call check(all(out==expct),'KAT pi')
  end subroutine kat_vectors
!
  subroutine range_and_determinism
    type(rnd_t) :: s1,s2
    real(dp) :: x,x1(20),x2(20)
    integer(i4) :: i4v
    integer :: i
    logical :: lok
    call rnd_init(2,0)
!-- draws strictly in (0,1)
    lok = .true.
    do i=1,10000
      call rnd_r(x,rnd_state)
      if(x<=0d0 .or. x>=1d0) lok = .false.
    enddo
    call check(lok,'rnd_r in (0,1)')
!-- rnd_i in [0,huge]
    lok = .true.
    do i=1,10000
      call rnd_i(i4v,rnd_state)
      if(i4v<0) lok = .false.
    enddo
    call check(lok,'rnd_i nonnegative')
!-- same (key0,epoch,ipart) gives the same sequence, whatever came before
    call rnd_seed_particle(s1,42)
    do i=1,20
      call rnd_r(x1(i),s1)
    enddo
    call rnd_seed_particle(s2,7)  !disturb the state
    call rnd_r(x,s2)
    call rnd_seed_particle(s2,42)
    do i=1,20
      call rnd_r(x2(i),s2)
    enddo
    call check(all(x1==x2),'per-particle keying deterministic')
!-- different particles give different sequences
    call rnd_seed_particle(s2,43)
    call rnd_r(x,s2)
    call check(x/=x1(1),'distinct particles distinct streams')
!-- epoch advance changes the sequence
    call rnd_advance_epoch
    call rnd_seed_particle(s2,42)
    call rnd_r(x,s2)
    call check(x/=x1(1),'epoch advance changes stream')
  end subroutine range_and_determinism
!
  subroutine uniformity
!-- sample mean of n uniform draws must be 0.5 within 4 sigma
    integer,parameter :: n = 100000
    real(dp) :: x,s,tol
    integer :: i
    type(rnd_t) :: st
    call rnd_seed_particle(st,1)
    s = 0d0
    do i=1,n
      call rnd_r(x,st)
      s = s + x
    enddo
    s = s/n
    tol = 4d0/sqrt(12d0*n)  !4 sigma of the mean
    call check(abs(s-0.5d0)<tol,'sample mean uniform')
  end subroutine uniformity
!
end program test_philox
! vim: fdm=marker
