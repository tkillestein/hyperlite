!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module randommod
  use kindmod
  implicit none
!***********************************************************************
! Counter-based random number generator: Philox-4x32-10 (Salmon et al.
! 2011, "Parallel random numbers: as easy as 1, 2, 3").
!
! Replaces the mzran recurrence (OVERHAUL Phase 5).  Draws are a pure
! function of (key, counter):
!   key     = (key0, domain)   key0 = impi + nmpi*in_rnd_seed
!   counter = (idraw, epoch, ipart, 0)
! Transport draws are keyed per particle (rnd_seed_particle): epoch is
! the iteration number and ipart the particle slot, so the draw sequence
! is independent of OpenMP thread count and loop scheduling.  Service
! draws (source sampling) use the stateful stream set up by rnd_init
! (domain 0); per-thread scratch streams use domain=ithread.
!
! All 32-bit words are stored in i8 with the high bits zero; the 32x32
! multiply uses a 16-bit limb decomposition to stay portable (no
! unsigned, no int128).
!***********************************************************************
  integer(i4) :: rnd_imax = 2147483579 !kept for API compatibility
!
!-- state: counter + key + buffered outputs of the current block
  type :: rnd_t
    integer(i8) :: ctr(4)
    integer(i8) :: key(2)
    integer(i8) :: buf(4)
    integer :: ibuf !next unused buffer slot; 5 = buffer empty
  end type rnd_t
!
  integer :: rnd_nstate
  type(rnd_t) :: rnd_state              !service stream (domain 0)
  type(rnd_t),allocatable :: rnd_states(:) !per-thread streams
!
!-- key base (rank/seed offset) and transport epoch (iteration counter)
  integer(i8) :: rnd_key0 = 0
  integer(i8) :: rnd_epoch = 0
!
!-- transport domain tag (distinct from service 0 and thread ids >= 1)
  integer(i8),parameter,private :: rnd_domain_transport = 2_i8**31
!
  integer(i8),parameter,private :: mask32 = 2_i8**32 - 1
  integer(i8),parameter,private :: mask16 = 2_i8**16 - 1
!-- Philox multipliers and Weyl key increments (as positive i8)
  integer(i8),parameter,private :: pm0 = int(z'D2511F53',i8)
  integer(i8),parameter,private :: pm1 = int(z'CD9E8D57',i8)
  integer(i8),parameter,private :: pw0 = int(z'9E3779B9',i8)
  integer(i8),parameter,private :: pw1 = int(z'BB67AE85',i8)
!
  save
!
  contains
!
!
  subroutine rnd_init(n,ioffset)
!-- -----------------------------
    implicit none
    integer,intent(in) :: n,ioffset
!***********************************************************************
! Initialize n per-thread states plus the service stream.  ioffset
! (rank + nmpi*seed) selects the key so every rank/ensemble member has
! an independent stream family.
!***********************************************************************
    integer :: i
!-- alloc
    rnd_nstate = n
    allocate(rnd_states(rnd_nstate))
!-- key base
    rnd_key0 = iand(int(ioffset,i8),mask32)
    rnd_epoch = 0
!-- service stream: domain 0
    rnd_state%key = [rnd_key0, 0_i8]
    rnd_state%ctr = 0
    rnd_state%ibuf = 5
!-- per-thread streams: domain = thread index
    do i=1,n
      rnd_states(i)%key = [rnd_key0, int(i,i8)]
      rnd_states(i)%ctr = 0
      rnd_states(i)%ibuf = 5
    enddo
  end subroutine rnd_init
!
!
  pure subroutine rnd_seed_particle(state,ipart)
!-- -------------------------------------------
    implicit none
    type(rnd_t),intent(inout) :: state
    integer,intent(in) :: ipart
!***********************************************************************
! Key a state for one particle in the current transport epoch: the draw
! sequence depends only on (key0, epoch, ipart), not on the thread that
! processes the particle.
!***********************************************************************
    state%key = [rnd_key0, rnd_domain_transport]
    state%ctr = [0_i8, rnd_epoch, iand(int(ipart,i8),mask32), 0_i8]
    state%ibuf = 5
  end subroutine rnd_seed_particle
!
!
  subroutine rnd_advance_epoch
!-- -------------------------
    implicit none
!-- one epoch per transport sweep (iteration); consistent on all ranks
    rnd_epoch = iand(rnd_epoch + 1_i8,mask32)
  end subroutine rnd_advance_epoch
!
!
  pure subroutine rnd_r(x,state)
!-- ----------------------------
    implicit none
    real(dp),intent(out) :: x
    type(rnd_t),intent(inout) :: state
!***********************************************************************
! Draws a uniform real number on (0,1).
!***********************************************************************
    integer(i8) :: w
    call rnd_word(w,state)
    x = (real(w,dp) + 0.5d0)*2d0**(-32)
  end subroutine rnd_r
!
!
  pure subroutine rnd_i(i,state)
!-- -------------------------------
    implicit none
    integer(i4),intent(out) :: i
    type(rnd_t),intent(inout) :: state
!***********************************************************************
! Draws a uniform integer on [0,huge(i4)].
!***********************************************************************
    integer(i8) :: w
    call rnd_word(w,state)
    i = int(iand(w,int(huge(i),i8)),i4)
  end subroutine rnd_i
!
!
  pure subroutine rnd_word(w,state)
!-- ------------------------------
    implicit none
    integer(i8),intent(out) :: w
    type(rnd_t),intent(inout) :: state
!-- return the next buffered 32-bit word, running Philox as needed
    if(state%ibuf > 4) then
      call philox4x32(state%ctr,state%key,state%buf)
!-- increment the block counter (with carry)
      state%ctr(1) = iand(state%ctr(1) + 1_i8,mask32)
      if(state%ctr(1) == 0_i8) &
        state%ctr(4) = iand(state%ctr(4) + 1_i8,mask32)
      state%ibuf = 1
    endif
    w = state%buf(state%ibuf)
    state%ibuf = state%ibuf + 1
  end subroutine rnd_word
!
!
  pure subroutine philox4x32(ctr,key,out)
!-- ------------------------------------
    implicit none
    integer(i8),intent(in) :: ctr(4),key(2)
    integer(i8),intent(out) :: out(4)
!***********************************************************************
! Philox-4x32 with 10 rounds.  All words are 32-bit values in i8.
!***********************************************************************
    integer(i8) :: c0,c1,c2,c3,k0,k1
    integer(i8) :: hi0,lo0,hi1,lo1
    integer :: r
!
    c0 = ctr(1)
    c1 = ctr(2)
    c2 = ctr(3)
    c3 = ctr(4)
    k0 = key(1)
    k1 = key(2)
    do r=1,10
      call mul32(pm0,c0,hi0,lo0)
      call mul32(pm1,c2,hi1,lo1)
      c0 = ieor(ieor(hi1,c1),k0)
      c1 = lo1
      c2 = ieor(ieor(hi0,c3),k1)
      c3 = lo0
      k0 = iand(k0 + pw0,mask32)
      k1 = iand(k1 + pw1,mask32)
    enddo
    out = [c0,c1,c2,c3]
  end subroutine philox4x32
!
!
  pure subroutine mul32(a,b,hi,lo)
!-- -----------------------------
    implicit none
    integer(i8),intent(in) :: a,b
    integer(i8),intent(out) :: hi,lo
!-- 32x32 -> 64 bit product of 32-bit words held in i8, via 16-bit limbs
    integer(i8) :: alo,ahi,blo,bhi,ll,mid
    alo = iand(a,mask16)
    ahi = ishft(a,-16)
    blo = iand(b,mask16)
    bhi = ishft(b,-16)
    ll = alo*blo
    mid = alo*bhi + ahi*blo + ishft(ll,-16)
    lo = ior(ishft(iand(mid,mask16),16),iand(ll,mask16))
    hi = iand(ahi*bhi + ishft(mid,-16),mask32)
  end subroutine mul32
!
!
end module randommod
! vim: fdm=marker
