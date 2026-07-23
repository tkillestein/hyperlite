!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module kindmod
!***********************************************************************
! Named kind parameters (iso_fortran_env), replacing the nonstandard
! real*8/integer*2/logical*2 declarations (OVERHAUL.md Phase 4).
!***********************************************************************
  use iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  implicit none

  integer, parameter :: sp = real32  !single precision (opacity tables)
  integer, parameter :: dp = real64  !double precision (default)

  integer, parameter :: i1 = int8   !eof-probe byte reads
  integer, parameter :: i2 = int16  !compact ion codes (bb_xs)
  integer, parameter :: i4 = int32  !mzran RNG words
  integer, parameter :: i8 = int64  !particle/vacancy counters

!-- 2-byte logical for the cache-packed group lumping masks; no standard
!-- named constant exists, kind=2 is supported by gfortran and ifx
  integer, parameter :: lk2 = 2

end module kindmod
! vim: fdm=marker
