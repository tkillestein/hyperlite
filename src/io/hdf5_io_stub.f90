!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module hdf5_io
!***********************************************************************
! No-op stand-in for the HDF5 wrappers, compiled when the build has no
! HDF5 Fortran library (HYPERLITE_HDF5=OFF, e.g. the ifx CI job whose
! distro HDF5 modules are gfortran-built).  output.h5 is not written;
! the legacy ASCII writers (in_io_ascii) remain available.
!***********************************************************************
  implicit none
  private

  public :: h5io_create, h5io_close, h5io_mkgroup
  public :: h5io_attr_i, h5io_attr_str
  public :: h5io_write_d1, h5io_write_d2, h5io_write_i3
  public :: h5io_append_d1, h5io_append_i1
  public :: h5io_append_d2, h5io_append_i2, h5io_append_r2
  public :: h5io_ropen, h5io_rclose, h5io_rexists
  public :: h5io_rattr_i
  public :: h5io_read_i1, h5io_read_d1, h5io_read_f1
  public :: h5io_read_d2, h5io_read_f2

!-- whether this build has real HDF5 support
  logical, parameter, public :: h5io_available = .false.
!-- runtime switch: read atomic data from atomic.h5 (set from
!-- in_io_atomdata by inputparmod; lives here to avoid module cycles)
  logical, public :: h5io_atomdata_h5 = .false.

contains

  subroutine h5io_create(fname)
    character(*), intent(in) :: fname
    write(6,*) 'hdf5_io: built without HDF5; ', fname, ' not written'
  end subroutine h5io_create

  subroutine h5io_close
  end subroutine h5io_close

  subroutine h5io_mkgroup(path)
    character(*), intent(in) :: path
    if(.false.) write(6,*) path
  end subroutine h5io_mkgroup

  subroutine h5io_attr_i(path, name, ival)
    character(*), intent(in) :: path, name
    integer, intent(in) :: ival
    if(.false.) write(6,*) path, name, ival
  end subroutine h5io_attr_i

  subroutine h5io_attr_str(path, name, sval)
    character(*), intent(in) :: path, name, sval
    if(.false.) write(6,*) path, name, sval
  end subroutine h5io_attr_str

  subroutine h5io_write_d1(path, arr)
    character(*), intent(in) :: path
    real*8, intent(in) :: arr(:)
    if(.false.) write(6,*) path, arr(1)
  end subroutine h5io_write_d1

  subroutine h5io_write_d2(path, arr)
    character(*), intent(in) :: path
    real*8, intent(in) :: arr(:,:)
    if(.false.) write(6,*) path, arr(1,1)
  end subroutine h5io_write_d2

  subroutine h5io_write_i3(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:,:,:)
    if(.false.) write(6,*) path, arr(1,1,1)
  end subroutine h5io_write_i3

  subroutine h5io_append_d1(path, arr)
    character(*), intent(in) :: path
    real*8, intent(in) :: arr(:)
    if(.false.) write(6,*) path, arr(1)
  end subroutine h5io_append_d1

  subroutine h5io_append_i1(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:)
    if(.false.) write(6,*) path, arr(1)
  end subroutine h5io_append_i1

  subroutine h5io_append_d2(path, arr)
    character(*), intent(in) :: path
    real*8, intent(in) :: arr(:,:)
    if(.false.) write(6,*) path, arr(1,1)
  end subroutine h5io_append_d2

  subroutine h5io_append_i2(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:,:)
    if(.false.) write(6,*) path, arr(1,1)
  end subroutine h5io_append_i2

  subroutine h5io_append_r2(path, arr)
    character(*), intent(in) :: path
    real*4, intent(in) :: arr(:,:)
    if(.false.) write(6,*) path, arr(1,1)
  end subroutine h5io_append_r2

  subroutine h5io_ropen(fname, ok)
    character(*), intent(in) :: fname
    logical, intent(out) :: ok
    ok = .false.
    if(.false.) write(6,*) fname
  end subroutine h5io_ropen

  subroutine h5io_rclose
  end subroutine h5io_rclose

  function h5io_rexists(path) result(lexist)
    character(*), intent(in) :: path
    logical :: lexist
    lexist = .false.
    if(.false.) write(6,*) path
  end function h5io_rexists

  subroutine h5io_rattr_i(path, name, ival)
    character(*), intent(in) :: path, name
    integer, intent(out) :: ival
    ival = 0
    stop 'hdf5_io stub: h5io_rattr_i called in a no-HDF5 build'
    if(.false.) write(6,*) path, name
  end subroutine h5io_rattr_i

  subroutine h5io_read_i1(path, arr)
    character(*), intent(in) :: path
    integer, intent(out) :: arr(:)
    arr = 0
    stop 'hdf5_io stub: h5io_read_i1 called in a no-HDF5 build'
    if(.false.) write(6,*) path
  end subroutine h5io_read_i1

  subroutine h5io_read_d1(path, arr)
    character(*), intent(in) :: path
    real*8, intent(out) :: arr(:)
    arr = 0d0
    stop 'hdf5_io stub: h5io_read_d1 called in a no-HDF5 build'
    if(.false.) write(6,*) path
  end subroutine h5io_read_d1

  subroutine h5io_read_f1(path, arr)
    character(*), intent(in) :: path
    real*4, intent(out) :: arr(:)
    arr = 0.
    stop 'hdf5_io stub: h5io_read_f1 called in a no-HDF5 build'
    if(.false.) write(6,*) path
  end subroutine h5io_read_f1

  subroutine h5io_read_d2(path, arr)
    character(*), intent(in) :: path
    real*8, intent(out) :: arr(:,:)
    arr = 0d0
    stop 'hdf5_io stub: h5io_read_d2 called in a no-HDF5 build'
    if(.false.) write(6,*) path
  end subroutine h5io_read_d2

  subroutine h5io_read_f2(path, arr)
    character(*), intent(in) :: path
    real*4, intent(out) :: arr(:,:)
    arr = 0.
    stop 'hdf5_io stub: h5io_read_f2 called in a no-HDF5 build'
    if(.false.) write(6,*) path
  end subroutine h5io_read_f2

end module hdf5_io
! vim: fdm=marker
