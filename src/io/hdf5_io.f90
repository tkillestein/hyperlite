!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module hdf5_io
  use kindmod
!***********************************************************************
! Thin wrappers over the HDF5 Fortran interface for the hyperlite
! output file (output.h5) and, later, atomic-data input.
!
! Two dataset flavors:
!  - h5io_write_*  : fixed-size datasets written once (grids, axes)
!  - h5io_append_* : datasets with an unlimited trailing (iteration)
!                    dimension, extended by one slice per call
!***********************************************************************
  use hdf5
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

!-- whether this build has real HDF5 support (stub sets .false.)
  logical, parameter, public :: h5io_available = .true.
!-- runtime switch: read atomic data from atomic.h5 (set from
!-- in_io_atomdata by inputparmod; lives here to avoid module cycles)
  logical, public :: h5io_atomdata_h5 = .false.

  integer(hid_t) :: fid = -1
  logical :: lopen = .false.
!-- separate handle for read-only files (atomic data)
  integer(hid_t) :: rfid = -1
  logical :: lropen = .false.

contains

  subroutine chkerr(ierr, msg)
    integer, intent(in) :: ierr
    character(*), intent(in) :: msg
    if(ierr /= 0) then
      write(0,*) 'hdf5_io error: ', msg
      stop 'hdf5_io: hdf5 call failed'
    endif
  end subroutine chkerr

!-- file lifecycle
!=================
  subroutine h5io_create(fname)
    character(*), intent(in) :: fname
    integer :: ierr
    call h5open_f(ierr)
    call chkerr(ierr, 'h5open '//fname)
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, ierr)
    call chkerr(ierr, 'create '//fname)
    lopen = .true.
  end subroutine h5io_create

  subroutine h5io_close
    integer :: ierr
    if(.not.lopen) return
    call h5fclose_f(fid, ierr)
    call chkerr(ierr, 'close file')
    call h5close_f(ierr)
    lopen = .false.
  end subroutine h5io_close

  subroutine h5io_mkgroup(path)
    character(*), intent(in) :: path
    integer(hid_t) :: gid
    integer :: ierr
    call h5gcreate_f(fid, path, gid, ierr)
    call chkerr(ierr, 'mkgroup '//path)
    call h5gclose_f(gid, ierr)
  end subroutine h5io_mkgroup

!-- attributes
!=============
  subroutine h5io_attr_i(path, name, ival)
    character(*), intent(in) :: path, name
    integer, intent(in) :: ival
    integer(hid_t) :: oid, sid, aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call h5oopen_f(fid, path, oid, ierr)
    call chkerr(ierr, 'attr open '//path)
    call h5screate_f(H5S_SCALAR_F, sid, ierr)
    call h5acreate_f(oid, name, H5T_NATIVE_INTEGER, sid, aid, ierr)
    call chkerr(ierr, 'attr create '//name)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, ival, dims, ierr)
    call chkerr(ierr, 'attr write '//name)
    call h5aclose_f(aid, ierr)
    call h5sclose_f(sid, ierr)
    call h5oclose_f(oid, ierr)
  end subroutine h5io_attr_i

  subroutine h5io_attr_str(path, name, sval)
    character(*), intent(in) :: path, name, sval
    integer(hid_t) :: oid, sid, aid, tid
    integer(hsize_t) :: dims(1)
    integer(size_t) :: slen
    integer :: ierr
    dims(1) = 1
    slen = max(len_trim(sval), 1)
    call h5oopen_f(fid, path, oid, ierr)
    call chkerr(ierr, 'attr open '//path)
    call h5tcopy_f(H5T_FORTRAN_S1, tid, ierr)
    call h5tset_size_f(tid, slen, ierr)
    call h5screate_f(H5S_SCALAR_F, sid, ierr)
    call h5acreate_f(oid, name, tid, sid, aid, ierr)
    call chkerr(ierr, 'attr create '//name)
    call h5awrite_f(aid, tid, sval, dims, ierr)
    call chkerr(ierr, 'attr write '//name)
    call h5aclose_f(aid, ierr)
    call h5sclose_f(sid, ierr)
    call h5tclose_f(tid, ierr)
    call h5oclose_f(oid, ierr)
  end subroutine h5io_attr_str

!-- fixed-size datasets
!======================
  subroutine h5io_write_d1(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(in) :: arr(:)
    integer(hid_t) :: sid, did
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = size(arr)
    call h5screate_simple_f(1, dims, sid, ierr)
    call h5dcreate_f(fid, path, H5T_NATIVE_DOUBLE, sid, did, ierr)
    call chkerr(ierr, 'write_d1 '//path)
    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call chkerr(ierr, 'write_d1 data '//path)
    call h5dclose_f(did, ierr)
    call h5sclose_f(sid, ierr)
  end subroutine h5io_write_d1

  subroutine h5io_write_d2(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(in) :: arr(:,:)
    integer(hid_t) :: sid, did
    integer(hsize_t) :: dims(2)
    integer :: ierr
    dims = shape(arr)
    call h5screate_simple_f(2, dims, sid, ierr)
    call h5dcreate_f(fid, path, H5T_NATIVE_DOUBLE, sid, did, ierr)
    call chkerr(ierr, 'write_d2 '//path)
    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call chkerr(ierr, 'write_d2 data '//path)
    call h5dclose_f(did, ierr)
    call h5sclose_f(sid, ierr)
  end subroutine h5io_write_d2

  subroutine h5io_write_i3(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:,:,:)
    integer(hid_t) :: sid, did
    integer(hsize_t) :: dims(3)
    integer :: ierr
    dims = shape(arr)
    call h5screate_simple_f(3, dims, sid, ierr)
    call h5dcreate_f(fid, path, H5T_NATIVE_INTEGER, sid, did, ierr)
    call chkerr(ierr, 'write_i3 '//path)
    call h5dwrite_f(did, H5T_NATIVE_INTEGER, arr, dims, ierr)
    call chkerr(ierr, 'write_i3 data '//path)
    call h5dclose_f(did, ierr)
    call h5sclose_f(sid, ierr)
  end subroutine h5io_write_i3

!-- appendable datasets: fixed leading dims + unlimited iteration dim
!====================================================================
!-- open-or-create an appendable dataset and select the slab for the
!-- next iteration slice; caller writes and then calls append_done
  subroutine append_slab(path, tid, rank, dfix, did, fspace, mspace)
    character(*), intent(in) :: path
    integer(hid_t), intent(in) :: tid
    integer, intent(in) :: rank            !total rank incl. iteration dim
    integer(hsize_t), intent(in) :: dfix(rank-1) !fixed leading dims
    integer(hid_t), intent(out) :: did, fspace, mspace
    integer(hid_t) :: sid, pid
    integer(hsize_t) :: dims(rank), maxdims(rank), offset(rank)
    integer :: ierr
    logical :: lexist
    call h5lexists_f(fid, path, lexist, ierr)
    call chkerr(ierr, 'append exists '//path)
    if(.not.lexist) then
!-- create with one iteration slice
      dims(:rank-1) = dfix
      dims(rank) = 1
      maxdims(:rank-1) = dfix
      maxdims(rank) = H5S_UNLIMITED_F
      call h5screate_simple_f(rank, dims, sid, ierr, maxdims)
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, ierr)
      call h5pset_chunk_f(pid, rank, dims, ierr)
      call h5dcreate_f(fid, path, tid, sid, did, ierr, pid)
      call chkerr(ierr, 'append create '//path)
      call h5pclose_f(pid, ierr)
      call h5sclose_f(sid, ierr)
      offset = 0
    else
!-- extend by one iteration slice
      call h5dopen_f(fid, path, did, ierr)
      call chkerr(ierr, 'append open '//path)
      call h5dget_space_f(did, sid, ierr)
      call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)
      call h5sclose_f(sid, ierr)
      offset = 0
      offset(rank) = dims(rank)
      dims(rank) = dims(rank) + 1
      call h5dset_extent_f(did, dims, ierr)
      call chkerr(ierr, 'append extend '//path)
    endif
!-- file-space slab for this slice + matching memory space
    dims(:rank-1) = dfix
    dims(rank) = 1
    call h5dget_space_f(did, fspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, offset, dims, &
         ierr)
    call chkerr(ierr, 'append slab '//path)
    call h5screate_simple_f(rank, dims, mspace, ierr)
  end subroutine append_slab

  subroutine append_done(did, fspace, mspace)
    integer(hid_t), intent(in) :: did, fspace, mspace
    integer :: ierr
    call h5sclose_f(mspace, ierr)
    call h5sclose_f(fspace, ierr)
    call h5dclose_f(did, ierr)
    call h5fflush_f(fid, H5F_SCOPE_LOCAL_F, ierr)
  end subroutine append_done

  subroutine h5io_append_d1(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(in) :: arr(:)
    integer(hid_t) :: did, fspace, mspace
    integer(hsize_t) :: dfix(1), dims(2)
    integer :: ierr
    dfix(1) = size(arr)
    dims(:1) = dfix
    dims(2) = 1
    call append_slab(path, H5T_NATIVE_DOUBLE, 2, dfix, did, fspace, &
         mspace)
    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr, mspace, &
         fspace)
    call chkerr(ierr, 'append_d1 '//path)
    call append_done(did, fspace, mspace)
  end subroutine h5io_append_d1

  subroutine h5io_append_i1(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:)
    integer(hid_t) :: did, fspace, mspace
    integer(hsize_t) :: dfix(1), dims(2)
    integer :: ierr
    dfix(1) = size(arr)
    dims(:1) = dfix
    dims(2) = 1
    call append_slab(path, H5T_NATIVE_INTEGER, 2, dfix, did, fspace, &
         mspace)
    call h5dwrite_f(did, H5T_NATIVE_INTEGER, arr, dims, ierr, mspace, &
         fspace)
    call chkerr(ierr, 'append_i1 '//path)
    call append_done(did, fspace, mspace)
  end subroutine h5io_append_i1

  subroutine h5io_append_d2(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(in) :: arr(:,:)
    integer(hid_t) :: did, fspace, mspace
    integer(hsize_t) :: dfix(2), dims(3)
    integer :: ierr
    dfix = shape(arr)
    dims(:2) = dfix
    dims(3) = 1
    call append_slab(path, H5T_NATIVE_DOUBLE, 3, dfix, did, fspace, &
         mspace)
    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr, mspace, &
         fspace)
    call chkerr(ierr, 'append_d2 '//path)
    call append_done(did, fspace, mspace)
  end subroutine h5io_append_d2

  subroutine h5io_append_i2(path, arr)
    character(*), intent(in) :: path
    integer, intent(in) :: arr(:,:)
    integer(hid_t) :: did, fspace, mspace
    integer(hsize_t) :: dfix(2), dims(3)
    integer :: ierr
    dfix = shape(arr)
    dims(:2) = dfix
    dims(3) = 1
    call append_slab(path, H5T_NATIVE_INTEGER, 3, dfix, did, fspace, &
         mspace)
    call h5dwrite_f(did, H5T_NATIVE_INTEGER, arr, dims, ierr, mspace, &
         fspace)
    call chkerr(ierr, 'append_i2 '//path)
    call append_done(did, fspace, mspace)
  end subroutine h5io_append_i2

  subroutine h5io_append_r2(path, arr)
    character(*), intent(in) :: path
    real(sp), intent(in) :: arr(:,:)
    integer(hid_t) :: did, fspace, mspace
    integer(hsize_t) :: dfix(2), dims(3)
    integer :: ierr
    dfix = shape(arr)
    dims(:2) = dfix
    dims(3) = 1
    call append_slab(path, H5T_NATIVE_REAL, 3, dfix, did, fspace, &
         mspace)
    call h5dwrite_f(did, H5T_NATIVE_REAL, arr, dims, ierr, mspace, &
         fspace)
    call chkerr(ierr, 'append_r2 '//path)
    call append_done(did, fspace, mspace)
  end subroutine h5io_append_r2

!-- read-only access (atomic data)
!================================
  subroutine h5io_ropen(fname, ok)
    character(*), intent(in) :: fname
    logical, intent(out) :: ok
    integer :: ierr
    call h5open_f(ierr)
    call chkerr(ierr, 'h5open '//fname)
    call h5eset_auto_f(0, ierr) !suppress error stack for missing file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, rfid, ierr)
    call h5eset_auto_f(1, ierr)
    ok = ierr == 0
    lropen = ok
  end subroutine h5io_ropen

  subroutine h5io_rclose
    integer :: ierr
    if(.not.lropen) return
    call h5fclose_f(rfid, ierr)
    call chkerr(ierr, 'rclose')
    lropen = .false.
!-- close the library only if no write file is open
    if(.not.lopen) call h5close_f(ierr)
  end subroutine h5io_rclose

  function h5io_rexists(path) result(lexist)
    character(*), intent(in) :: path
    logical :: lexist
    integer :: ierr, i, n
    character(len(path)) :: sub
!-- h5lexists on a nested path requires each parent link to exist
    n = len_trim(path)
    lexist = .true.
    do i = 2, n
      if(path(i:i) == '/') then
        sub = path(:i-1)
        call h5lexists_f(rfid, trim(sub), lexist, ierr)
        if(ierr /= 0 .or. .not.lexist) then
          lexist = .false.
          return
        endif
      endif
    enddo
    call h5lexists_f(rfid, path, lexist, ierr)
    if(ierr /= 0) lexist = .false.
  end function h5io_rexists

  subroutine h5io_rattr_i(path, name, ival)
    character(*), intent(in) :: path, name
    integer, intent(out) :: ival
    integer(hid_t) :: oid, aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call h5oopen_f(rfid, path, oid, ierr)
    call chkerr(ierr, 'rattr open '//path)
    call h5aopen_f(oid, name, aid, ierr)
    call chkerr(ierr, 'rattr '//name)
    call h5aread_f(aid, H5T_NATIVE_INTEGER, ival, dims, ierr)
    call chkerr(ierr, 'rattr read '//name)
    call h5aclose_f(aid, ierr)
    call h5oclose_f(oid, ierr)
  end subroutine h5io_rattr_i

  subroutine read_open(path, rank, did, dims)
    character(*), intent(in) :: path
    integer, intent(in) :: rank
    integer(hid_t), intent(out) :: did
    integer(hsize_t), intent(out) :: dims(rank)
    integer(hid_t) :: sid
    integer(hsize_t) :: maxdims(rank)
    integer :: ierr
    call h5dopen_f(rfid, path, did, ierr)
    call chkerr(ierr, 'read open '//path)
    call h5dget_space_f(did, sid, ierr)
    call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)
    call h5sclose_f(sid, ierr)
  end subroutine read_open

  subroutine h5io_read_i1(path, arr)
    character(*), intent(in) :: path
    integer, intent(out) :: arr(:)
    integer(hid_t) :: did
    integer(hsize_t) :: dims(1)
    integer :: ierr
    call read_open(path, 1, did, dims)
    if(dims(1) /= size(arr)) call chkerr(1, 'read_i1 size '//path)
    call h5dread_f(did, H5T_NATIVE_INTEGER, arr, dims, ierr)
    call chkerr(ierr, 'read_i1 '//path)
    call h5dclose_f(did, ierr)
  end subroutine h5io_read_i1

  subroutine h5io_read_d1(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(out) :: arr(:)
    integer(hid_t) :: did
    integer(hsize_t) :: dims(1)
    integer :: ierr
    call read_open(path, 1, did, dims)
    if(dims(1) /= size(arr)) call chkerr(1, 'read_d1 size '//path)
    call h5dread_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call chkerr(ierr, 'read_d1 '//path)
    call h5dclose_f(did, ierr)
  end subroutine h5io_read_d1

  subroutine h5io_read_f1(path, arr)
    character(*), intent(in) :: path
    real(sp), intent(out) :: arr(:)
    integer(hid_t) :: did
    integer(hsize_t) :: dims(1)
    integer :: ierr
    call read_open(path, 1, did, dims)
    if(dims(1) /= size(arr)) call chkerr(1, 'read_f1 size '//path)
    call h5dread_f(did, H5T_NATIVE_REAL, arr, dims, ierr)
    call chkerr(ierr, 'read_f1 '//path)
    call h5dclose_f(did, ierr)
  end subroutine h5io_read_f1

  subroutine h5io_read_d2(path, arr)
    character(*), intent(in) :: path
    real(dp), intent(out) :: arr(:,:)
    integer(hid_t) :: did
    integer(hsize_t) :: dims(2)
    integer :: ierr
    call read_open(path, 2, did, dims)
    if(any(dims /= shape(arr))) call chkerr(1, 'read_d2 size '//path)
    call h5dread_f(did, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call chkerr(ierr, 'read_d2 '//path)
    call h5dclose_f(did, ierr)
  end subroutine h5io_read_d2

  subroutine h5io_read_f2(path, arr)
    character(*), intent(in) :: path
    real(sp), intent(out) :: arr(:,:)
    integer(hid_t) :: did
    integer(hsize_t) :: dims(2)
    integer :: ierr
    call read_open(path, 2, did, dims)
    if(any(dims /= shape(arr))) call chkerr(1, 'read_f2 size '//path)
    call h5dread_f(did, H5T_NATIVE_REAL, arr, dims, ierr)
    call chkerr(ierr, 'read_f2 '//path)
    call h5dclose_f(did, ierr)
  end subroutine h5io_read_f2

end module hdf5_io
! vim: fdm=marker
