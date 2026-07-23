*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine sourcetally
c     ------------------------
      use gridmod
      use groupmod
      use sourcemod
      use particlemod
      use miscmod
      implicit none
************************************************************************
* This routine tallies the source particles for output
************************************************************************
      integer :: ipart
      integer :: ix,iy,iz,ic,ig
c
c-- init
      src_energy = 0.d0
      src_number = 0
      if(.not.allocated(src_energy)) stop "src_energy not allocated"
c
!!c$omp do schedule(static,1) !round-robin
      do ipart=1,prt_npartmax
c-- check vacancy
      if(prt_isvacant(ipart)) cycle
c
c-- retrieving group and cell id
       ig = binsrch(prt_wl(ipart),grp_wl,grp_ng+1,.false.)
       ix = binsrch(prt_x(ipart),grd_xarr,grd_nx+1,.false.)
       iy = binsrch(prt_y0,grd_yarr,grd_ny+1,.false.)
       iz = binsrch(prt_z0,grd_zarr,grd_nz+1,.false.)
       ic = grd_icell(ix,iy,iz)
c
c-- tallying source packets
       if(ic<1.or.ic>grd_ncell) stop 'ic outside range'
       if(ig<1.or.ic>grp_ng) stop 'ig outside range'
       src_energy(ic,ig) = src_energy(ic,ig)+prt_e(ipart)
       src_number(ic,ig) = src_number(ic,ig)+1
c
      enddo !ipart
!!c$omp end do nowait
c
      end subroutine sourcetally
c vim: fdm=marker
