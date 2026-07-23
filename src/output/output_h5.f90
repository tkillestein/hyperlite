!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine output_h5
  use kindmod
!***********************************************************************
! Write the self-describing HDF5 output file (output.h5), one group per
! physical quantity (see OVERHAUL.md section 3).  Static geometry/axes
! are written on the first call; per-iteration fields are appended along
! an unlimited trailing dimension.  Mirrors the legacy ASCII writers
! (output_grid/output_grp/output_flux/output_source + tot_energy).
!***********************************************************************
  use hdf5_io
  use versionmod
  use mpimod, only: nmpi
  use inputparmod
  use timingmod
  use gridmod
  use groupmod
  use fluxmod
  use totalsmod
  use sourcemod
  implicit none
!
  logical, save :: lfirst = .true.
  real(dp) :: t0, t1
  real(dp) :: arr(grd_ncell), tot(4)
!
  t0 = t_time()
!
  if(lfirst) then
    call h5io_create('output.h5')
!-- meta
    call h5io_mkgroup('/meta')
    call h5io_attr_str('/meta', 'code', 'hyperlite')
    call h5io_attr_str('/meta', 'version', ver_release)
    call h5io_attr_str('/meta', 'git_revision', ver_git)
    call h5io_attr_i('/meta', 'schema_version', 1)
    call h5io_attr_i('/meta', 'nmpi', nmpi)
    call h5io_attr_i('/meta', 'nomp', in_nomp)
    call h5io_attr_i('/meta', 'rnd_seed', in_rnd_seed)
    call h5io_attr_i('/meta', 'niter', in_niter)
!-- grid: geometry + static axes
    call h5io_mkgroup('/grid')
    call h5io_attr_i('/grid', 'igeom', grd_igeom)
    call h5io_attr_i('/grid', 'nx', grd_nx)
    call h5io_attr_i('/grid', 'ny', grd_ny)
    call h5io_attr_i('/grid', 'nz', grd_nz)
    call h5io_write_d1('/grid/xarr', grd_xarr)
    call h5io_write_d1('/grid/yarr', grd_yarr)
    call h5io_write_d1('/grid/zarr', grd_zarr)
    call h5io_write_i3('/grid/icell', grd_icell)
!-- group: wavelength grid
    call h5io_mkgroup('/group')
    call h5io_attr_i('/group', 'ng', grp_ng)
    call h5io_write_d1('/group/wl', grp_wl)
!-- flux: axes
    call h5io_mkgroup('/flux')
    call h5io_attr_i('/flux', 'ng', flx_ng)
    call h5io_attr_i('/flux', 'nmu', flx_nmu)
    call h5io_attr_i('/flux', 'nom', flx_nom)
    call h5io_write_d1('/flux/wl', flx_wl)
    call h5io_write_d1('/flux/mu', flx_mu)
    call h5io_write_d1('/flux/om', flx_om)
!-- per-iteration groups
    call h5io_mkgroup('/total')
    call h5io_mkgroup('/source')
    lfirst = .false.
  endif
!
!-- energy totals [eout,evelo,sflux,sthermal] (legacy output.tot_energy)
  tot(1) = tot_eout
  tot(2) = tot_evelo
  tot(3) = tot_sflux
  tot(4) = tot_sthermal
  call h5io_append_d1('/total/energy', tot)
!
!-- grid fields (legacy output.grd_*)
  if(.not.in_io_nogriddump) then
    call h5io_append_i1('/grid/nvol', grd_nvol)
    arr = 1d0/grd_tempinv
    call h5io_append_d1('/grid/temp', arr)
    call h5io_append_d1('/grid/radtemp', grd_radtemp)
    if(trim(in_io_opacdump)/='off') then
      call h5io_append_d1('/grid/capgrey', grd_capgrey)
      call h5io_append_d1('/grid/capemitgrey', grd_capemitgrey)
      call h5io_append_d1('/grid/capross', grd_capross)
      call h5io_append_d1('/grid/sig', grd_sig)
    endif
    arr = grd_tally/grd_vol
    call h5io_append_d1('/grid/eraddens', arr)
    if(in_io_dogrdtally) then
      call h5io_append_i1('/grid/methodswap', grd_methodswap)
    endif
  endif
!
!-- multigroup opacities (legacy output.grp_*)
  if(trim(in_io_opacdump)/='off') then
    call h5io_append_r2('/group/cap', grd_cap(:grp_ng,:))
    call h5io_append_r2('/group/capemit', grd_capemit(:grp_ng,:))
    call h5io_append_r2('/group/emiss', grd_emiss(:grp_ng,:))
  endif
!
!-- flux (legacy output.flx_*): [ng, nmu*nom, iter]
  call h5io_append_d2('/flux/luminos', &
       reshape(flx_luminos, [flx_ng, flx_nmu*flx_nom]))
  call h5io_append_i2('/flux/lumnum', &
       reshape(flx_lumnum, [flx_ng, flx_nmu*flx_nom]))
  call h5io_append_d2('/flux/lumdev', &
       reshape(flx_lumdev, [flx_ng, flx_nmu*flx_nom]))
!
!-- source (legacy output.src_*): [ncell, ng, iter]
  call h5io_append_i2('/source/number', src_number)
  call h5io_append_d2('/source/energy', src_energy)
!
!-- timing
  t1 = t_time()
  call timereg(t_output, t1-t0)
!
end subroutine output_h5


subroutine output_h5_profile
  use kindmod
!***********************************************************************
! Write the final ejecta-structure profile to /profile (legacy
! output.profile), one dataset per column plus the mass-fraction table.
!***********************************************************************
  use hdf5_io
  use elemdatamod
  use miscmod, only: lcase
  use gridmod
  use gasmod, only: gas_nelem
  implicit none
  integer :: l
  character(4) :: el_name
  character(3*gas_nelem) :: el_list
  real(dp) :: arr(grd_ncell)
!
  call h5io_mkgroup('/profile')
  call h5io_write_d1('/profile/x_left', grd_xarr(:grd_ncell))
  call h5io_write_d1('/profile/x_right', grd_xarr(2:grd_ncell+1))
  call h5io_write_d1('/profile/vx_left', grd_vxarr(:grd_ncell))
  call h5io_write_d1('/profile/vx_right', grd_vxarr(2:grd_ncell+1))
  call h5io_write_d1('/profile/mass', grd_mass)
  call h5io_write_d1('/profile/rho', grd_rho)
  call h5io_write_d1('/profile/vol', grd_vol)
  arr = 1d0/grd_tempinv
  call h5io_write_d1('/profile/avg_temp', arr)
  call h5io_write_d1('/profile/rad_temp', grd_radtemp)
  call h5io_write_d1('/profile/ye', grd_ye)
  arr = grd_nelec*grd_natom/grd_vol
  call h5io_write_d1('/profile/n_e', arr)
  arr = grd_natom/grd_vol
  call h5io_write_d1('/profile/n_atom', arr)
!-- mass fractions for all elements, with the symbol list as attribute
  call h5io_write_d2('/profile/massfr', grd_massfr(:gas_nelem,:))
  el_list = ''
  do l = 1, gas_nelem
    el_name = lcase(trim(elem_data(l)%sym))
    el_list = trim(el_list)//' '//trim(el_name)
  enddo
  call h5io_attr_str('/profile/massfr', 'elements', &
       trim(adjustl(el_list)))
!
end subroutine output_h5_profile


subroutine output_h5_close
  use kindmod
!-- close output.h5 (call once after the iteration loop, master rank)
  use hdf5_io
  implicit none
  call h5io_close
end subroutine output_h5_close
! vim: fdm=marker
