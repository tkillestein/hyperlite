*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module gridmod
      use kindmod
c     --------------
      implicit none
c
      integer :: grd_igeom=0
c
      integer,private :: ng=0
c
      integer :: grd_nep=0 !number of emission probability bins
      integer :: grd_nepg=0 !number of groups per emission probability bin
c
c-- complete domain
      integer :: grd_nx=0
      integer :: grd_ny=0
      integer :: grd_nz=0
c-- spatial arrays
      real(dp),allocatable :: grd_xarr(:)  !(nx+1), left cell edge values
      real(dp),allocatable :: grd_yarr(:)  !(ny+1), left cell edge values
      real(dp),allocatable :: grd_zarr(:)  !(nz+1), left cell edge values
c-- velocity array
      real(dp),allocatable :: grd_vxarr(:)  !(nx+1), left cell edge values
c
c-- pointer into compressed domain
      integer,allocatable :: grd_icell(:,:,:) !(nx,ny,nz)
c
c-- domain decomposition
      integer :: grd_idd1=0
      integer :: grd_ndd=0
c
c-- compressed domain
      integer :: grd_ncell=0  !number of cells
c
c-- Probability of emission in a given zone and group
      real(dp),allocatable :: grd_emitprob(:,:) !(nep,ncell)

c-- Line+Cont extinction coeff
      real(sp),allocatable :: grd_cap(:,:) !(ng,ncell)
c-- Line+Cont "emission opacity" for NLTE calculations
      real(sp),allocatable :: grd_capemit(:,:) !(ng,ncell)
c-- Line emissivity for output
      real(sp),allocatable :: grd_emiss(:,:) !(ng,ncell)
c-- leakage opacities
      real(dp),allocatable :: grd_opaclump(:,:) !(11,ncell) leak(6),speclump,caplump,capemitlump,igemitmax,doplump
c-- DDMC per-cell transport tables, rebuilt once per transport sweep by
c-- ddmc_tables (opacity/temperature are constant during the sweep)
      real(dp),allocatable :: grd_specarr(:,:) !(ng,ncell) Planck spectral weights
      integer,allocatable :: grd_nlump(:) !(ncell) number of lumped groups
      integer(i2),allocatable :: grd_glumps(:,:) !(ng,ncell) group order: lumped first
      logical(lk2),allocatable :: grd_llumps(:,:) !(ng,ncell) group is in the lump
      real(dp),allocatable :: grd_emitlump(:) !(ncell)
      real(dp),allocatable :: grd_caplump(:) !(ncell)
      real(dp),allocatable :: grd_capemitlump(:) !(ncell)
      real(dp),allocatable :: grd_doplump(:) !(ncell)
      real(dp),allocatable :: grd_capgreyinv(:) !(ncell)
      real(dp),allocatable :: grd_capemitgreyinv(:) !(ncell)
c-- cumulative group-sampling tables: the exact partial sums the
c-- diffusion11 scans accumulate, in scan order (see ddmc_tables)
      real(dp),allocatable :: grd_leaklcum(:,:) !(ng,ncell) left leakage
      real(dp),allocatable :: grd_leakrcum(:,:) !(ng,ncell) right leakage
      real(dp),allocatable :: grd_scatcum(:,:) !(ng,ncell) effective scattering
      real(dp),allocatable :: grd_dopcum(:,:) !(ng,ncell) doppler shift
      integer(i2),allocatable :: grd_dopidx(:,:) !(ng,ncell) group tested per doppler entry
      integer,allocatable :: grd_ndop(:) !(ncell) number of doppler entries
      real(dp),allocatable :: grd_tempinv(:) !(ncell)
c-- scattering coefficient
      real(dp),allocatable :: grd_sig(:) !(ncell) !grey scattering opacity
c-- Planck opacity (gray)
      real(dp),allocatable :: grd_capgrey(:) !(ncell)
c-- "Emission" opacity (gray) for NLTE calcualtions
      real(dp),allocatable :: grd_capemitgrey(:) !(ncell)
c-- Rosseland mean opacity
      real(dp),allocatable :: grd_capross(:) !(ncell)
c-- For output
      real(dp),allocatable :: grd_natom(:) !(ncell) !number of atoms
      real(dp),allocatable :: grd_nelec(:) !(ncell) !number of electrons
      real(dp),allocatable :: grd_rho(:) !(ncell)
      real(dp),allocatable :: grd_mass(:) !(ncell)
      real(dp),allocatable :: grd_radtemp(:) !(ncell)
      real(dp),allocatable :: grd_lum(:) !(ncell)
      real(dp),allocatable :: grd_ye(:) !(ncell) !electron fraction
      real(dp),allocatable :: grd_massfr(:,:) !(nelem,ncell)
c-- cell centered radii
      real(dp),allocatable :: grd_rcell(:) !(ncell)
c-- divergence of velocity
      real(dp),allocatable :: grd_divv(:) !(ncell)
c-- radial velocity gradient dvx/dx (permanent; = interp slope = dvdr)
      real(dp),allocatable :: grd_dvdx(:) !(ncell)
c-- radiation energy density for tally
      real(dp),allocatable :: grd_tally(:)   !(ncell) (eraddens)
c-- amplification factor excess
      real(dp),allocatable :: grd_eamp(:)   !(ncell)


c-- number of IMC-DDMC method changes per cell per time step
      integer,allocatable :: grd_methodswap(:) !(ncell)
c-- number of census prt_particles per cell

c
c-- packet number and energy distribution
c========================================
c
      real(dp),allocatable :: grd_vol(:)  !(ncell)
c
      integer,allocatable :: grd_nvol(:) !(ncell) number of thermal source particles generated per cell
c
      real(dp),allocatable :: grd_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
      real(dp),allocatable :: grd_emitex(:) !(ncell) amount of external energy emitted per cell in a time step
      real(dp),allocatable :: grd_evolinit(:) !(ncell) amount of initial energy per cell per group
c
      real(dp),allocatable :: grd_jrad(:,:) !(ng,ncell) radiation intensity per cell per group
c
      interface
      pure function emitgroup(r,ic) result(ig)
      use kindmod
      integer :: ig
      real(dp),intent(in) :: r
      integer,intent(in) :: ic
      end function emitgroup
      end interface
c
      save
c
      contains
c
      subroutine gridmod_init(ltalk,ngin,ncell,idd1,ndd)
c     --------------------------------------------------
      implicit none
      logical,intent(in) :: ltalk
      integer,intent(in) :: ngin
      integer,intent(in) :: ncell,idd1,ndd
************************************************************************
* Allocate grd variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      ng = ngin
c
c-- emission probability
      grd_nep = nint(sqrt(dble(ng)))
      grd_nepg = ceiling(ng/(grd_nep + 1d0))
c
c-- number of non-void cells, plus one optional dummy cell if void cells exist
      grd_ncell = ncell
      grd_idd1 = idd1
      grd_ndd = ndd
c
      allocate(grd_xarr(grd_nx+1))
      allocate(grd_yarr(grd_ny+1))
      allocate(grd_zarr(grd_nz+1))
c
      allocate(grd_vxarr(grd_nx+1))
c
c-- complete domain
      allocate(grd_icell(grd_nx,grd_ny,grd_nz))
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = int((int(grd_ncell,8)*(8*(28+11) + 5*4))/1024) !kB
       write(6,*) 'ALLOC grd      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
       n = int((int(grd_ncell,8)*4*(2+ng))/1024) !kB
       write(6,*) 'ALLOC grd_cap  :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif
c
c-- ndim=3 alloc
      allocate(grd_tally(grd_ncell))
      allocate(grd_eamp(grd_ncell))
      allocate(grd_capgrey(grd_ncell))
      allocate(grd_capemitgrey(grd_ncell))
      allocate(grd_capross(grd_ncell))
      allocate(grd_sig(grd_ncell))
      allocate(grd_tempinv(grd_ncell))
      allocate(grd_natom(grd_ncell))
      allocate(grd_nelec(grd_ncell))
      allocate(grd_vol(grd_ncell))
      allocate(grd_rho(grd_ncell))
      allocate(grd_mass(grd_ncell))
      allocate(grd_radtemp(grd_ncell))
      allocate(grd_lum(grd_ncell))
      allocate(grd_ye(grd_ncell))
      allocate(grd_rcell(grd_ncell))
      allocate(grd_divv(grd_ncell))
      allocate(grd_dvdx(grd_ncell))
c
      allocate(grd_emit(grd_ncell))
      grd_emit = 0d0
      allocate(grd_emitex(grd_ncell))
      allocate(grd_evolinit(grd_ncell))
      allocate(grd_jrad(ng,grd_ncell))
c
c-- ndim=3 integer
      allocate(grd_nvol(grd_ncell))
      grd_nvol = 0
c
      allocate(grd_methodswap(grd_ncell))
c
c-- DDMC per-cell tables
      allocate(grd_specarr(ng,grd_ncell))
      allocate(grd_nlump(grd_ncell))
      allocate(grd_glumps(ng,grd_ncell))
      allocate(grd_llumps(ng,grd_ncell))
      allocate(grd_emitlump(grd_ncell))
      allocate(grd_caplump(grd_ncell))
      allocate(grd_capemitlump(grd_ncell))
      allocate(grd_doplump(grd_ncell))
      allocate(grd_capgreyinv(grd_ncell))
      allocate(grd_capemitgreyinv(grd_ncell))
      allocate(grd_leaklcum(ng,grd_ncell))
      allocate(grd_leakrcum(ng,grd_ncell))
      allocate(grd_scatcum(ng,grd_ncell))
      allocate(grd_dopcum(ng,grd_ncell))
      allocate(grd_dopidx(ng,grd_ncell))
      allocate(grd_ndop(grd_ncell))
c
c-- ndim=4 alloc
      allocate(grd_opaclump(11,grd_ncell))
      allocate(grd_emitprob(grd_nep,grd_ncell))
c-- ndim=4 alloc
      allocate(grd_cap(ng,grd_ncell))
      allocate(grd_capemit(ng,grd_ncell))
      allocate(grd_emiss(ng,grd_ncell))
c
      end subroutine gridmod_init
c
c
      subroutine grid_dealloc
      deallocate(grd_xarr)
      deallocate(grd_yarr)
      deallocate(grd_zarr)
c-- velocity arrays
      deallocate(grd_vxarr)
c-- complete domain
      deallocate(grd_icell)
c-- gasmod
      deallocate(grd_tally)
      deallocate(grd_eamp)
      deallocate(grd_capgrey)
      deallocate(grd_capemitgrey)
      deallocate(grd_capross)
      deallocate(grd_sig)
      deallocate(grd_tempinv)
      deallocate(grd_natom)
      deallocate(grd_nelec)
      deallocate(grd_vol)
      deallocate(grd_rho)
      deallocate(grd_mass)
      deallocate(grd_radtemp)
      deallocate(grd_lum)
      deallocate(grd_ye)
      deallocate(grd_rcell)
      deallocate(grd_divv)
      deallocate(grd_dvdx)
c
      ! deallocate(grd_capgam)
      deallocate(grd_emit)
      deallocate(grd_emitex)
      deallocate(grd_evolinit)
      deallocate(grd_jrad)
c-- ndim=3 integer
      deallocate(grd_nvol)
      deallocate(grd_methodswap)
c-- ndim=4 alloc
      deallocate(grd_opaclump)
      deallocate(grd_specarr,grd_nlump,grd_glumps,grd_llumps)
      deallocate(grd_emitlump,grd_caplump,grd_capemitlump,grd_doplump)
      deallocate(grd_capgreyinv,grd_capemitgreyinv)
      deallocate(grd_leaklcum,grd_leakrcum,grd_scatcum)
      deallocate(grd_dopcum,grd_dopidx,grd_ndop)
      deallocate(grd_emitprob)
c-- ndim=4 alloc
      deallocate(grd_cap)
      deallocate(grd_capemit)
      deallocate(grd_emiss)
      end subroutine grid_dealloc
c
      end module gridmod
c vim: fdm=marker
