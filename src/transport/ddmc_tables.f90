!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine ddmc_tables
  use kindmod
  use gridmod
  use groupmod
  use transportmod
  use inputparmod, only:in_nlte
  implicit none
!##################################################
! Precompute the per-cell DDMC quantities that are constant during one
! transport sweep: Planck spectral weights (specintv), the group-lump
! partition, and the lumped opacities.  Replaces the per-thread
! grp_t_cache that was rebuilt on nearly every cell change inside
! diffusion11 (~one specintv call per DDMC event; see OVERHAUL Phase 6).
! Values are identical to the old per-particle rebuild: same
! expressions, evaluated once per cell.
!##################################################
  integer :: i,j,k,l,iig,glump,gunlump
  real(dp) :: dist
  real(dp) :: emitlump,caplump,capemitlump,doplump,speclump

!$omp parallel do schedule(static) default(shared) &
!$omp    private(i,j,k,l,iig,glump,gunlump,dist, &
!$omp    emitlump,caplump,capemitlump,doplump,speclump)
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     dist = grd_xarr(i+1) - grd_xarr(i)
!
     grd_capgreyinv(l) = max(1d0/grd_capgrey(l),0d0) !catch nans
     grd_capemitgreyinv(l) = max(1d0/grd_capemitgrey(l),0d0)
!
!-- Planck spectral weights
     call specintv(grd_tempinv(l),grp_ng,grd_specarr(:,l))
!
!-- lump testing ---------------------------------------------
     glump = 0
     gunlump = grp_ng
     grd_glumps(:,l) = 0
!
!-- find lumpable groups
     speclump = grd_opaclump(7,l)
     if(speclump==0d0) then
        glump=0
        grd_llumps(:,l) = .false.
     else
        do iig=1,grp_ng
           if(grd_cap(iig,l)*dist >= trn_taulump .and. &
                (grd_sig(l) + grd_cap(iig,l))*dist >= trn_tauddmc) then
              grd_llumps(iig,l) = .true.
              glump=glump+1
              grd_glumps(glump,l) = int(iig,2)
           else
              grd_llumps(iig,l) = .false.
              grd_glumps(gunlump,l) = int(iig,2)
              gunlump=gunlump-1
           endif
        enddo
     endif
!
!-- calculate lumped values
     if(glump==grp_ng) then
        emitlump = 1d0
        caplump = grd_capgrey(l)
        capemitlump = grd_capemitgrey(l)
        doplump = 0d0
     else
!-- Planck x-section lump
        caplump = grd_opaclump(8,l)*speclump
        capemitlump = grd_opaclump(9,l)*speclump
        if(in_nlte) then
          emitlump = grd_opaclump(9,l)*grd_capemitgreyinv(l)
        else
          emitlump = grd_opaclump(8,l)*grd_capgreyinv(l)
        endif
        emitlump = min(emitlump,1d0)
        doplump = grd_opaclump(11,l)*speclump
     endif
!
     grd_nlump(l) = glump
     grd_emitlump(l) = emitlump
     grd_caplump(l) = caplump
     grd_capemitlump(l) = capemitlump
     grd_doplump(l) = doplump
  enddo
  enddo
  enddo
!$omp end parallel do

end subroutine ddmc_tables
! vim: fdm=marker
