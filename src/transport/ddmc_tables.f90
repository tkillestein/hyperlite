!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine ddmc_tables
  use kindmod
  use gridmod
  use groupmod
  use transportmod
  use physconstmod
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
!
! Also precompute the cumulative group-sampling tables (left/right
! leakage, effective scattering, doppler shift): the exact partial sums
! that the diffusion11 linear scans accumulate per event, in scan order.
! An upper-bound search on these returns the same group the scan exits
! at, with the same draw sequence -> byte-identical results.
!##################################################
  integer :: i,j,k,l,iig,glump,gunlump
  real(dp) :: dist
  real(dp) :: emitlump,caplump,capemitlump,doplump,speclump
!-- cumulative sampling tables
  integer :: iiig,lnb,ndop
  logical :: lhelp
  real(dp) :: dist3,dxnb,help,cum,specig,mfphelp,pp
  real(dp) :: resopacleak,resdopleak

!$omp parallel do schedule(static) default(shared) &
!$omp    private(i,j,k,l,iig,glump,gunlump,dist, &
!$omp    emitlump,caplump,capemitlump,doplump,speclump, &
!$omp    iiig,lnb,ndop,lhelp,dist3,dxnb,help,cum,specig,mfphelp,pp, &
!$omp    resopacleak,resdopleak)
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
!
!-- cumulative group-sampling tables ---------------------------
!-- (expressions replicate the diffusion11 scans verbatim)
     if(glump>0) then
        dist3 = grd_xarr(i+1)**3 - grd_xarr(i)**3
!
!-- left leakage (unused at i==1: diffusion11 returns before sampling)
        if(i>1 .and. grd_opaclump(1,l)>0d0) then
           lnb = grd_icell(i-1,j,k)
           dxnb = grd_xarr(i) - grd_xarr(i-1)
           help = 1d0/grd_opaclump(1,l)
           cum = 0d0
           do iig=1,glump
              iiig = grd_glumps(iig,l)
              specig = grd_specarr(iiig,l)
              if((grd_cap(iiig,lnb)+ &
                   grd_sig(lnb))*dxnb<trn_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,l)+grd_sig(l))*dist
                 pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*pp*(grd_xarr(i))**2/dist3
              else
!-- IMC interface
                 mfphelp = ((grd_sig(l)+grd_cap(iiig,l))*dist+ &
                      (grd_sig(lnb)+grd_cap(iiig,lnb))*dxnb)
                 resopacleak = 2.0d0*(grd_xarr(i))**2/ &
                      (mfphelp*dist3)
              endif
              cum = cum + specig*resopacleak*speclump*help
              grd_leaklcum(iig,l) = cum
           enddo
        else
           grd_leaklcum(1:glump,l) = 0d0
        endif
!
!-- right leakage (i==grd_nx uses the DDMC-interface formula, as at the
!-- outer-boundary call site)
        if(grd_opaclump(2,l)>0d0) then
           lhelp = i==grd_nx
           if(lhelp) then
              lnb = l
              dxnb = dist
           else
              lnb = grd_icell(i+1,j,k)
              dxnb = grd_xarr(i+2) - grd_xarr(i+1)
           endif
           help = 1d0/grd_opaclump(2,l)
           cum = 0d0
           do iig=1,glump
              iiig = grd_glumps(iig,l)
              specig = grd_specarr(iiig,l)
              if(lhelp .or. (grd_cap(iiig,lnb)+ &
                   grd_sig(lnb))*dxnb<trn_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,l)+grd_sig(l))*dist
                 pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*pp*(grd_xarr(i+1))**2/dist3
              else
!-- IMC interface
                 mfphelp = ((grd_sig(l)+grd_cap(iiig,l))*dist+ &
                      (grd_sig(lnb)+grd_cap(iiig,lnb))*dxnb)
                 resopacleak = 2.0d0*(grd_xarr(i+1))**2/ &
                      (mfphelp*dist3)
              endif
              cum = cum + specig*resopacleak*speclump*help
              grd_leakrcum(iig,l) = cum
           enddo
        else
           grd_leakrcum(1:glump,l) = 0d0
        endif
!
!-- effective-scattering resample over the unlumped groups, scan order
!-- iig=grp_ng,glump+1,-1; entry j maps to position grp_ng-j+1
        if(glump<grp_ng) then
           help = 1d0/(1d0-emitlump)
           cum = 0d0
           do iig=grp_ng,glump+1,-1
              iiig = grd_glumps(iig,l)
              if(in_nlte) then !NLTE
                 cum = cum + grd_specarr(iiig,l)*grd_capemit(iiig,l)* &
                      grd_capemitgreyinv(l)*help
              else
                 cum = cum + grd_specarr(iiig,l)*grd_cap(iiig,l)* &
                      grd_capgreyinv(l)*help
              endif
              grd_scatcum(grp_ng-iig+1,l) = cum
           enddo
        endif
!
!-- doppler shift, compacted to the entries the scan tests (skipped
!-- entries never test the cumulative sum); grd_divv indexed by the
!-- x position as in diffusion11
        ndop = 0
        if(doplump>0d0) then
           help = 1d0/doplump
           cum = 0d0
           if(grd_divv(i).ge.0) then ! redshift
              do iig=1,glump
                 iiig = grd_glumps(iig,l)
                 if(iiig == grp_ng) cycle
                 if(grd_cap(iiig+1,l)*dist >= trn_taulump) cycle
                 resdopleak = dopspeccalc(grd_tempinv(l),iiig)* &
                      grd_divv(i)/(3*pc_c)
                 cum = cum + resdopleak*speclump*help
                 ndop = ndop + 1
                 grd_dopidx(ndop,l) = int(iiig,2)
                 grd_dopcum(ndop,l) = cum
              enddo
           else ! blueshift
              do iig=glump,1,-1
                 iiig = grd_glumps(iig,l)
                 if(iiig == 1) cycle
                 if(grd_cap(iiig-1,l)*dist >= trn_taulump) cycle
                 resdopleak = dopspeccalc(grd_tempinv(l),iiig-1)* &
                      grd_divv(i)/(3*pc_c)
                 cum = cum - resdopleak*speclump*help
                 ndop = ndop + 1
                 grd_dopidx(ndop,l) = int(iiig,2)
                 grd_dopcum(ndop,l) = cum
              enddo
           endif
        endif
        grd_ndop(l) = ndop
     endif !glump>0
  enddo
  enddo
  enddo
!$omp end parallel do

end subroutine ddmc_tables
! vim: fdm=marker
