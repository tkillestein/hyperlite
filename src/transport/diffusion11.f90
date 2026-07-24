!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
pure subroutine diffusion11(ptcl,ptcl2,vx,vy,vz,rndstate,&
  eraddens,jrad,totevelo,ierr)
  use kindmod

  use randommod
  use miscmod
  use groupmod
  use gridmod
  use physconstmod
  use particlemod
  use transportmod
  use fluxmod
  use inputparmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real(dp),intent(inout) :: vx,vy,vz
  real(dp),intent(out) :: eraddens
  real(dp),intent(out) :: jrad
  real(dp),intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real(dp),parameter :: cinv = 1d0/pc_c
!
  integer :: l
  integer :: iig, iiig
  logical :: lhelp
  real(dp) :: r1, r2
  real(dp) :: denom
  real(dp) :: ddmct, tau, pa, pdop
!-- lumped quantities -----------------------------------------

  real(dp) :: emitlump, caplump, capemitlump, doplump
  real(dp) :: mfphelp, ppl, ppr
  real(dp) :: opacleak(2)
  real(dp) :: probleak(2)
  integer :: glump
  real(dp) :: capgreyinv
  real(dp) :: capemitgreyinv
  real(dp) :: speclump
  real(dp) :: dist, help
  real(dp) :: dummy

  integer,pointer :: ix, ic, ig
  integer,parameter :: iy=1, iz=1
  real(dp),pointer :: x, mu, e, e0, wl
  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!
  dummy = vy
  dummy = vz
  iiig = 0 !deterministic; misuse traps under -fcheck rather than UB

!-- no error by default
  ierr = 0
!-- init
  eraddens = 0d0
  jrad = 0d0
!
  dist = dx(ix)

!
!-- per-cell DDMC tables (precomputed once per sweep in ddmc_tables)
  capgreyinv = grd_capgreyinv(ic)
  capemitgreyinv = grd_capemitgreyinv(ic)
  speclump = grd_opaclump(7,ic)

!
!-- in lump?
  if(grd_cap(ig,ic)*dist >= trn_taulump) then
     glump = grd_nlump(ic)
  else
     glump = 0
  endif
!
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic))*dist < trn_tauddmc) then
     ierr = 100
     return
  endif
!
!-- retrieve from the per-cell tables
  if(glump>0) then
     emitlump = grd_emitlump(ic)
     caplump = grd_caplump(ic)
     capemitlump = grd_capemitlump(ic)
     doplump = grd_doplump(ic)
  else
!-- outside the lump
     if(in_nlte) then !NLTE
       emitlump = grd_specarr(ig,ic)* &
                  capemitgreyinv*grd_capemit(ig,ic)
     else ! LTE
       emitlump = grd_specarr(ig,ic)* &
                  capgreyinv*grd_cap(ig,ic)
     endif
     caplump = grd_cap(ig,ic)
     capemitlump = grd_capemit(ig,ic)
     if(ig/=grp_ng.and.grd_divv(ix).ge.0) then ! redshift
       doplump = dopspeccalc(grd_tempinv(ic),ig)*grd_divv(ix) &
                   /(grd_specarr(ig,ic)*3*pc_c)
     elseif(ig/=1.and.grd_divv(ix).lt.0) then ! blueshift
       doplump = -1*dopspeccalc(grd_tempinv(ic),ig-1)*grd_divv(ix) &
                   /(grd_specarr(ig,ic)*3*pc_c)
     else
       doplump = 0d0
     endif
  endif
!
!-- calculate lumped values
  if(glump>0) then
!-- leakage opacities
     opacleak = grd_opaclump(1:2,ic)
!-- calculating unlumped values
  else

!-- inward
     if(ix/=1) l = grd_icell(ix-1,iy,iz)
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,l)+ &
        grd_sig(l))*dx(ix-1)<trn_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= 1.5d0*ppl*(grd_xarr(ix))**2/dx3(ix)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))
        opacleak(1)=2.0d0*(grd_xarr(ix))**2/(mfphelp*dx3(ix))
     endif
!
!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*dx(ix+1)<trn_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=1.5d0*ppr*(grd_xarr(ix+1))**2/dx3(ix)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))
        opacleak(2)=2.0d0*(grd_xarr(ix+1))**2/(mfphelp*dx3(ix))
     endif
  endif
!
!-------------------------------------------------------------
!

!-- calculate time to event
  denom = sum(opacleak) + (1d0-emitlump)*caplump + doplump
  denom = 1d0/denom

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)*denom*cinv)
  ddmct = tau

!
!-- calculating energy depostion and density
  eraddens  = e*ddmct
!

!-- updating radiation intensity
  jrad = pc_c*eraddens


!-- perform event
  call rnd_r(r1,rndstate)

!-- leak probability
  probleak = opacleak*denom

!-- absorption probability
  pa = 0d0

!-- Doppler shift
  pdop = doplump*denom

!-- doppler shift
  if (r1>=pa .and. r1<pa+pdop) then

     if(glump==0) then
        iiig = ig
     else
!-- sample group (precomputed cumulative table; see ddmc_tables)
        call rnd_r(r1,rndstate)
        iig = firstgt(r1,grd_ndop(ic),grd_dopcum(:,ic))
        if(iig<=grd_ndop(ic)) then
           iiig = grd_dopidx(iig,ic)
        elseif(grd_divv(ix).ge.0) then ! scan fall-through: last iterated
           iiig = grd_glumps(glump,ic)
        else
           iiig = grd_glumps(1,ic)
        endif
     endif

!-- reshift/blueshift particle in this group
     if(grd_divv(ix).ge.0) then ! redshift
       ig = iiig+1
       wl = grp_wl(ig)
       ig = min(ig,grp_ng)
     else ! blueshift
       ig = iiig-1
       wl = grp_wl(ig+1)
       ig = max(ig,1)
     endif

!-- method changes to IMC
     if((grd_sig(ic)+grd_cap(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- direction sampled isotropically
        call rnd_r(r1,rndstate)
        mu = 1d0-2d0*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(ix+1)**3 + (1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- determine velocity coordinates (currently uses linear interpolation)
!-- x
        vx = grd_vxarr(ix) + grd_dvdx(ix)*(x - grd_xarr(ix))
!
!-- velocity effects accounting
        mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
        wl = wl*(1d0-vx*mu*cinv)
        help = 1d0/(1d0-vx*mu*cinv)
        totevelo = totevelo+e*(1d0 - help)
        e = e*help
        e0 = e0*help
     endif

!-- left leakage sample
  elseif (r1>=pa+pdop .and. r1<pa+pdop+probleak(1)) then
     ptcl2%idist = -3

!-- checking if at inner bound
     if(ix==1) then
!       stop 'diffusion11: non-physical inward leakage'
!        ierr = 101
        return

!-- sample adjacent group (assumes aligned ig bounds)
     else

        l = grd_icell(ix-1,iy,iz)
        if(glump==0) then
           iiig = ig
        else
!-- sample group (precomputed cumulative table; see ddmc_tables)
           call rnd_r(r1,rndstate)
           iig = min(firstgt(r1,glump,grd_leaklcum(:,ic)),glump)
           iiig = grd_glumps(iig,ic)
        endif
!
!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix-1) < trn_tauddmc) then
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set right bound of left cell
           x = grd_xarr(ix)
           vx = grd_vxarr(ix)
!-- particle angle sampled from isotropic b.c. inward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
!-- doppler and aberration corrections
           mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1d0-vx*mu*cinv)
        endif
!
!-- update particle
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
        ig = iiig

     endif


!-- right leakage sample
  elseif (r1>=pa+pdop+probleak(1) .and. r1<pa+pdop+sum(probleak)) then
     ptcl2%idist = -4
!
!-- checking if at outer bound
     if(ix==grd_nx) then
        ptcl2%stat = 'flux'
!-- outbound luminosity tally
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        if(glump==0) then
        else
!-- sample group (precomputed cumulative table; see ddmc_tables)
           call rnd_r(r1,rndstate)
           iig = min(firstgt(r1,glump,grd_leakrcum(:,ic)),glump)
           iiig = grd_glumps(iig,ic)
           ig = iiig
        endif
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(ig+1) + (1d0-r1)*grp_wlinv(ig))
!-- position
        x=grd_xarr(grd_nx+1)
        vx = grd_vxarr(grd_nx+1)
!-- changing from comoving frame to observer frame
        mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
        mu = min(mu,1d0)
        help = 1d0/(1d0-mu*vx*cinv)
!-- velocity effects accounting
        totevelo = totevelo+e*(1d0 - help)
        wl = wl/help
        e = e*help
        e0 = e0*help
        return
!
!
     else

        l = grd_icell(ix+1,iy,iz)
!-- sample adjacent group (assumes aligned ig bounds)
        if(glump==0) then
           iiig = ig
        else
!-- sample group (precomputed cumulative table; see ddmc_tables)
           call rnd_r(r1,rndstate)
           iig = min(firstgt(r1,glump,grd_leakrcum(:,ic)),glump)
           iiig = grd_glumps(iig,ic)
        endif

!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix+1)< trn_tauddmc) then
!
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set left bound of right cell
           x = grd_xarr(ix+1)
           vx = grd_vxarr(ix+1)
!-- particle angle sampled from isotropic b.c. outward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1d0-vx*mu*cinv)
        endif
!
!-- update particle
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
        ig = iiig
!
     endif!


!-- effective scattering sample
  else
     ptcl2%idist = -2

     if(glump==grp_ng) then
!       stop 'diffusion11: effective scattering with glump==ng'
        ierr = 102
        return
     endif

     if(glump==0) then
!-- sample group
        call rnd_r(r1,rndstate)
        iiig = emitgroup(r1,ic)
     else
!-- sample group (precomputed cumulative table; see ddmc_tables)
        call rnd_r(r1,rndstate)
        iig = grp_ng + 1 - &
             min(firstgt(r1,grp_ng-glump,grd_scatcum(:,ic)),grp_ng-glump)
        iiig = grd_glumps(iig,ic)
     endif
     ig = iiig
     if((grd_sig(ic)+grd_cap(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
!-- direction sampled isotropically
        call rnd_r(r1,rndstate)
        mu = 1d0-2d0*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(ix+1)**3 + (1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- determine velocity coordinates (currently uses linear interpolation)
!-- x
        vx = grd_vxarr(ix) + grd_dvdx(ix)*(x - grd_xarr(ix))
!
!-- doppler and aberration corrections
        mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
!-- velocity effects accounting
        help = 1d0/(1d0-vx*mu*cinv)
        totevelo = totevelo+e*(1d0 - help)
!
        e = e*help
        e0 = e0*help
        wl = wl*(1d0-vx*mu*cinv)
     endif

  endif

!-- eliminate particles crossing inner boundary
  if (x  < grd_xarr(1)) then
    ptcl2%stat = 'dead'
    return
  endif

!-- Doppler shift energy weights and wavelengths ! added - gaw
  help = exp(-grd_divv(ix)*ddmct/3)
  totevelo = totevelo + e*(1d0-help)
  e = e*help
  e0 = e0*help


contains

  pure real(dp) function dx(l)
    use kindmod
    integer, intent(in) :: l
    dx = grd_xarr(l+1) - grd_xarr(l)
  end function dx

  pure real(dp) function dx3(l)
    use kindmod
    integer, intent(in) :: l
    dx3 = grd_xarr(l+1)**3 - grd_xarr(l)**3
  end function dx3

  pure integer function firstgt(r,n,cum)
!-- first k in [1,n] with cum(k)>r, or n+1 if none: the index the old
!-- sequential scan exits at.  cum is nondecreasing except for a
!-- possible NaN tail from degenerate (Inf) per-cell tables; a NaN last
!-- entry falls back to the sequential scan, which reproduces the
!-- pre-table comparisons exactly.
    use kindmod
    real(dp), intent(in) :: r
    integer, intent(in) :: n
    real(dp), intent(in) :: cum(:)
    integer :: klo, khi, kmid
    if(n<1) then
       firstgt = n+1
       return
    endif
    if(cum(n)/=cum(n)) then
       do kmid=1,n
          if(cum(kmid)>r) exit
       enddo
       firstgt = kmid
       return
    endif
    klo = 1
    khi = n+1
    do while(klo<khi)
       kmid = (klo+khi)/2
       if(cum(kmid)>r) then
          khi = kmid
       else
          klo = kmid+1
       endif
    enddo
    firstgt = klo
  end function firstgt

end subroutine diffusion11
! vim: fdm=marker
