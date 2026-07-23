*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module transportmod
      use kindmod
c     -------------------
      use physconstmod, only:pc_c,pc_pi2
      implicit none
c
      private pc_c,pc_pi2
      real(dp),private,parameter :: cinv=1d0/pc_c
c
      logical :: trn_errorfatal     !stop on transport error, disable for production runs
      logical :: trn_noampfact      !don't use the particle amplification factor
c
      real(dp) :: trn_tauddmc
      real(dp) :: trn_taulump
      logical :: trn_thm_aniso !experimental !evaluates anisotropic scattering angle for Thomson scattering
c
c-- explicit interfaces
      interface
c
c-- transport
      pure subroutine transport11(ptcl,ptcl2,vx,vy,vz,
     &  rndstate,eraddens,eamp,jrad,totevelo,ierr)
      use kindmod
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      real(dp),intent(inout) :: vx,vy,vz
      type(rnd_t),intent(inout) :: rndstate
      real(dp),intent(out) :: eraddens, eamp
      real(dp),intent(out) :: jrad
      real(dp),intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport11
c
c-- diffusion
c
      pure subroutine diffusion11(ptcl,ptcl2,vx,vy,vz,
     &  rndstate,eraddens,jrad,totevelo,ierr)
      use kindmod
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      real(dp),intent(inout) :: vx,vy,vz
      type(rnd_t),intent(inout) :: rndstate
      real(dp),intent(out) :: eraddens
      real(dp),intent(out) :: jrad
      real(dp),intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion11

      end interface
c
c-- abstract interfaces
      abstract interface
      subroutine direction2lab_(x0,y0,z0,mu0,om0)
      use kindmod
c     -------------------------------------------
      real(dp),intent(in) :: x0,y0,z0
      real(dp),intent(inout) :: mu0,om0
      end subroutine direction2lab_
c
      pure subroutine transport_(ptcl,ptcl2,vx,vy,vz,
     &  rndstate,eraddens,eamp,jrad,totevelo,ierr)
      use kindmod
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      real(dp),intent(inout) :: vx,vy,vz
      type(rnd_t),intent(inout) :: rndstate
      real(dp),intent(out) :: eraddens, eamp
      real(dp),intent(out) :: jrad
      real(dp),intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport_
c
      pure subroutine diffusion_(ptcl,ptcl2,vx,vy,vz,
     &  rndstate,eraddens,jrad,totevelo,ierr)
      use kindmod
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      real(dp),intent(inout) :: vx,vy,vz
      type(rnd_t),intent(inout) :: rndstate
      real(dp),intent(out) :: eraddens
      real(dp),intent(out) :: jrad
      real(dp),intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion_

      end interface
c
c-- procedure pointers
      procedure(direction2lab_),pointer :: direction2lab => null()
      procedure(transport_),pointer :: transport => null()
      procedure(diffusion_),pointer :: diffusion => null()
c
c-- private interfaces
      private diffusion11
      private transport11
c
      save
c
      contains
c
c
c
      subroutine transportmod_init(igeom)
c     -----------------------------------------------
      integer,intent(in) :: igeom
c
c-- set procedure pointers
      select case(igeom)
      case(11)
       direction2lab => direction2lab1
       transport => transport11
       diffusion => diffusion11
      case default
       stop 'transportmod_init: invalid igeom'
      end select
      end subroutine transportmod_init
c
c
c
      subroutine direction2lab1(vx0,vy0,vz0,mu0,om0)
c     -------------------------------------------
      implicit none
      real(dp),intent(in) :: vx0,vy0,vz0
      real(dp),intent(inout) :: mu0,om0
      real(dp) :: cmffact,mu,dummy
c
      dummy = vy0
      dummy = vz0
      dummy = om0
c
      cmffact = 1d0+mu0*vx0*cinv
      mu = (mu0+vx0*cinv)/cmffact
      mu = min(mu,1d0)
      mu = max(mu,-1d0)
      mu0 = mu
      end subroutine direction2lab1
c
      end module transportmod
c vim: fdm=marker
