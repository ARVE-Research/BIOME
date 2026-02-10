module soilphysicsmod

implicit none

public :: Ksat_env
public :: infiltration

contains

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Ksat_env(ki,P,Tair)  ! (mm h-1)

use parametersmod, only : sp,g
use physicsmod,    only : pwater,muwater

implicit none

! Function to calculate instantaneous saturated conductivity as influenced by ambient environmental conditions
! (temperature and pressure)

! arguments

real(sp), intent(in) :: ki    ! soil intrinsic permeability (m2)
real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)

! ----
! NB the 3.6e6 factor converts Ksat (m s-2) to Ksat (mm h-1)

Ksat_env = 3.6e6 * ki * pwater(P,Tair) * g / muwater(P,Tair)

end function Ksat_env

! ----------------------------------------------------------------------------------------------------------------

subroutine infiltration(pixel,dmet,soilstate)

! infiltration following Sandoval et al. (2024) section 2.2.4

use parametersmod, only : sp
use typesmod,      only : pixeltype,metvars_daily,soilstatetype

implicit none

! arguments

type(pixeltype),     intent(in)    :: pixel
type(metvars_daily), intent(in)    :: dmet
type(soilstatetype), intent(inout) :: soilstate  ! state variables of the soil (top layer only)

! parameters

real(sp), parameter :: td = 6.   ! precipitation duration (hr)

! local variables

real(sp) :: P       ! atmospheric pressure (Pa)
real(sp) :: srad    ! terrain slope (radians)

real(sp) :: rain    ! liquid precipitation (mm d-1)
real(sp) :: melt    ! snowmelt (mm d-1)

real(sp) :: Tair

real(sp) :: Tsat    ! soil porosity (fraction)
real(sp) :: Ksat    ! saturated hydraulic conductivity (mm h-1)
real(sp) :: ki      ! intrinsic permeability (m2)
real(sp) :: lambda  ! pore size distribution (unitless)
real(sp) :: psi_e   ! soil water potential at air entry (mm)

real(sp) :: theta  ! volumetric water content (fraction)

real(sp) :: psi_f   ! capillary head at the wetting front (mm)
real(sp) :: tp      ! ponding time (?)

real(sp) :: r       ! total water flux at the surface (mm h-1)

real(sp) :: cos2s   ! cos2(slope)
real(sp) :: tp1     ! adjusted ponding time = ponding time / cos2(slope), shorter time under steeper slope

real(sp) :: infil   ! soil water infiltration flux (mm h-1)

real(sp) :: dtheta

! ----

P    = pixel%P
srad = pixel%srad

Tair = (dmet%tmax + dmet%tmin) / 2.  ! estimate of 24-hr mean temperature
rain = dmet%rain
melt = dmet%melt

ki    = soilstate%ki
Tsat  = soilstate%Tsat
theta = soilstate%theta
lambda = soilstate%lambda
psi_e = soilstate%psi_e

! ----

dtheta = 0.1

r = (rain + melt) / td  ! convert daily total to mm h-1

Ksat = Ksat_env(ki,P,Tair)  ! calculate current saturated conductivity based on temperature and pressure

if (r <= Ksat) then

  infil = r * td   ! eqn 31b
  
  ! write(0,*)'infil a',rain,melt,r,Ksat,infil
  
else

  cos2s = cos(srad)**2

  psi_f = (2. + 3. * lambda) / (1. + 3. * lambda) * psi_e / 2.  ! eqn 29
  
  tp = Ksat * psi_f * (dtheta) / (r * (r - Ksat))        ! eqn 30

  tp1 = tp / cos2s

  write(0,*)'infil b',rain,melt,r,Ksat,srad,cos2s,psi_e,psi_f,tp,tp1,r * tp1 / (psi_f * dtheta)

  infil = r * tp1 + Ksat * (td - tp1) - psi_f * dtheta * log(1. - r * tp1 / (psi_f * dtheta))

end if

end subroutine infiltration

! ----------------------------------------------------------------------------------------------------------------

! subroutine soilwaterflux(dmet,soilstate)
! 
! real(sp), allocatable, dimension(:) :: zposmm  ! depth of soil layer midpoint, zero = surface, positive down (mm)
! 
! real(sp), allocatable, dimension(:) :: psi     ! soil matric potential at the layer midpoint (mm)
! real(sp), allocatable, dimension(:) :: theta
! 
! 
! real(dp), dimension(:) :: Ku           !soil water instantaneous (unsaturated) conductivity at layer boundary (mm s-1)
! 
! 
! ! timestepping stage
! ! set timestep dynamically based on infiltration at the soil surface
! 
! 
!   ! -------------------------------------------------
!   ! Set up r, a, b, and c vectors for tridiagonal solution of the Richards equation
! 
!   ! ---------
!   ! top layer
! 
!   l = 1  ! layer index number
! 
!   Fin = qinfl  ! water flux at the soil surface
! 
!   dz = zposmm(l+1) - zposmm(l)  ! distance between soil layer midpoints (mm)
! 
!   dpsi_eq = psi_eq(l+1) - psi_eq(l)  ! gradient in 
! 
!   dpsi = (psi(l+1) - psi(l)) - dpsi_eq  ! pressure gradient minus equilibrium water content
! 
!   Fout = -Ku(l) * nterm / dz
! 
!   ddFoutTliq1 = -(-Ku(l) * ddpsiTliq(l) + nterm * ddKuTliq(l)) / dterm
! 
!   ddFoutTliq2 = -( Ku(l) * ddpsiTliq(l+1) + nterm * ddKuTliq(l)) / dterm
! 
!   rvect(l) =  Fin - Fout - rootfeff(l)
! 
!   avect(l) =  0._dp
! 
!   bvect(l) =  dzmm(l) * (1._dp/dt) + ddFoutTliq1
! 
!   cvect(l) =  ddFoutTliq2
!   
!   
! 
! end subroutine soilwaterflux

end module soilphysicsmod
