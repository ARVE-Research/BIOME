module soilphysicsmod

implicit none

public :: Ksat_env
public :: infiltration

contains

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Ksat_env(Ksat,P,Tair)

use parametersmod, only : sp,g
use physicsmod,    only : pwater,muwater

implicit none

! Function to adjust the saturated conductivity of soil to ambient environmental conditions
! based on temperature and pressure

! arguments

real(sp), intent(in) :: Ksat  ! saturated hydraulic conductivity under standard conditions (mm h-1)
real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)

! parameters

real(sp), parameter :: mustd = 8.87255592e-4  ! viscosity of water at SATP
real(sp), parameter :: pstd  = 997.045532     ! density of water at SATP (kg m-3)

! local variables

real(sp) :: ki

! ----

ki = Ksat * mustd / (pstd * g)

Ksat_env = ki * pwater(P,Tair) * g / muwater(P,Tair)

end function Ksat_env

! ----------------------------------------------------------------------------------------------------------------

subroutine infiltration(srad,dmet,soilstate)

! infiltration following Sandoval et al. (2024) section 2.2.4

use parametersmod, only : sp
use typesmod,      only : metvars_daily,soilstatetype

implicit none

! arguments

real(sp),            intent(in)    :: srad       ! terrain slope (radians)
type(metvars_daily), intent(in)    :: dmet
type(soilstatetype), intent(inout) :: soilstate  ! state variables of the soil (top layer only)


! local variables

real(sp) :: rain    ! liquid precipitation (mm d-1)
real(sp) :: melt    ! snowmelt (mm d-1)

real(sp) :: Tsat    ! soil porosity (fraction)
real(sp) :: Ksat    ! saturated hydraulic conductivity (mm h-1)
real(sp) :: lambda  ! pore size distribution (unitless)
real(sp) :: psi_e   ! soil water potential at air entry (mm)

real(sp) :: theta  ! volumetric water content (fraction)

real(sp) :: psi_f   ! capillary head at the wetting front (mm)
real(sp) :: tp      ! ponding time (?)

real(sp) :: r       ! total water flux at the surface (mm h-1)

real(sp) :: cos2s   ! cos2(slope)

! ----

rain = dmet%rain
melt = dmet%melt

Tsat = soilstate%Tsat
theta = soilstate%theta

! ----

! ki = Ksat0 * mustd / pstd
! 
! Ksat = ki * pg / muw()
! 
! cos2s = cos(srad)**2
! 
! r = (rain + melt) / 24.  ! convert to mm h-1
! 
! psi_f = (2. + 3. * lambda) / (1. + 3 * lambda) * psi_e / 2.
! 
! tp = Ksat * psi_f * (Tsat - theta) / (r * (r - Ksat))
! 
! rtd = 
! 
! if (r <= Ksat) then
! 
!   infil = r * 

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
