module physicsmod

use parametersmod, only : sp

implicit none

real(sp), parameter :: Tfreeze = 273.15

public :: pw
public :: lvap
public :: gamma
public :: esat
public :: desdT
public :: Econ
public :: Dp
public :: PET

contains

! ----------------------------------------------------------------------------------------------------------------

real(sp) function stdP(elv)  ! (Pa)

! Function to calculate air pressure typical of the standard atmosphere, based on:
! Stull, R. (2017). Practical Meteorology: An Algebra-based Survey of Atmospheric Science. 
! Vancouver, Canada: Dept. of Earth, Ocean & Atmospheric Sciences, University of British Columbia.

! valid for elevations below 11 km a.s.l.

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: elv  ! elevation (meters above sea level)

! local variable

real(sp) :: H  ! elevation (km)
real(sp) :: T  ! standard atmosphere temperature

! ----

H = elv / 1000.

T = 288.15 - 6.5 * H  ! eqn 1.16

stdP = 101325. * (288.15 / T)**(-5.255877)  ! eqn 1.17

end function stdP

! ----------------------------------------------------------------------------------------------------------------

real(sp) function pw(T)  ! (kg m-3)

! Function to calculate the temperature-dependent density of liquid water, based on:
! Kell, G. S. (1975). Density, Thermal Expansivity, and Compressibility of Liquid Water from 0° to 150°C: 
! Correlations and Tables for Atmospheric Pressure and Saturation Reviewed and Expressed on 1968 Temperature Scale. 
! Journal of Chemical and Engineering Data, 20(1), 97-105. 

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: T  ! water temperature (degC)

! parameters

real(sp), parameter :: a = 999.83952
real(sp), parameter :: b =  16.945176
real(sp), parameter :: c =   7.9870401e-3
real(sp), parameter :: d =  46.170461e-6
real(sp), parameter :: e = 105.56302e-9
real(sp), parameter :: f = 280.54253e-12
real(sp), parameter :: g =  16.879850e-3

! ----

pw = (a + b * T - c * T**2 - d * T**3 + e * T**4 + f * T**5) / (1. + g * T)  ! Kell (1975) eqn 16

end function pw

! ----------------------------------------------------------------------------------------------------------------

real(sp) function lvap(Tair)  ! (kJ kg-1)

! function to calculate the temperature-dependent latent heat of vaporization of water, based on:
! Henderson-Sellers, B. (1984). A new formula for latent heat of vaporization of water as a function of temperature.
! Quarterly Journal of the Royal Meteorological Society, 110(466), 1186-1190. doi:10.1002/qj.49711046626

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variable

real(sp) :: T  ! air temperature (K)

! ----

T = Tair + Tfreeze

lvap = 0.001 * 1.91846e6 * (T / (T - 33.91))**2  ! eqn 8

end function lvap

! ----------------------------------------------------------------------------------------------------------------

real(sp) function gamma(P,Tair)  ! (Pa K-1)

! function to calculate the psychrometric constant
! there are many sources for the formula to calculate temperature- and pressure-dependent gamma
! one of the first is in the reference below, although this paper does not provide a thermodynamic rationale
! Tanner, C. B., & Fuchs, M. (1968). Evaporation from unsaturated surfaces: A generalized combination method. 
! Journal of Geophysical Research, 73(4), 1299-1304. doi:10.1029/JB073i004p01299
! constants from AMS Glossary and:
! Gatley, D. P., Herrmann, S., & Kretzschmar, H.-J. (2011). A Twenty-First Century Molar Mass for Dry Air. 
! HVAC&R Research, 14(5), 655-662. doi:10.1080/10789669.2008.10391032

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)

real(sp), parameter  :: eps = 0.6219453  ! ratio of molar mass of water vapor to dry air (Gatley et al., 2011)
real(sp), parameter  :: Cp  = 1.0057     ! specific heat capacity of air at constant pressure (kJ kg-1 K-1) (AMS Glossary)

! ----

gamma = Cp * P / (eps * lvap(Tair)) ! eqn 8

end function gamma

! ----------------------------------------------------------------------------------------------------------------

real(sp) function esat(Tair)   ! (Pa)

! Function to calculate saturation vapor pressure in water and ice, based on:
! Ambaum, M. H. P. (2020). Accurate, simple equation for saturated vapour pressure over water and ice. 
! Quarterly Journal of the Royal Meteorological Society, 146(733), 4252-4258. doi:10.1002/qj.3899

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: Tair  ! temperature (degC)

! parameters

real(sp), parameter :: Rv  =  461.52  ! specific gas constant for water vapor (J kg-1 K-1)
real(sp), parameter :: e0  =  611.655 ! saturation vapor pressure at the triple point of water (Pa)
real(sp), parameter :: T0  =  273.16  ! temperature of water at its triple point (K)

real(sp), parameter :: cl  = 2180.    ! difference in isobaric heat capacity: liquid - vapor (cpl - cpv) (J kg-1 K-1)
real(sp), parameter :: ci  =  212.    ! difference in isobaric heat capacity: ice - vapor (cpi - cpv) (J kg-1 K-1)

real(sp), parameter :: Lv0 = 2.501e6  ! latent heat of vaporization at the triple point of water (J kg-1)
real(sp), parameter :: Ls0 = 2.834e6  ! latent heat of sublimation at the triple point of water (J kg-1)

real(sp), parameter :: clRv = cl / Rv ! combination of above terms for liquid water
real(sp), parameter :: ciRv = ci / Rv ! combination of above terms for ice

real(sp), parameter :: LRTv = Lv0 / (Rv * T0) ! combination of above terms for liquid water
real(sp), parameter :: LRTs = Ls0 / (Rv * T0) ! combination of above terms for ice

! variables

real(sp) :: T  ! air temperature (K)
real(sp) :: L  ! latent heat of vaporization of water (temperature dependent) (J kg-1)

! ----

T = Tair + Tfreeze

if (Tair >= 0.) then  ! equation for liquid water

  L = Lv0 - cl * (T - T0)
  
  esat = e0 * (T0 / T)**clRv * exp(LRTv - L / (Rv * T))  ! eqn 13
  
else  ! equation for ice

  L = Ls0 - ci * (T - T0)
  
  esat = e0 * (T0 / T)**ciRv * exp(LRTs - L / (Rv * T))  ! eqn 17
  
end if

end function esat

! ----------------------------------------------------------------------------------------------------------------

real(sp) function desdT(Tair)  ! (Pa K-1)

! Function to calculate the rate of chance of saturation vapor pressure as a function of temperature, based on:
! Yang, Y., & Roderick, M. L. (2019). Radiation, surface temperature and evaporation over wet surfaces.
! Quarterly Journal of the Royal Meteorological Society, 145(720), 1118-1129. doi:10.1002/qj.3481
! the origin of this function appears to be from Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
! Crop evapotranspiration - Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56. 
! Rome: FAO - Food and Agriculture Organization of the United Nations.
! comparison with other algorithms yields similar results
! JOK Jan 2025

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: Tair  ! temperature (degC)

! variable

real(sp) :: T  ! air temperature (K)

! ----

T = Tair + Tfreeze

desdT = 4098. * esat(Tair) / (T - 35.8)**2  ! eqn 3

end function desdT

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Econ(P,Tair)  ! (m3 kJ-1)

! Function to calculate the energy-to-water conversion factor, based on:
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024
! because lvap is calculated in kJ, the result of this equation is kJ

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variables

real(sp) :: ss  ! rate of change of saturation vapor pressure with temperature (Pa K-1)

! ----

ss = desdT(Tair)

Econ = ss / (lvap(Tair) * pw(Tair) * (ss + 0.24 * gamma(P,Tair)))  ! eqn 51

end function Econ

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Dp(P,Tair,Inet)  ! (mm h-1)

! Function to calculated hourly evaporative demand, based on:
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)
real(sp), intent(in) :: Inet  ! net radiation (W m-2)

Dp = 1000. * 3600. * Econ(P,Tair) * Inet  ! eqn 50

end function Dp

! ----------------------------------------------------------------------------------------------------------------

real(sp) function PET(P,Tair,netrad)  ! (mm)

! Function to calculate integrated daily potential evapotranspiration, based on:
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: P       ! air pressure (Pa)
real(sp), intent(in) :: Tair    ! air temperature (degC)
real(sp), intent(in) :: netrad  ! net radiation (kJ m-2 d-1)

! parameter

real(sp), parameter :: omega = 0.26  ! entrainment factor, dimensionless (Sandoval et al., 2024; Priestley and Taylor, 1972)

! ----

pet = (1. + omega) * 1000. * Econ(P,Tair) * max(netrad,0.)  ! 1000 converts m to mm

end function PET

! ----------------------------------------------------------------------------------------------------------------

end module physicsmod
