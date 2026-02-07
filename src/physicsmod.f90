module physicsmod

use parametersmod, only : sp

implicit none

real(sp), parameter :: Tfreeze = 273.15

public :: pw
public :: lvap
public :: psygamma
public :: esat
public :: desdT
public :: Econ
public :: Dp
public :: pet
public :: dewpoint
public :: rhum
public :: abshum

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

real(sp) function psygamma(P,Tair)  ! (Pa K-1)

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

psygamma = Cp * P / (eps * lvap(Tair)) ! eqn 8

end function psygamma

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

! NB this function is differs between the Davis et al (2017) and Sandoval et al. (2024) formulations
! because of the addition of the 0.24 coefficient to gamma in the denominator in the newer publication.
! This constant allows the elimination of the entrainment factor omega in calculating potential evapotranspiration.
! Compare Davis eqns 19 and 22 with Sandoval eqn 51 and 50. See also
! Yang, Y., & Roderick, M. L. (2019). Radiation, surface temperature and evaporation over wet surfaces. 
! Quarterly Journal of the Royal Meteorological Society, 145(720), 1118-1129. doi:10.1002/qj.3481


use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: P     ! air pressure (Pa)
real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variables

real(sp) :: ss  ! rate of change of saturation vapor pressure with temperature (Pa K-1)

! ----

ss = desdT(Tair)

Econ = ss / (lvap(Tair) * pw(Tair) * (ss + 0.24 * psygamma(P,Tair)))  ! eqn 51

end function Econ

! ----------------------------------------------------------------------------------------------------------------

subroutine pet(P,Tair,HNpos,lw_rad,ru,rv,rw,dpet,dpmax)

! Estimation of daily total equilibrium evapotranspiration (dpet), mm
! and hourly maximum evapotranspiration rate (mm h-1)
! based on equations in Sandoval et al. and Davis et al. and source code

! Although Sandoval provides an equation for evaporative demand (Dp), it is not used directly in the calculation of  
! actual evapotranspiration. AET is instead calculated as a function of mean daytime shortwave and longwave fluxes,
! solar angle parameters, and evaporative supply rate (Sw). Sw is calculated in-turn as a function of the maximum 
! hourly evaporative demand (DpMAX), which is a function of solar angle variables and mean daytime longwave.

! 

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in)  :: P       ! mean air pressure (Pa)
real(sp), intent(in)  :: Tair    ! air temperature (degC)
real(sp), intent(in)  :: HNpos   ! daytime accumulated net radiation (J m-2 d-1)

real(sp), intent(in)  :: lw_rad  ! daily mean longwave radiation (W m-2)
real(sp), intent(in)  :: ru      ! simplification variables related to radiation
real(sp), intent(in)  :: rv      ! 
real(sp), intent(in)  :: rw      ! 

real(sp), intent(out) :: dpet    ! total daily potential evapotranspiration (mm d-1)
real(sp), intent(out) :: dpmax   ! maximum potential evapotranspiration rate (mm h-1)

! local variables

real(sp) :: Ec   ! energy-to-water conversion factor (m3 kJ-1)
real(sp) :: rx   ! simplification variable (mm m2 W-1 h-1)

! ----

Ec = Econ(P,Tair)   ! Sandoval eqn 51, see function  (m3 kJ-1)

dpet = Ec * HNpos   ! Sandoval eqn 50, but using HNpos so total per day, as per EVAP.cpp code (mm d-1)

rx = 3600. * Ec     ! m3 kJ-1 -> mm m2 W-1 h-1 [1000 (mm/m) / (1000 (J / kJ) / 3600 (s/h))]

dpmax = rx * ((rw * (ru + rv)) - lw_rad)  ! Sandoval eqn 53

end subroutine pet

! ----------------------------------------------------------------------------------------------------------------

real(sp) function dewpoint(tmin,tmax,dpet,pann)

! Estimate dewpoint temperature
! Based on Kimball, J. S., Running, S. W., & Nemani, R. (1997). 
! An improved method for estimating surface humidity from daily minimum temperature. 
! Agricultural and Forest Meteorology, 85(1), 87-98. doi:https://doi.org/10.1016/S0168-1923(96)02366-0

use parametersmod, only : sp,tfreeze

implicit none

! arguments

real(sp), intent(in) :: tmin  ! daily minimum temperature (degC)
real(sp), intent(in) :: tmax  ! daily maximum temperature (degC)
real(sp), intent(in) :: dpet  ! daily evapotranspiration (mm)
real(sp), intent(in) :: pann  ! total annual precipitation (mm)

! local variables

real(sp) :: tminK
real(sp) :: tmaxK
real(sp) :: EF
real(sp) :: tdewK

! ---

tminK = tmin + tfreeze
tmaxK = tmax + tfreeze

EF = dpet / pann

tdewK = tminK * (-0.127 + 1.121 * (1.003 - 1.444 * EF + 12.312 * EF**2 - 32.766 * EF**3) + 0.0006 * (tmaxK - tminK))  ! eqn 4

dewpoint = tdewK - tfreeze

end function dewpoint

! ----------------------------------------------------------------------------------------------------------------

real(sp) function rhum(tair,tdew)  ! (%)

! Estimate relative humidity from air temperature and dewpoint temperature
! Based on Lawrence, M. G. (2005). The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air: 
! A Simple Conversion and Applications. 
! Bulletin of the American Meteorological Society, 86(2), 225-234. doi:10.1175/bams-86-2-225

use parametersmod, only : sp,tfreeze,Rw

implicit none

! arguments

real(sp), intent(in)  :: tair  ! air temperature (degC)
real(sp), intent(in)  :: tdew  ! dewpoint temperature (degC)

! local variables

real(sp) :: Ta      ! air temperature (K)
real(sp) :: Td      ! dewpoint temperature (K)
real(sp) :: L       ! enthalpy of vaporization of water (J kg-1)

! ------

Ta = tair + tfreeze
Td = tdew + tfreeze

L = 1.91846e6 * (Ta / (Ta))**2  ! (J kg-1) Eqn. from Henderson-Sellers (1984)

rhum = 100. * exp(-L / (Rw * Ta * Td) * (Ta - Td)) ! eqn 12

end function rhum

! ----------------------------------------------------------------------------------------------------------------

real(sp) function abshum(tair,RH)

! Function to calculate absolute humidity, based on:
! Stull, R. (2017). Practical Meteorology: An Algebra-based Survey of Atmospheric Science. 
! Vancouver, Canada: Dept. of Earth, Ocean & Atmospheric Sciences, University of British Columbia.

use parametersmod, only : sp,tfreeze,Rw

implicit none

! arguments

real(sp), intent(in)  :: Tair  ! air temperature (degC)
real(sp), intent(in)  :: RH    ! relative humidity (%)

! local variables

real(sp) :: Ta      ! air temperature (K)

! ------

Ta = tair + tfreeze

Pvs = esat(Tair) / (Rw * Ta)

abshum = Pvs * RH / 100. 

end function abshum

! ----------------------------------------------------------------------------------------------------------------

end module physicsmod
