module radiationmod

! consolidated module for calculation for surface downwelling shortwave, net longwave, net radiation and PET.

implicit none

! module subroutines and functions

! public  :: radpet
! private :: dewpoint
! private :: surf_sw
! private :: surf_lw
! private :: netrad_pet

contains

! ----------------------------------------------------------------------------------------------------------------

subroutine radpet(pixel,dmet)

! shortwave and longwave radiation and dew point all depend on potential evapotranspiration (PET)
! we therefore pick an initial value for PET in this routine and iterate to a stable solution

use parametersmod, only : sp,dp,pi => pi_sp,pir => pir_sp
use typesmod,      only : orbitpars,airmasspars,pixeltype,solarpars,metvars_daily
use airmassmod,    only : airmass

implicit none

! arguments

type(pixeltype),     intent(in)    :: pixel
! type(solarpars),     intent(in)    :: solar
type(metvars_daily), intent(inout) :: dmet    ! daily meteorological variables

! parameter

real(sp) :: pisec = 86400. / pi

! local variables

real(dp) :: lat     ! latitude (degrees)
real(sp) :: elv     ! elevation (m)
real(sp) :: tday    ! daytime mean temperature (C)
real(sp) :: tnight  ! nighttime mean temperature (C)
real(sp) :: prec    ! daily total precipitation
real(sp) :: cldf    ! cloud cover fraction
real(sp) :: tdew    ! estimated dewpoint temperature (degC)
  
real(sp) :: tcm     ! temperature of the coldest month
real(sp) :: Pann    ! total annual precipitation (mm) 
real(sp) :: Pjj     ! precipitation equitability index 
real(sp) :: Ratm    ! relative atmospheric pressure (based on elevation)
real(sp) :: P       ! mean atmospheric pressure (based on elevation)

real(sp) :: albedo 

type(airmasspars) :: air

real(sp) :: toa_sw  ! mean daily total top-of-the-atmosphere downwelling shortwave rad (W m-2)
real(sp) :: delta   ! solar declination (degrees)
real(sp) :: direct  ! mean daily direct beam surface downwelling shortwave (W m-2)
real(sp) :: diffuse ! mean daily diffuse surface downwelling shortwave (W m-2)

real(sp) :: dayl    ! day length (h)

real(sp) :: sw_rad  ! total surface downwelling shortwave (W m-2)
real(sp) :: lw_rad  ! net longwave from surf_lw2 (Sandoval)(W m-2)
real(sp) :: lw_rad2 ! net longwave from surf_lw (Josey)(W m-2)

real(sp) :: lw_day   ! longwave radiation, daytime mean (W m-2)
real(sp) :: lw_night ! longwave radiation, nighttime mean (W m-2)

real(sp) :: lw_down ! downwelling longwave (W m-2)
real(sp) :: lw_up   ! upwelling longwave (W m-2)
real(sp) :: netrad  ! net radiation (kJ m-2 d-1)

real(sp) :: pet0    ! previous value for PET (mm d-1)
real(sp) :: dpet    ! day potential evapotranspiraton (mm)
real(sp) :: aet     ! actual evapotranspiration of the previous day (mm)

real(sp) :: tmin
real(sp) :: tmax

real(sp) :: phi     ! latitude (rad)
real(sp) :: slope   ! slope (rad)
real(sp) :: aspect  ! aspect (rad)
real(sp) :: hs      ! hour angle of sunset (rad)
real(sp) :: hn      ! net radiation cross-over hour angle (rad)
real(sp) :: sinhs   ! sin of the hour angle of sunset
real(sp) :: ru
real(sp) :: rv
real(sp) :: rw
real(sp) :: sw

! real(sp) :: Ilw
real(sp) :: HNpos   ! daytime accumulated net radiation (J m-2 d-1)
real(sp) :: HNneg

real(sp) :: sunf

real(sp) :: lwterm

real(sp) :: hour_sw
real(sp) :: hour_net

! counters

integer :: i

! ----------------------------------------------------------------------------------

toa_sw = dmet%rad0
dayl   = dmet%dayl
delta  = dmet%delta

albedo = dmet%Bsw

tmin   = dmet%tmin
tmax   = dmet%tmax
tday   = dmet%tday
tnight = dmet%tnight
cldf   = dmet%cldf
prec   = dmet%prec * dayl / 24. ! distribute 24-hr precipitation over the day and night (mm hr-1)
aet    = dmet%aet

lat    = pixel%lat
elv    = pixel%elv
phi    = pixel%phi    ! latitude (rad)
slope  = pixel%srad   ! slope inclination (rad)
aspect = pixel%gamma  ! slope orientation (rad), 0 = S, values increasing clockwise

tcm  = pixel%tcm
Pjj  = pixel%Pjj
Pann = pixel%Pann
Ratm = pixel%Ratm
P    = pixel%P

! ----
! hour angle when the solar radiation flux reaches the horizon or sunset hour (hs)

call hourangle(delta,phi,slope,aspect,ru,rv,hs,sinhs)

! ----
! atmospheric transparency (airmass) based on latitude, solar declination, daylength, and relative pressure

call airmass(lat,delta/pir,dayl,Ratm,air)

! Note: in the Yin (1997,1998) formulation for downwelling shortwave, radiation depends on PET, 
! presumably because evapotranspiration leads to atmospheric haze, which reduces surface radiation.
! This presents a problem in very dry environments with high PET but low soil moisture, so AET is
! relatively low. We solve this here by using the previous day's AET in place of the original PET.
! (January 2026)

! shortwave flux using the Yin method that empirically accounts for thickness of the atmosphere and haze

call surf_sw(Pjj,Ratm,toa_sw,cldf,air,albedo,prec,tcm,aet,direct,diffuse,sw_rad)

! Davis-Sandoval r_w term

rw = sw_rad * pi * (1. - albedo) / (ru * hs + rv * sin(hs))    ! Sandoval eqn. 8

! longwave flux - using Josey method

!   tdew = dewpoint(tmin,tmax,dpet,Pann)
!   call surf_lw(tday,tdew,cldf,lw_rad)

! longwave flux -- Sandoval method

sunf = sf(elv,toa_sw,sw_rad)         ! bright sunshine fraction as a function of elevation and sunlight attenuation Sandoval eqn 12

call surf_lw2(sunf,tday,lw_day)      ! daytime longwave radiation

call surf_lw2(sunf,tnight,lw_night)  ! nighttime longwave radiation

! hour angle of net radiation crossover (shortwave = longwave) eqn 13

lwterm = (lw_day - rw * ru) / (rw * rv)

if (lwterm <= -1.) then
  hn = pi
else if (lwterm >= 1.) then
  hn = 0
else
  hn = acos(lwterm)
end if

! positive net radiation (daytime) eqn 14

HNpos = pisec * ((rw * ru - lw_day) * hn + rw * rv * sin(hn))

! negative net radiation (nighttime) eqn 15

HNneg = pisec * (rw * rv * (sin(hs) - sin(hn)) + rw * ru * (hs - hn) - lw_night * (pi - hn))

! potential evapotranspiration

call pet(HNpos,dpet)



hour_sw  = 24. * hs / (2. * pi)
hour_net = 24. * hn / (2. * pi)

! write(0,'(2f7.1,3f7.3,4f7.1,2f7.3,2f12.1)')toa_sw,sw_rad,albedo,cldf,sunf,tday,lw_day,tnight,lw_night,hour_sw,hour_net,HNpos,HNneg



lw_rad = lw_day


! Isw = sw_rad * (1. - albedo)

! netrad = Isw - lw_rad

call pet(P,tday,lw_rad,ru,rv,rw,hn,dpet)

! write(0,'(f8.2,2f8.3,6f8.2)')toa_sw,cldf,sunf,direct,diffuse,sw_rad,lw_rad,dpet,aet

! call surf_lw(tday,tdew,cldf,lw_rad)

! call netrad_pet(sw_rad,lw_rad,albedo,P,tday,netrad,dpet)

! write(0,*)i,hs,cldf,sunf,tday,rad0,direct,diffuse,sw_rad,lw_rad,Hnpos,Hnneg,dpet

! write(0,*)i,Pann,tcm,cldf,tmin,tmax,tday,tdew,dpet,pet0,rhum(tday,tdew)


! ! also calculate josey method for comparison output
!   tdew = dewpoint(tmin,tmax,dpet,Pann)
!   call surf_lw(tday,tdew,cldf,lw_rad)

dmet%rdirect  = direct
dmet%rdiffuse = diffuse
dmet%dpet     = dpet
dmet%HNpos    = HNpos
dmet%HNneg    = HNneg
dmet%sunf     = sunf
dmet%hour_sw  = hour_sw
dmet%hour_net = hour_net
dmet%lwday    = lw_day
dmet%lwnight  = lw_night
! radiation outputs for diagnostics
dmet%swrad   = sw_rad   ! total surface downwelling shortwave (W m-2)
dmet%lw_rad  = lw_rad    ! net longwave from surf_lw2 / Sandoval (W m-2) - USED IN PHYSICS
!dmet%lw_rad2 = lw_rad2   ! net longwave from surf_lw / Josey (W m-2) - for comparison


! night timestep

end subroutine radpet

! ----------------------------------------------------------------------------------------------------------------

subroutine surf_sw(Pjj,Ratm,rad0,cldf,air,albedo,prec,tcm,pet,direct,diffuse,sw_rad)

! This code is based on:
! Yin, X. (1997). Optical air mass: Daily integration and its applications.
!   Meteorology and Atmospheric Physics, 63(3-4), 227-233. doi:10.1007/Bf01027387
! Yin, X. W. (1998). Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!   Systems development and application. Meteorology and Atmospheric Physics, 68(1-2), 99-113. doi:10.1007/Bf01025387

use parametersmod, only : sp,pi => pi_sp
use typesmod,      only : airmasspars

implicit none

! arguments

real(sp),          intent(in)  :: Pjj      ! precipitation equitability index (see function)
real(sp),          intent(in)  :: Ratm     ! relative atmospheric pressure (see function)
real(sp),          intent(in)  :: rad0     ! top-of-atmospere insolation (W m-2)
real(sp),          intent(in)  :: cldf     ! total cloud cover (fraction)
type(airmasspars), intent(in)  :: air      ! airmass parameters
real(sp),          intent(in)  :: albedo   ! surface shortwave albedo
real(sp),          intent(in)  :: prec     ! precipitation mm/day
real(sp),          intent(in)  :: tcm      ! temperature of the coldest month (used as tropics indicator)
real(sp),          intent(in)  :: pet      ! potential evapotranspiration mm/day
real(sp),          intent(out) :: direct   ! direct-beam downwelling shortwave (W m-2)
real(sp),          intent(out) :: diffuse  ! diffuse downwelling shortwave W m-2)
real(sp),          intent(out) :: sw_rad   ! total downwelling shortwave W m-2)

! parameters

real(sp), parameter :: kp  = 0.500  ! curve of the relationship relating absorption coeff. to trans. coeff.
real(sp), parameter :: kag = 3.300  ! ground albedo parameter
real(sp), parameter :: kan = 2.320  ! cloud albedo parameter
real(sp), parameter :: kn  = 0.686  ! cloud albedo parameter

! local variables

real(sp) :: mbar   ! daytime mean optical air mass (unitless, 1 at equatorial noon), calculated in airmassmod
real(sp) :: mo     ! air mass at cosine zenith angle maximum
real(sp) :: mc     ! air mass at cosine zenith angle medium
real(sp) :: ml     ! air mass at cosine zenith angle bottom quarter range point

real(sp) :: tau    ! transmission coefficient for direct radiation
real(sp) :: zeta   ! transmission coefficient for diffuse radiation
real(sp) :: fm     ! airmass function

real(sp) :: sunf   ! bright sunshine duration fraction, n/N (fraction)
real(sp) :: sunf2   ! bright sunshine duration fraction, n/N (fraction)

! ----------------------------------------------------------------------------

mbar = air%mbar
mo   = air%mo
mc   = air%mc
ml   = air%ml

sunf = 1. - cldf  ! as defined by Yin (1998), pg. 101

! ----
! direct radiation

tau = tau0(tcm,Ratm,prec,pet,Pjj)

fm = 0.01452 * (mbar + ml) * exp(1.403 * tau) - 0.1528 * mo + mc + 0.4870 * (mc - ml) + 0.2323   ! Yin (1997) eqn 3.3

direct = sunf * tau**kp * rad0 * tau**fm   ! eqn. 2.4

! ----
! diffuse radiation

zeta = zeta0(Ratm,prec,pet,albedo,cldf)

! this was my misinterpretation of eqn 2.5, which led to unphysical values for diffuse radiation 
! with surface shortwave > top-of-the-atmosphere forcing, especially under conditions of high snow and cloud cover

! diffuse = zeta * kag**albedo * kan**cldf * (1. - kn * cldf) * (tau**kp * rad0 - direct)   ! Eqn. 2.5 NOT USED

diffuse = zeta * (tau**kp * rad0 - direct)  ! first part of eqn 2.5

! ----

sw_rad = direct + diffuse

! if (sw_rad > rad0) then
!   write(0,*)cldf,albedo,rad0,sw_rad,direct,diffuse,zeta
! end if

end subroutine surf_sw

! ----------------------------------------------------------------------------------------------------------------

real(sp) function tau0(tcm,Ratm,prec,pet,Pjj)

! Estimate the transmission coefficient of the atmosphere for direct solar irradiance,
! which is also called the clear-sky atmospheric transmittance. Based on:
! Yin, X. (1998). Temporally-aggregated atmospheric optical properties as a function of common climatic information:
! Systems development and application. Meteorology and Atmospheric Physics, 68(1-2), 99-113. doi:10.1007/Bf01025387
 
use parametersmod, only : sp,pi => pi_sp

implicit none

! arguments

real(sp), intent(in) :: tcm  ! temperature of the coldest month (used as tropics indicator)
real(sp), intent(in) :: Ratm ! relative atmospheric pressure (see function)
real(sp), intent(in) :: prec ! precipitation (mm/day)
real(sp), intent(in) :: pet  ! potential evapotranspiration (mm/day)
real(sp), intent(in) :: Pjj  ! precipitation equitability index (see function)

! local variables

real(sp) :: x  ! tropics indicator (1 if tropical, 0 if extratropical)

! ----
! determine if the region of interest is "tropical" or not based on temperature of the coldest month

if (tcm < 10.) then
  x = 0.
else if (tcm > 20.) then
  x = 1.
else
  x = sin(pi / 2. * (tcm / 10. - 1.))
end if

! calculate tau

tau0 = exp(-0.115 * Ratm * ((2.15 - 0.713 * x + exp(-6.74 / (prec + 1.))) * exp(0.0971 * pet) - 0.650 * (1. - x) * Pjj))  ! eqn 4.1

end function tau0

! ----------------------------------------------------------------------------------------------------------------

real(sp) function zeta0(Ratm,prec,pet,albedo,cldf)

! Estimate the transmission coefficient of the atmosphere for diffuse solar irradiance, based on
! Yin, X. (1998). Temporally-aggregated atmospheric optical properties as a function of common climatic information: Systems development and application. 
! Meteorology and Atmospheric Physics, 68(1-2), 99-113. doi:10.1007/Bf01025387

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: Ratm   ! relative atmospheric pressure (see function)
real(sp), intent(in) :: prec   ! precipitation (mm/day)
real(sp), intent(in) :: pet    ! evapotranspiration (mm/day) (Yin calls for PET but we will use AET)
real(sp), intent(in) :: albedo ! surface shortwave albedo (fraction)
real(sp), intent(in) :: cldf   ! cloud cover (fraction)

! parameters

real(sp), parameter :: kag = 3.30   ! surface albedo multiplier (Fig. 6a)
real(sp), parameter :: kan = 2.32   ! cloud albedo reflection multiplier (Fig. 6b)
real(sp), parameter :: kn  = 0.686  ! cloud albedo filtering multiplier

! ----
! eqn 4.2

zeta0 = 0.503 * exp(-1.20 * Ratm * exp(-0.633 / (prec + 1.) - 0.226 * pet)) * kag**albedo * kan**cldf * (1. - kn * cldf)
        
end function zeta0

! ----------------------------------------------------------------------------------------------------------------

subroutine surf_lw(tair,tdew,cldf,lw_rad)

! estimate net longwave radiation flux (W m-2) (negative down) based on
! Josey, S. A. (2003). A new formula for determining the atmospheric longwave flux at the ocean surface at mid-high latitudes. 
! Journal of Geophysical Research, 108(C4). doi:10.1029/2002jc001418

use parametersmod, only : sp,tfreeze

implicit none

! arguments

real(sp), intent(in)  :: tair   ! 2m air temperature integrated over the reference period (degC)
real(sp), intent(in)  :: tdew   ! mean daily (24-hr) dewpoint temperature (degC)
real(sp), intent(in)  :: cldf   ! cloud cover fraction
real(sp), intent(out) :: lw_rad

! parameters

real(sp), parameter :: sb = 5.6704e-8  ! Stefan-Bolzmann constant (W m-2 K-4)
real(sp), parameter :: e  = 0.98       ! emissivity (unitless)
real(sp), parameter :: al = 0.045      ! longwave reflectivity (lw albedo), Josey et al., pg 5-9

real(sp), parameter :: a  =  10.77     ! parameters in Josey et al.
real(sp), parameter :: b  =   2.34
real(sp), parameter :: c  = -18.44


! local variables

real(sp) :: Ts  ! ground surface (skin) temperature (K)
real(sp) :: Ta  ! surface air temperature (K)
real(sp) :: Td  ! dewpoint temperature (K)
real(sp) :: D   ! dew point depression (Td - Ta) (K)

real(sp) :: Qdn ! downwelling 
real(sp) :: Qup ! upwelling

! ----

Ta = tair + tfreeze
Td = tdew + tfreeze
D  = Td - Ta

Ts = Ta  ! simplify skin temperature = air temperature

Qup = e * sb * Ts**4

Qdn = sb * (Ta + a*cldf**2 + b*cldf + c + 0.84 * (D + 4.01))**4  ! eqn 14,J2

lw_rad = Qup - (1. - al) * Qdn  ! eqn 1

! write(0,*)'longwave',Qup,Qdn

end subroutine surf_lw

! ----------------------------------------------------------------------------------------------------------------

subroutine surf_lw2(sunf,Tair,lw_rad)

! longwave radiation (sf) from eqn 12 of
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in)  :: sunf    ! bright sunshine fraction
real(sp), intent(in)  :: Tair    ! air temperature (degC)
real(sp), intent(out) :: lw_rad

! parameters

real(sp), parameter :: k1 = 91.86  ! (degC)
real(sp), parameter :: k2 =  1.96
real(sp), parameter :: k3 =  0.20
real(sp), parameter :: k4 =  0.088

! ----

lw_rad = (k4 + (1. - k3) * sunf) * (k1 + k2 * Tair)

end subroutine surf_lw2

! ----------------------------------------------------------------------------------------------------------------

subroutine pet(P,Tair,HNpos,lw_rad,ru,rv,rw,dpet,dpmax)

! Estimation of daily total equilibrium evapotranspiration (dpet), mm
! and hourly maximum evapotranspiration rate (mm h-1)
! based on equations in Sandoval et al. and Davis et al. and source code

! Although Sandoval provides an equation for evaporative demand (Dp), it is not used directly in the calculation of  
! actual evapotranspiration. AET is instead is calculated as a function of mean daytime shortwave and longwave fluxes,
! solar angle parameters, and evaporative supply rate (Sw). Sw is calculated in-turn as a function of the maximum 
! hourly evaporative demand (DpMAX), which is a function of solar angle variables and mean daytime longwave.

! 

use parametersmod, only : sp
use physicsmod,    only : Econ

implicit none

! arguments

real(sp), intent(in)  :: P       ! mean air pressure (Pa)
real(sp), intent(in)  :: Tair    ! air temperature (degC)
real(sp), intent(in)  :: HNpos   ! daytime accumulated net radiation (J m-2 d-1)

real(sp), intent(in)  :: lw_rad  ! mean longwave radiation (W m-2)
real(sp), intent(in)  :: ru      ! aggregate variables related to radiation
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

rx = 3.6e6 * Ec / 1000.  ! double check units

dpmax = rx * ((rw * (ru + rv)) - lw_rad)

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

real(sp) function rhum(tair,tdew)

! Estimate relative humidity from air temperature and dewpoint temperature
! Based on Lawrence, M. G. (2005). The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air: 
! A Simple Conversion and Applications. 
! Bulletin of the American Meteorological Society, 86(2), 225-234. doi:10.1175/bams-86-2-225

use parametersmod, only : sp,tfreeze

implicit none

! arguments

real(sp), intent(in)  :: tair  ! air temperature (degC)
real(sp), intent(in)  :: tdew  ! dewpoint temperature (degC)

! parameter

real(sp), parameter :: Rw = 461.5  ! gas constant for water vapor (J K-1 kg-1)

! local variables

real(sp) :: Ta
real(sp) :: Td
real(sp) :: L       ! enthalpy of vaporization of water (J kg-1)

! ------

Ta = tair + tfreeze
Td = tdew + tfreeze

L = 1.91846e6 * (Ta / (Ta))**2  ! (J kg-1) Eqn. from Henderson-Sellers (1984)

rhum = 100. * exp(-L / (Rw * Ta * Td) * (Ta - Td)) ! eqn 12

end function rhum

! ----------------------------------------------------------------------------------------------------------------

real(sp) function sf(elv,rad0,sw_rad)  ! bright sunshine fraction (fraction)

! estimation of bright sunshine fraction based on cloud cover (sf) from eqn 11 of
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: elv      ! elevation (m) 
real(sp), intent(in) :: rad0     ! top of the atmosphere shortwave radiation 
real(sp), intent(in) :: sw_rad   ! surface shortwave radiation

! parameters

real(sp), parameter :: k5 = 0.1898
real(sp), parameter :: k6 = 0.7410
real(sp), parameter :: c = 1. / k6  ! exponent Sandoval et al. eqn 12

! local variables

real(sp) :: tau   ! total atmospheric transmittance: (ratio of total surface SW to TOA SW, includes clouds)
real(sp) :: tau0  ! clear-sky atmospheric transmittance (function of elevation)

real(sp) :: a
real(sp) :: b
real(sp) :: ab

! ----
! The equation for clear-sky transmittance (tau0) comes from line 170 of SOLAR.cpp in the rsplash github code.
! I cannot find any source for the exact formulation. The paper cites Allen (1996), which does contain a 
! similar formula for a variable called KT (eqn 2): [0.75 + 2.e-5 * elv], which has its original source in Allen et al. (1994):
!   Allen, R. G., Smith, M., Pereira, L. S., & Perrier, A. (1994). An update for the calculation of reference evapotranspiration. 
!   ICID Bulletin, 43(2), 35-92. 

tau0 = 0.75 * 1. + 2.67e-5 * elv

if (rad0 > 0) then

  tau = sw_rad / rad0

else

  tau = tau0

end if

a = tau - tau0 * k5   ! numerator Sandoval et al. eqn 12

b = tau0 * (1. - k5)  ! denominator Sandoval et al. eqn 12

! write(0,*)tau0,tau,a,b,c

ab = max(a / b,0.)

sf = ab**c       ! Sandoval et al. eqn 12

sf = max(min(sf,1.),0.)

end function sf

! ----------------------------------------------------------------------------------------------------------------

subroutine hourangle(delta,phi,slope,aspect,ru,rv,hs,sinh)

! calculate the hour angle and subcomponents of downwelling shortwave radiation
! Sandoval et al eqns 3-6

use parametersmod, only : sp,pi=>pi_sp

implicit none

! arguments

real(sp), intent(in)  :: delta  ! solar declination angle (rad)
real(sp), intent(in)  :: phi    ! geodetic latitude (rad)
real(sp), intent(in)  :: slope  ! slope inclination (rad)
real(sp), intent(in)  :: aspect ! slope orientation (rad), 0 for slopes oriented due south with values increasing clockwise)
real(sp), intent(out) :: ru 
real(sp), intent(out) :: rv
real(sp), intent(out) :: hs

! local variables

real(sp) :: sind
real(sp) :: sinp
real(sp) :: sins
real(sp) :: sing

real(sp) :: cosd
real(sp) :: cosp
real(sp) :: coss
real(sp) :: cosg

real(sp) :: a
real(sp) :: b
real(sp) :: c
real(sp) :: sinh

real(sp) :: t1
real(sp) :: t2
real(sp) :: t3

! ----

sind = sin(delta)
sinp = sin(phi)
sins = sin(slope)
sing = sin(aspect)

cosd = cos(delta)
cosp = cos(phi)
coss = cos(slope)
cosg = cos(aspect)

! --

a = sind * cosp * sins * cosg - sind * sinp * coss   ! eqn 6a

b = cosd * cosp * coss + cosd * sinp * sins * cosg   ! eqn 6b

c = cosd * sing * sins                               ! eqn 6c

sinh = sinhs(a,b,c)                                  ! eqn 5

! --

t1 = sind * sinp * coss           ! eqn 3, numerator first term
t2 = sind * cosp * sins * cosg    ! eqn 3, numerator second term
t3 = cosd * sing * sins * sinh    ! eqn 3, numerator third term

ru = t1 - t2 + t3                 ! eqn 3, numerator

rv = cosd * cosp * coss + cosd * sinp * sins *cosg   ! eqn 3, denominator

if (ru / rv >= 1.) then
  hs = pi
else if (ru / rv <= -1.) then
  hs = 0.
else
  hs = acos(-ru / rv)       ! eqn 4
end if

end subroutine hourangle

! ----------------------------------------------------------------------------------------------------------------

real(sp) function sinhs(a,b,c)

! Sandoval et al. (2024) eqn 5

use parametersmod, only : sp

implicit none

real(sp), intent(in)  :: a
real(sp), intent(in)  :: b
real(sp), intent(in)  :: c

! ----

sinhs = (a * c + b * sqrt(b**2 + c**2 - a**2)) / (b**2 + c**2)

end function sinhs

! ----------------------------------------------------------------------------------------------------------------

end module radiationmod
