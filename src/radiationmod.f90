module radiationmod

! consolidated module for calculation airmass, surface downwelling shortwave, net longwave, net radiation and PET.

use parametersmod, only : sp

implicit none

! module subroutines and functions

public  :: initairmass
public  :: Pjj
public  :: elev_corr
public  :: radpet
public  :: dewpoint

private :: airmass
private :: surf_sw
private :: surf_lw
private :: netrad_pet

private :: m
private :: F

! module calculated parameters

real(sp), dimension(3) :: c00  ! air mass coefficients for solar zenith angle <=80 degrees
real(sp), dimension(3) :: c80  ! air mass coefficients for solar zenith angle  >80 degrees

contains

! ----------------------------------------------------------------------------------------------------------------

subroutine initairmass()

! calculate parameters used in the airmass calculations

use parametersmod, only : sp,pir

implicit none

! parameters

real(sp), parameter :: m0  =  1.   ! air mass at  0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  ! air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  ! air mass at 90 degree solar zenith angle

real(sp), parameter :: cos80 = real(cos(80. * pir))  ! (degrees)

! ---

c00(1) = 0.008307
c00(2) = (m0 - m80) * (c00(1) + 1.) * (c00(1) + cos80) / (cos80 - 1.)
c00(3) = m0 - c00(2) / (c00(1) + 1.)

c80(1) = 0.037160
c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
c80(3) = m90 - c80(2) / c80(1)

end subroutine initairmass

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Pjj(mtemp,mprec)

! calculation of the "precipitation equitability index"

use parametersmod, only : sp
use utilitymod,    only : imaxloc,iminloc

implicit none

! arguments

real(sp), dimension(:), intent(in) :: mtemp  ! annual time series of monthly temperature
real(sp), dimension(:), intent(in) :: mprec  ! annual time series of monthly precipitation

! local variables

integer :: wm  ! index position of the warmest month
integer :: cm  ! index position of the warmest month

real(sp) :: p_wm  ! precipitation in the warmest month 
real(sp) :: p_cm  ! precipitation in the coldest month

! ------

wm   = imaxloc(mtemp)
cm   = iminloc(mtemp)
p_wm = mprec(wm)
p_cm = mprec(cm)   

if (p_wm + p_cm > 0.) then

  Pjj = 2. * (p_wm - p_cm) / (p_wm + p_cm)
  Pjj = max(Pjj,0.)

else

  Pjj = 0.

end if

end function Pjj

! ----------------------------------------------------------------------------------------------------------------

real(sp) function elev_corr(elevation)

implicit none

real(sp), intent(in)  :: elevation 

real(sp), parameter :: z0 = 1. / 8000.

! ----

elev_corr = exp(-elevation * z0)

end function elev_corr

! ----------------------------------------------------------------------------------------------------------------

subroutine radpet(pixel,solar,albedo,met)

! calculate surface radiation budget and potential evapotranspiration

use parametersmod, only : sp,dp,pir
use typesmod,      only : orbitpars,airmasspars,pixeltype,solarpars,metvars_daily

implicit none

! arguments

type(pixeltype),   intent(in)    :: pixel
type(solarpars),   intent(in)    :: solar
real(sp),          intent(in)    :: albedo  ! surface shortwave albedo
type(metvars_daily), intent(inout) :: met

! local variables

real(dp) :: lat   ! latitude (degrees)
real(sp) :: elv   ! elevation (m)
real(sp) :: tday  ! daytime mean temperature (C)
real(sp) :: prec  ! daily total precipitation
real(sp) :: cldf  ! cloud cover fraction
real(sp) :: tdew  ! estimated dewpoint temperature (degC)

real(sp) :: tcm   ! temperature of the coldest month
real(sp) :: Pann  ! total annual precipitation (mm) 
real(sp) :: Pjj   ! precipitation equitability index 
real(sp) :: Ratm  ! relative atmospheric pressure (based on elevation)
real(sp) :: P     ! mean atmospheric pressure (based on elevation)

type(airmasspars) :: air

real(sp) :: toa_sw  ! top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
real(sp) :: delta   ! solar declination (degrees)
real(sp) :: direct  ! direct beam surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: diffuse ! diffuse surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: lw_rad  ! net longwave (kJ m-2 d-1)

real(sp) :: rad0
real(sp) :: dayl    ! day length (h)
real(sp) :: sw_rad  ! total surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: netrad  ! net radiation (kJ m-2 d-1)

real(sp) :: pet0    ! previous value for PET (mm d-1)
real(sp) :: dpet    ! day potential evapotranspiraton (mm)

real(sp) :: tmin
real(sp) :: tmax

! counters

integer :: i

! ----------------------------------------------------------------------------------

rad0  = solar%rad0
dayl  = solar%dayl
delta = solar%delta / pir    ! convert radians to degrees

tmin = met%tmin
tmax = met%tmax
tday = met%tday
cldf = met%cldf
prec = met%prec * dayl / 24. ! distribute 24-hr precipitation over the day and night

lat  = pixel%lat
elv  = pixel%elv
tcm  = pixel%tcm
Pjj  = pixel%Pjj
Pann = pixel%Pann
Ratm = pixel%Ratm
P    = pixel%P

toa_sw = rad0 * dayl * 3.6   ! convert W m-2 to kJ m-2 d-1

! calculate the airmass based on latitude, solar declination, daylength, and relative pressure

call airmass(lat,delta,dayl,Ratm,air)

! iterate to estimate dewpoint, surface radiation budget, potential evapotranspiration

i = 1

pet0 = 1.
dpet = pet0

do

  tdew = dewpoint(tmin,tmax,dpet,Pann)

  call surf_sw(Pjj,Ratm,toa_sw,cldf,air,albedo,prec,tcm,dpet,direct,diffuse)

  sw_rad = direct + diffuse

  lw_rad = surf_lw(tday,tdew,cldf,dayl)
  
  ! lw_rad = surf_lw2(elv,toa_sw,direct,diffuse,tday)

  call netrad_pet(P,tday,sw_rad,lw_rad,albedo,netrad,dpet)

  ! write(0,*)i,Pann,tcm,cldf,tmin,tmax,tday,tdew,dpet,pet0,rhum(tday,tdew)

  if (abs(dpet - pet0) < 0.01 .or. i > 100) exit

  pet0 = dpet

  i = i + 1

end do

met%dayl = dayl
met%tdew = tdew
met%rdirect = direct
met%rdiffuse = diffuse
met%dpet = dpet

end subroutine radpet

! ----------------------------------------------------------------------------------------------------------------

subroutine airmass(lat,delta,dayl,Ratm,air)

! This code is based on the paper:
! X. Yin (1997) Optical air mass: Daily integration and its applications, Meteorol. Atmos. Phys. 63, 227-233
! Jed Kaplan, EPFL, 2008

use parametersmod, only : sp,dp,pir
use typesmod,      only : airmasspars

implicit none

! arguments

real(dp),          intent(in)  :: lat    ! latitude (degrees)
real(sp),          intent(in)  :: delta  ! solar declination (degrees)
real(sp),          intent(in)  :: dayl   ! day length (hours)
real(sp),          intent(in)  :: Ratm   ! relative atmospheric pressure
type(airmasspars), intent(out) :: air    ! airmass parameters

! parameters

real(sp), parameter :: cos80 = real(cos(80._dp * pir))  ! (degrees)
real(sp), parameter :: mindayl = 2. * tiny(0._sp)

real(sp), parameter :: m0  =  1.   ! air mass at  0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  ! air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  ! air mass at 90 degree solar zenith angle

real(sp), parameter :: w   = 15.     ! solar angular velocity (degrees hr-1)

! local variables

real(sp) :: mbar    ! daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      ! air mass at cosine zenith angle maximum
real(sp) :: mc      ! air mass at cosine zenith angle medium
real(sp) :: ml      ! air mass at cosine zenith angle bottom quarter range point

real(sp) :: rlat    ! latitude (radians)
real(sp) :: rdelta  ! solar declination (radians)

real(sp) :: t1      ! number of hours between sunrise/sunset and solar noon (hr)
real(sp) :: t80     ! solar hour corresponding to the 80 degree zenith angle

! real(sp) :: t       ! solar hour (hr)

real(sp) :: Z       ! solar zenith angle (degrees)
real(sp) :: Zn      ! lesser of solar zenith angle at sunset or at midnight (degrees)
real(sp) :: Z0      ! zenith angle at solar noon (degrees)
real(sp) :: cosZ    ! cosine solar zenith angle (fraction), used in calculation of instantaneous air mass

! real(sp) :: l

! integer :: steps    ! integer number of time steps
! integer :: i        ! counter

real(sp) :: sinlat
real(sp) :: coslat
real(sp) :: sindel
real(sp) :: cosdel

real(sp)               :: a   ! values in equation 2.6b
real(sp)               :: b
real(sp), dimension(3) :: c

real(sp) :: tmp1
real(sp) :: tmp2
real(sp) :: tmp3

real(sp) :: tinv

real(sp) :: rZ0
real(sp) :: rZn

! -------------------------------------
! calculate daily mean air mass (mbar)

if (dayl == 0.) then

  mbar = m90
  mc   = m90
  ml   = m90
  mo   = m90

else
  
  ! basic setup

  rlat   = pir * lat
  rdelta = pir * delta
  
  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  ! ------
  
  ! calculate the half-daylength
  
  if (dayl > mindayl) then
    t1 = 0.5 * dayl
  else
    t1 = dayl
  end if
  
  ! Eqn. 2.9
  if (abs(lat + delta) >= 90.) then
    Zn = acos(sinlat * sindel - coslat * cosdel) / pir
  else
    Zn = 90.
  end if

  ! Eqn. 2.10
  if (abs(lat - delta) >= 90.) then
    Z0 = 90.
  else
    Z0 = lat - delta
  end if
  
  rZ0 = Z0 * pir  ! convert to radians
  rZn = Zn * pir
  
  ! --------------------------

  b = coslat * cosdel
  
  if (t1 == 0.) then

    mbar = m90
    
  else if (abs(Zn) <= 80.) then
  
    tinv = 1. / t1

    c = c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80.) then
  
    tinv = 1. / t1

    c = c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else
    
    t80 = 1. / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / pir  ! Eqn. 2.8

    c = c00
    a = c(1) + sinlat * sindel
        
    tmp1 = F(t80,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)

    tinv = 1. / t1

    
    mbar = tinv * (tmp1 + tmp2 - tmp3)
    
  end if

  ! --------------------------
  ! calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)
  
  Z = Z0

  cosZ = cos(Z * pir)

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if
  
  mo = m(cosZ,c)

  ! --

  Z = (Z0 + Zn) / 2.
  
  cosz = (cos(rZ0) + cos(rZn)) / 2.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  mc = m(cosZ,c)

  ! --

  Z = (Z0 + 3. * Zn) / 4.
  
  cosz = (cos(rZ0) + 3. * cos(rZn)) / 4.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  ml = m(cosZ,c)

end if

! -------------------------------------
! correct calculated air mass for elevation

air%mbar = Ratm * mbar
air%mo   = Ratm * mo
air%mc   = Ratm * mc
air%ml   = Ratm * ml

end subroutine airmass

! ----------------------------------------------------------------------------------------------------------------

subroutine surf_sw(Pjj,Ratm,rad0,cldf,air,albedo,prec,tcm,pet,direct,diffuse)

! This code is based on the papers:
! Yin, X. (1997). Optical air mass: Daily integration and its applications.
!   Meteorology and Atmospheric Physics, 63(3-4), 227-233. doi:10.1007/Bf01027387
! Yin, X. W. (1998). Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!   Systems development and application. Meteorology and Atmospheric Physics, 68(1-2), 99-113. doi:10.1007/Bf01025387

use parametersmod, only : sp,pi
use typesmod,      only : airmasspars

implicit none

! arguments

real(sp),          intent(in)  :: Pjj      ! precipitation equitability index
real(sp),          intent(in)  :: Ratm     ! relative atmospheric pressure (see function)
real(sp),          intent(in)  :: rad0     ! top-of-atmospere insolation (kJ m-2 d-1)
real(sp),          intent(in)  :: cldf     ! total cloud cover (fraction)
type(airmasspars), intent(in)  :: air      ! airmass parameters
real(sp),          intent(in)  :: albedo   ! surface shortwave albedo
real(sp),          intent(in)  :: prec     ! precipitation mm/day
real(sp),          intent(in)  :: tcm      ! temperature of the coldest month (used as tropics indicator)
real(sp),          intent(in)  :: pet      ! potential evapotranspiration mm/day
real(sp),          intent(out) :: direct   ! direct-beam downwelling shortwave (kJ m-2 d-1)
real(sp),          intent(out) :: diffuse  ! diffuse downwelling shortwave (kJ m-2 d-1)

! local variables

real(sp) :: mbar   ! daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo     ! air mass at cosine zenith angle maximum
real(sp) :: mc     ! air mass at cosine zenith angle medium
real(sp) :: ml     ! air mass at cosine zenith angle bottom quarter range point

real(sp) :: tau    ! direct insolation atmospheric turbidity factor
real(sp) :: zeta0  ! diffuse insolation atmospheric turbidity factor
real(sp) :: x      ! tropics indicator (tropical = 1, else 0)
real(sp) :: fm     ! atmospheric transmittance function, this is also called tau0 in some papers

! real(sp) :: j2w
! real(sp) :: fdif
! real(sp) :: stmp
real(sp) :: sunf   ! bright sunshine duration fraction, n/N (fraction)

! -----------------------------------
! parameters

real(sp), parameter :: kp  = 0.500  ! links absorption coeff. to trans. coeff.
real(sp), parameter :: kag = 3.300
real(sp), parameter :: kan = 2.320
real(sp), parameter :: kn  = 0.686  ! cloud parameter

! ----------------------------------------------------------------------------

mbar = air%mbar
mo   = air%mo
mc   = air%mc
ml   = air%ml

! ------

sunf = 1. - cldf

if (tcm < 10.) then
  x = 0.
else if (tcm > 20.) then
  x = 1.
else
  x = sin(pi / 2. * (tcm / 10. - 1.))
end if

! Yin Eqn. 4.1
tau = exp(-0.115 * Ratm * ((2.15 - 0.713 * x + exp(-6.74 / (prec + 1.))) * exp(0.0971 * pet) - 0.650 * (1. - x) * Pjj))

fm = 0.01452 * (mbar + ml) * exp(1.403 * tau) - 0.1528 * mo + mc + 0.4870 * (mc - ml) + 0.2323   ! Yin (1997) eqn 3.3

direct = sunf * tau**kp * rad0 * tau**fm   ! Eqn. 2.4

! Yin Eqn. 4.2
zeta0 = 0.503 * exp(-1.20 * Ratm * exp(-0.633 / (prec + 1.) - 0.226 * pet)) * kag**albedo * kan**(1. - sunf) * & 
        (1. - kn * (1. - sunf))

diffuse = zeta0 * kag**albedo * kan**(1. - sunf) * (1 - kn * (1. - sunf)) * (tau**kp * rad0 - direct)   ! Eqn. 2.5

end subroutine surf_sw

! ----------------------------------------------------------------------------------------------------------------

real(sp) function surf_lw(tair,tdew,cldf,time)

! calculate daytime net longwave radiation (kJ m-2 d-1)
! This code is based on the paper:
! A. Haxeltine and Prentice, I.C., BIOME3..., Glob. Biogeochem. Cycles, 10, 693-709
! With a new calculations of:
! downwelling longwave   (Josey et al., 2003. J. Geophys. Res., 108(C4), 3108, doi:10.1029/2002JC001418)
! dEsat/dT               (Oleson et al., 2004, CLM 3.0 technical note)
! lvap                   (Henderson-Sellers, 1984. Quart. J. R. Met. Soc., 110, 1186-1190)
! Other references:
! Linacre (1968) Agr. Meteorol., 5, 49-63
! Prentice et al. (1993) Ecological Modelling, 65, 51-70.
! Jed Kaplan, EPFL, 2008, 2011

use parametersmod, only : tfreeze

implicit none

! arguments

real(sp), intent(in)  :: tair    ! 2m air temperature integrated over the reference period (degC)
real(sp), intent(in)  :: tdew    ! dewpoint temperature (degC)
real(sp), intent(in)  :: cldf    ! cloud cover fraction 
real(sp), intent(in)  :: time    ! reference period (h) (typically length of the day or night)

! parameters

real(sp), parameter :: sb = 5.6704e-8  ! Stefan-Bolzmann constant (W m-2 K-4)
real(sp), parameter :: e  = 0.98       ! emissivity (unitless)
real(sp), parameter :: al = 0.045      ! longwave reflectivity (lw albedo), Josey et al., pg 5-9

real(sp), parameter :: a  =  10.77     ! parameters in Josey et al.
real(sp), parameter :: b  =   2.34
real(sp), parameter :: c  = -18.44

real(sp), parameter :: cs = 1.5 ! shape parameter for the curve relating fractional cloud cover to fractional sunshine duration

! local variables

real(sp) :: tairK  ! surface air temperature (K)
real(sp) :: tskin  ! ground surface (skin) temperature (K)
real(sp) :: tdewK  ! dewpoint temperature (K)
real(sp) :: D      ! dew point depression (K)
real(sp) :: es     ! saturation vapor pressure

real(sp) :: f      ! Linacre parameter (function of sunshine fraction)

real(sp) :: Ql     ! net longwave radiation (W m-2)
real(sp) :: Ql_up  ! upwelling longwave radiation (W m-2)
real(sp) :: Ql_dn  ! downwelling longwave radiation (W m-2)

real(sp) :: sunf   ! bright sunshine duration fraction, n/N (fraction)

! -------------------------------------------------

sunf = 1. - cldf

tairK = tair + tfreeze

f = 0.2 + 0.8 * sunf  ! Linacre Eqn. 7

! -------------------------------------------------
! calculate longwave radiation

tskin = tairK ! approximation that over vegetated surfaces mean daily skin temperature equals air temp

! black body upwelling longwave (W m-2)  ! various sources e.g., Oleson et al.

Ql_up = e * sb * tskin**4

! --
! Josey et al. (2003) formulation for downwelling longwave

tdewK = tdew + tfreeze

D = tdewK - tairK

Ql_dn = sb * (tairK + a*cldf**2 + b*cldf + c + 0.84 * (D + 4.01))**4  ! downwelling longwave (W m-2) Josey et al. Eqn. 14,J2

Ql = (1. - al) * Ql_dn - Ql_up   ! Josey et al., Eqn 1

! ----

surf_lw = 0.001 * 3600. * time * Ql  ! daytime net longwave (kJ m-2 d-1)

end function surf_lw

! ----------------------------------------------------------------------------------------------------------------

subroutine netrad_pet(P,Tair,sw_rad,lw_rad,albedo,netrad,pet)

use parametersmod, only : sp,tfreeze
use physicsmod,    only : desdT,gamma,lvap,Econ

implicit none

! arguments

real(sp), intent(in)  :: P       ! mean air pressure (Pa)
real(sp), intent(in)  :: Tair    ! air temperature (degC)
real(sp), intent(in)  :: sw_rad  ! downwelling shortwave radiation (kJ m-2 d-1)
real(sp), intent(in)  :: lw_rad  ! net longwave radiation (kJ m-2 d-1)
real(sp), intent(in)  :: albedo  ! surface shortwave albedo (fraction)
real(sp), intent(out) :: netrad  ! net radiation (kJ m-2 d-1)
real(sp), intent(out) :: pet     ! potential evapotranspiration (mm d-1)

! parameter

real(sp), parameter :: omega = 0.26  ! entrainment factor, dimensionless (Sandoval et al., 2024; Priestley and Taylor, 1972)

! local variables

real(sp) :: Tk     ! surface air temperature (K)
real(sp) :: ss     ! rate of increase of saturated vapor pressure with temperature (desdT) (Pa K-1)

! functions used

! real(sp) :: lvap   ! Latent heat of vaporization of water (temperature dependent) (kJ kg-1)
! real(sp) :: gamma  ! psychrometer constant (Pa K-1)
! real(sp) :: Econ   ! energy-to-water conversion factor (m3 kJ-1)

! ----
! calculate gamma, lvap, and ss

! Tk = temp + Tfreeze
! 
! ss = desdT(Tk)

! calculate net radiation and PET

netrad = (1. - albedo) * sw_rad - lw_rad             ! (kJ m-2 d-1)

pet = (1. + omega) * 1000. * Econ(P,Tair) * max(netrad,0.)  ! 1.e3 converts m to mm

! pet = max((ss / (ss + gamma(temp))) * netrad / lvap(temp), 0.)   ! (mm d-1)

end subroutine netrad_pet

! ----------------------------------------------------------------------------------------------------------------

real(sp) function m(cosZ,c)

! Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(sp),               intent(in) :: cosZ
real(sp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

! ----------------------------------------------------------------------------------------------------------------

real(sp) function F(t1,a,b,c)

! integral air mass function F, equation 2.6b in Yin, 1997
! section inside curly braces only - multiply result by 1/t1 to get mbar

use parametersmod, only : sp,pi,pir

implicit none

! arguments

real(sp),               intent(in) :: t1
real(sp),               intent(in) :: a
real(sp),               intent(in) :: b
real(sp), dimension(:), intent(in) :: c

! parameters

real(sp), parameter :: w   = 15.      ! solar angular velocity (degrees hr-1)
real(sp), parameter :: rw  = pir * w  ! solar angular velocity (radians hr-1)

! local variables

real(sp) :: wt1
real(sp) :: wpi

real(sp) :: e1
real(sp) :: e2

! ----

wpi  = 180. / (pi * w)
wt1  = rw * t1

if (a > b) then
  
  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then
  
  e1 = sqrt((b + a) * (1. + cos(wt1))) + sqrt((b - a) * (1. - cos(wt1)))
  e2 = sqrt((b + a) * (1. + cos(wt1))) - sqrt((b - a) * (1. - cos(wt1)))
  
  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2.) + c(3) * t1

end if

end function F

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

real(sp) :: tairK
real(sp) :: tdewK
real(sp) :: L       ! enthalpy of vaporization of water (J kg-1)

! ------

tairK = tair + tfreeze
tdewK = tdew + tfreeze

L = 1.91846e6 * (tairK / (tairK))**2  ! (J kg-1) Eqn. from Henderson-Sellers (1984)

! tdewK = min(tdewK,tairK)  ! limit dewpoint temperature to air temperature

rhum = 100. * exp(-L / (Rw * tairK * tdewK) * (tairK - tdewK)) ! eqn 12

end function rhum

! ----------------------------------------------------------------------------------------------------------------

real(sp) function sf(elv,rad0,direct,diffuse)  ! bright sunshine fraction (fraction)

! sunshine fraction (sf) from eqn 11 of
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: elv      ! elevation (m) 
real(sp), intent(in) :: rad0     ! 
real(sp), intent(in) :: direct   ! 
real(sp), intent(in) :: diffuse  ! 

! parameters

real(sp), parameter :: k5 = 0.1898
real(sp), parameter :: k6 = 0.7410
real(sp), parameter :: c = 1. / k6

! local variables

real(sp) :: tau   ! total atmospheric transmittance: (ration of total surface SW to TOA SW)
real(sp) :: tau0  ! clear-sky atmospheric transmittance: (ratio of surface direct SW to TOA SW)

real(sp) :: a
real(sp) :: b

! ----

tau0 = 0.75 * 1. + 2.67e-5 * elv

if (rad0 > 0) then

  tau = (direct + diffuse) / rad0

else

  tau = tau0

end if

a = tau - tau0 * k5

b = tau0 * (1. - k5)

sf = (a / b)**c

end function sf

! ----------------------------------------------------------------------------------------------------------------

real(sp) function surf_lw2(elv,rad0,direct,diffuse,Tair)

! longwave radiation (sf) from eqn 12 of
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
! Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
! Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in)  :: elv      ! elevation (m)
real(sp), intent(in)  :: rad0     ! 
real(sp), intent(in)  :: direct   ! 
real(sp), intent(in)  :: diffuse  ! 
real(sp), intent(in)  :: Tair     ! air temperature (degC)

! parameters

real(sp), parameter :: k1 = 91.86  ! (degC)
real(sp), parameter :: k2 =  1.95
real(sp), parameter :: k3 =  0.20
real(sp), parameter :: k4 =  0.088

! local variables

real(sp) :: sunf

! ----

sunf = sf(elv,rad0,direct,diffuse)

surf_lw2 = (k4 + (1. - k3) * sunf) * (k1 + k2 * Tair)

end function surf_lw2

! ----------------------------------------------------------------------------------------------------------------

end module radiationmod
