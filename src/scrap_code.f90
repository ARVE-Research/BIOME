


! local variables

real(sp) :: Tk     ! surface air temperature (K)
real(sp) :: ss     ! rate of increase of saturated vapor pressure with temperature (desdT) (Pa K-1)

subroutine netrad_pet(P,Tair,sw_rad,lw_rad,albedo,netrad,pet)

use parametersmod, only : sp,tfreeze
use physicsmod,    only : desdT,psygamma,lvap,Econ

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
! real(sp) :: psygamma  ! psychrometer constant (Pa K-1)
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

real(sp) function sf(elv,rad0,direct,diffuse)  ! bright sunshine fraction (fraction)

! estimation of bright sunshine fraction based on cloud cover (sf) from eqn 11 of
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
real(sp), parameter :: c = 1. / k6  ! exponent Sandoval et al. eqn 12

! local variables

real(sp) :: tau   ! total atmospheric transmittance: (ration of total surface SW to TOA SW, includes clouds)
real(sp) :: tau0  ! clear-sky atmospheric transmittance (function of elevation)

real(sp) :: a
real(sp) :: b
real(sp) :: ab

! ----
! the equation for clear-sky transmittance (tau0) comes from line 170 of SOLAR.cpp in the rsplash github code.
! I cannot find any source for the exact formulation. The paper cites Allen (1996), which does contain a 
! similar formula for a variable called KT (eqn 2): [0.75 + 2.e-5 * elv], which has its original source in Allen et al. (1994):
!   Allen, R. G., Smith, M., Pereira, L. S., & Perrier, A. (1994). An update for the calculation of reference evapotranspiration. 
!   ICID Bulletin, 43(2), 35-92. 

tau0 = 0.75 * 1. + 2.67e-5 * elv

if (rad0 > 0) then

  tau = (direct + diffuse) / rad0

else

  tau = tau0

end if

a = tau - tau0 * k5   ! numerator Sandoval et al. eqn 12

b = tau0 * (1. - k5)  ! denominator Sandoval et al. eqn 12

! write(0,*)tau0,tau,a,b,c

ab = max(a / b,0.)

sf = (ab)**c       ! Sandoval et al. eqn 12

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

subroutine hourangle(delta,phi,slope,aspect,ru,rv,hs,sinh)

! calculate the hour angle and subcomponents of downwelling shortwave radiation
! Sandoval et al eqns 3-4

use parametersmod, only : sp

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

a = sind * cosp * sins * cosg - sind * sinp * coss

b = cosd * cosp * coss + cosd * sinp * sins * cosg

c = cosd * sing * sins

sinh = sinhs(a,b,c)

! --

t1 = sind * sinp * coss
t2 = sind * cosp * sins * cosg
t3 = cosd * sing * sins * sinh

ru = t1 - t2 + t3

rv = cosd * cosp * coss + cosd * sinp * sins *cosg

hs = acos(-ru / rv)

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

real(sp) function rwf(sw_rad,albedo,ru,hs,rv,sinhs)

! calculate rw, Sandoval et al. (2024) eqn 8

use parametersmod, only : sp,pi => pi_sp

implicit none

real(sp), intent(in)  :: sw_rad
real(sp), intent(in)  :: albedo
real(sp), intent(in)  :: ru
real(sp), intent(in)  :: hs
real(sp), intent(in)  :: rv
real(sp), intent(in)  :: sinhs

! ----

rwf = (sw_rad * pi * (1. - albedo)) / (ru * hs + rv * sinhs)

end function rwf

! ----------------------------------------------------------------------------------------------------------------

subroutine calcnetrad(Ilw,rw,ru,rv,hs,Hnpos,Hnneg)

! calculate net radiation, Sandoval et al. (2024) eqns 13-15

use parametersmod, only : sp,pi => pi_sp

implicit none

! arguments

real(sp), intent(in)  :: Ilw
real(sp), intent(in)  :: rw
real(sp), intent(in)  :: ru
real(sp), intent(in)  :: rv
real(sp), intent(in)  :: hs
real(sp), intent(out) :: Hnpos
real(sp), intent(out) :: Hnneg

! parameter

real(sp), parameter :: pisec = 86400. / pi

! local variable

real(sp) :: netrad
real(sp) :: hn

! ----

netrad = (Ilw - rw * ru) / (rw * rv)

if (netrad <= -1.) then
  hn = pi
else if (netrad >= 1.) then
  hn = 0.
else
  hn = acos(netrad)
end if

Hnpos = pisec * ((rw * ru - Ilw) * hn + rw * rv * sin(hn))

Hnneg = pisec * (rw * rv * (sin(hs) - sin(hn)) + rw * ru * (hs - hn) - Ilw * (pi - hn))

end subroutine calcnetrad





















toa_sw = rad0 * dayl * 3.6   ! convert W m-2 to kJ m-2 d-1

! calculate the airmass based on latitude, solar declination, daylength, and relative pressure

call airmass(lat,delta,dayl,Ratm,air)

! calculate daytime and nighttime mean longwave radiation (W m-2)


! calculate daytime and nighttime mean shortwave radiation (W m-2)


! iterate to estimate dewpoint, shortwave and longwave radiation, net radiation, and potential evapotranspiration

i = 1

pet0 = 1.
dpet = pet0

! do

  tdew = dewpoint(tmin,tmax,dpet,Pann)

  call surf_sw(Pjj,Ratm,toa_sw,cldf,air,albedo,prec,tcm,dpet,direct,diffuse)

!  sw_rad = direct + diffuse

  call surf_lw(tday,tdew,cldf,dayl,Ilw,lw_rad)
  
  ! lw_rad = surf_lw2(elv,toa_sw,direct,diffuse,tday)

!   call netrad_pet(P,tday,sw_rad,lw_rad,albedo,netrad,dpet)

  ! write(0,*)i,Pann,tcm,cldf,tmin,tmax,tday,tdew,dpet,pet0,rhum(tday,tdew)

!   if (abs(dpet - pet0) < 0.01 .or. i > 100) exit

!   pet0 = dpet

!   i = i + 1

! end do

delta = solar%delta
phi = real(pixel%lat,sp) * pir
slope = 0.
aspect = 0.
sw = sw_rad / (3.6 * dayl)

call hourangle(delta,phi,slope,aspect,ru,rv,hs,sinh)

rw = rwf(sw,albedo,ru,hs,rv,sinh)

call calcnetrad(Ilw,rw,ru,rv,hs,Hnpos,Hnneg)

write(0,*)hs,ru,rv,rw,Ilw,Hnpos,Hnneg

dmet%dayl = dayl
dmet%tdew = tdew
! dmet%rhum = rhum
dmet%rdirect = direct
dmet%rdiffuse = diffuse
dmet%dpet = dpet








! unused variables taken from surf_sw routine
! real(sp) :: j2w
! real(sp) :: fdif
! real(sp) :: stmp








      ! calculate 24-hr mean dewpoint
      
      dmet(i)%tdew = dewpoint(dmet(i)%tmin,dmet(i)%tmax,dpet,pann)
      
      ! calculate day- and night-time mean downwelling longwave radiation (W m-2)
      
      dmet(i)%lwday   = lw_rad(dmet(i)%tday,tdew,cldf)
      dmet(i)%lwnight = 
      

      ! write(*,'(i5,5f8.2)')doy,dmet(i)%tmin,dmet(i)%tday,dmet(i)%tmax,dmet(i)%tnight,dmet1(i)%tmin
      
      ! estimate dewpoint temperature 












subroutine surf_lw(tair,tdew,cldf,time,Ql,lw_rad)

! calculate daytime net longwave radiation (kJ m-2 d-1)
! This code is based on the paper:
! A. Haxeltine and Prentice, I.C., BIOME3..., Glob. Biogeochem. Cycles, 10, 693-709
! With a new calculations of:
! downwelling longwave   (Josey et al., 2003. J. Geophys. Res., 108(C4), 3108, doi:10.1029/2002JC001418)
! lvap                   (Henderson-Sellers, 1984. Quart. J. R. Met. Soc., 110, 1186-1190)
! Other references:
! Linacre (1968) Agr. Meteorol., 5, 49-63
! Prentice et al. (1993) Ecological Modelling, 65, 51-70.
! Jed Kaplan, EPFL, 2008, 2011

use parametersmod, only : tfreeze

implicit none

! arguments

real(sp), intent(in)  :: tair   ! 2m air temperature integrated over the reference period (degC)
real(sp), intent(in)  :: tdew   ! mean daily (24-hr) dewpoint temperature (degC)
real(sp), intent(in)  :: cldf   ! cloud cover fraction 
real(sp), intent(in)  :: time   ! reference period (h) (typically length of the day or night)
real(sp), intent(out) :: Ql     ! daily mean net longwave radiation (W m-2)
real(sp), intent(out) :: lw_rad ! daily total net longwave radiation (kJ m-2 d-1)

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

! real(sp) :: f      ! Linacre parameter (function of sunshine fraction)

real(sp) :: Ql_up  ! upwelling longwave radiation (W m-2)
real(sp) :: Ql_dn  ! downwelling longwave radiation (W m-2)

! real(sp) :: sunf   ! bright sunshine duration fraction, n/N (fraction)

! -------------------------------------------------

! sunf = 1. - cldf

tairK = tair + tfreeze

! f = 0.2 + 0.8 * sunf  ! Linacre Eqn. 7

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

! write(0,*)'lw',tairK,cldf,D,Ql_dn

Ql = (1. - al) * Ql_dn ! - Ql_up   ! Josey et al., Eqn 1

! ----

lw_rad = 0.001 * 3600. * time * Ql  ! daytime net longwave (kJ m-2 d-1)

end subroutine surf_lw








! 10 format(2i4,f8.3,f6.1,2f7.1,f7.3,f7.1,3f7.1,2f7.1,f7.1,2f10.1,f7.1)





subroutine getmonlen(orbit,yrlen,ndm)

! initialize baseline month lengths for 0ka

use parametersmod, only : dp,nmos,nd_365,nd_366,veqday_365,present_mon_noleap,veqday_366,present_mon_leap
use typesmod,      only : monthinfotype,orbitpars

implicit none

type(orbitpars),        intent(in)  :: orbit
real(dp),               intent(in)  :: yrlen
real(dp), dimension(:), intent(out) :: ndm

real(dp) :: veqday

integer,  dimension(nmos) :: imonlen
real(dp), dimension(nmos) :: rmonbeg 
real(dp), dimension(nmos) :: rmonmid
real(dp), dimension(nmos) :: rmonend

! ----
! calculate real month lengths

yrlen = real(nd_365,dp)

veqday = veqday_365

imonlen = nint(present_mon_noleap)

call monlen(yrlen,veqday,imonlen,orbit,ndm%noleapyr,rmonbeg,rmonmid,rmonend)

! ----
! calculate real month lengths for 0ka for leap year

yrlen = real(nd_366,dp)

veqday = veqday_366

imonlen = nint(present_mon_leap)

call monlen(yrlen,veqday,imonlen,orbit,ndm%leapyear,rmonbeg,rmonmid,rmonend)

! ----

end subroutine initmonlen







! calculate insolation based on this true longitude

call cwj(wd,e,pibar,eps,phi,ww,dayl)

w = ww

! write(0,*)'wjcal',datecal,e,eps,pibarh,phi,w


!(datecal,e,eps,pibarh,phi,w,dayl)


real(dp), intent(in)        :: lat



real(dp), intent(in)  :: datecal ! mean longitude of the day relative to the vernal equinox (degrees)


real(dp), intent(in)  :: e       ! orbital eccentricity (unitless)
real(dp), intent(in)  :: pibarh  ! longitude of perhelion relative to the equinox (radians)
real(dp), intent(in)  :: eps     ! orbital obliquity (radians)
real(dp), intent(in)  :: phi     ! latitude (radians)

real(sp), intent(out) :: w       ! mean daily top-of-the-atmosphere insolation (W m-2)
real(dp), intent(out) :: dayl    ! day length (h)




real(dp) :: pibar  ! longitude of perihelion 
real(dp) :: hlm0
real(dp) :: hlm1
real(dp) :: wd
real(dp) :: ww




      ! call insol(slon,orbit,latr,met_out(i,2)%rad0,met_out(i,2)%dayl)

! real(sp) :: tcm   ! temperature of the coldest month
! real(dp) :: rad0  ! TOA radiation
! real(dp) :: dayl  ! daylength (h)

! type(dayinfotype), dimension(2) :: dayinfo  ! for the present and next day
! integer,  dimension(nmos) :: imonlen
! real(dp), dimension(nmos) :: rmonlen
! real(dp), dimension(nmos) :: rmonbeg
! real(dp), dimension(nmos) :: rmonmid
! real(dp), dimension(nmos) :: rmonend

! logical, allocatable, dimension(:,:) :: valid
! logical :: leapyr



real(sp), parameter :: w   = 15.     ! solar angular velocity (degrees hr-1)
real(sp), parameter :: rw  = pir * w ! solar angular velocity (radians hr-1)

real(sp), parameter :: m0  =  1.   ! air mass at  0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  ! air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  ! air mass at 90 degree solar zenith angle

real(sp), parameter :: cos80 = real(cos(80._dp * pir))  ! (degrees)

real(sp), parameter :: albedo = 0.17  ! surface shortwave albedo (fraction)

! module shared variables

real(sp), dimension(3) :: c00  ! air mass coefficients for solar zenith angle <=80 degrees
real(sp), dimension(3) :: c80  ! air mass coefficients for solar zenith angle  >80 degrees

type airmasspars
  real(sp) :: Ratm    ! relative atmospheric pressure 1=sea level
  real(sp) :: mbar    ! daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(sp) :: mo      ! air mass at cosine zenith angle maximum
  real(sp) :: mc      ! air mass at cosine zenith angle medium
  real(sp) :: ml      ! air mass at cosine zenith angle bottom quarter range point
end type airmasspars



      !call insol(slon,orbit,latr,met_out(i,1)%rad0,met_out(i,1)%dayl)
      ! write(*,*)'2023',m,d,j,met_out(i)%prec,met_out(i)%tmin,met_out(i)%tmax,met_out(i)%cldf,met_out(i)%wind




!       write(0,*)'doy',doy,m,d
!       write(0,*)met_out(i,:)%dayl
!       write(0,*)met_out(i,:)%prec
!       write(0,*)met_out(i,:)%tmin
!       write(0,*)met_out(i,:)%tmax
!       write(0,*)met_out(i,:)%tday
!       write(0,*)met_out(i,:)%tnight
!       write(0,*)
!       write(*,*)'2022',m,d,j,met_out(i)%prec,met_out(i)%tmin,met_out(i)%tmax,met_out(i)%cldf,met_out(i)%wind


! rhum = 100. * exp(((1. - (tairK / tdewK)) * (L / Rw)) / tairK) ! eqn. 12










! ----------------------------------------------------------------------------------------------------------------

subroutine surf_lw(tair,tmin,cldf,time,lw_rad,tdew)

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
real(sp), intent(in)  :: tmin    ! 2m air temperature daily minimum (degC)
real(sp), intent(in)  :: cldf    ! cloud cover fraction 
real(sp), intent(in)  :: time    ! reference period (h) (typically length of the day or night)

real(sp), intent(out) :: lw_rad  ! daytime net longwave radiation (kJ m-2 d-1)
real(sp), intent(out) :: tdew    ! dew point temperature (based on input tmin) (C)

! parameters

real(sp), parameter :: sb = 5.6704e-8  ! Stefan-Bolzmann constant (W m-2 K-4)
real(sp), parameter :: e  = 0.98       ! emissivity (unitless)
real(sp), parameter :: al = 0.045      ! longwave reflectivity (lw albedo), Josey et al., pg 5-9

real(sp), parameter :: a  =  10.77     ! parameters in Josey et al.
real(sp), parameter :: b  =   2.34
real(sp), parameter :: c  = -18.44

real(sp), parameter :: cs = 1.5 ! shape parameter for the curve relating fractional cloud cover to fractional sunshine duration

! local variables

real(sp) :: Tk     ! surface air temperature (K)
real(sp) :: Ts     ! ground surface (skin) temperature (K)
real(sp) :: TdewK  ! dewpoint temperature (K)
real(sp) :: D      ! dew point depression (K)
real(sp) :: es     ! saturation vapor pressure

real(sp) :: f      ! Linacre parameter (function of sunshine fraction)

real(sp) :: Ql     ! net longwave radiation (W m-2)
real(sp) :: Ql_up  ! upwelling longwave radiation (W m-2)
real(sp) :: Ql_dn  ! downwelling longwave radiation (W m-2)

real(sp) :: sunf   ! bright sunshine duration fraction, n/N (fraction)

! -------------------------------------------------

sunf = 1. - cldf

Tk = temp + Tfreeze

f = 0.2 + 0.8 * sunf  ! Linacre Eqn. 7

! -------------------------------------------------
! calculate longwave radiation

Ts = Tk ! approximation that mean daily surface temperature equals air temp.

! black body upwelling longwave (W m-2)  ! various sources e.g., Oleson et al.

Ql_up = e * sb * Ts**4

! --
! Josey et al. (2003) formulation for downwelling longwave

! To estimate dewpoint temperature we use the day's minimum temperature
! this makes the asumption that there is a close correlation between Tmin and dewpoint
! see, e.g., Glassy & Running, Ecological Applications, 1994

es = 0.01 * esat(Tk) ! saturation vapor pressure (mbar)  

TdewK = 34.07 + 4157. / log(2.1718e8 / es)  ! Josey et al., Eqn. 10

D = TdewK - Tk

Ql_dn = sb * (Tk + a*cldf**2 + b*cldf + c + 0.84 * (D + 4.01))**4  ! downwelling longwave (W m-2) Josey et al. Eqn. 14,J2

Ql = Ql_up - (1. - al) * Ql_dn   ! Josey et al., Eqn 1

! ----

lw_rad = 0.001 * 3600. * time * Ql  ! daytime net longwave (kJ m-2 d-1)

tdew = TdewK - Tfreeze

end subroutine surf_lw




! do m = 1,nmos
!   write(0,*)m,ndm0noleap(m),ndm1noleap(m),cal%noleap%ndmr(m),cal%noleap%ndmi(m)
! end do
! 
! write(0,*)
! 
! do m = 1,nmos
!   write(0,*)m,ndm0leapyr(m),ndm1leapyr(m),cal%leapyr%ndmr(m),cal%leapyr%ndmi(m)
! end do
! 
! write(0,*)sum(cal%noleap%ndmi),sum(cal%leapyr%ndmi)



! --------------------
! Note: input values for temperature that are very close to zero can cause underflow in the 
! polynomial calculations in this routine. We therefore round the input values to the nearest 0.01.
! Also, NB Fortran sets the value of 0**0 to 1 (November 2024)

dm%tmin_mn = roundto(dm%tmin_mn,2)
dm%tmax_mn = roundto(dm%tmax_mn,2)

do i=1,4

  if (i == 4 .or. dm%tmin_mn <= tmin_sd_breaks(i)) then

    if (pday) then

!       write(stdout,*)'A1',i,dm%tmin_mn,tmin_sd_w(i,:)
!       flush(stdout)

      dm%tmin_sd = sum(tmin_sd_w(i,:) * dm%tmin_mn**exponents)  !the vector 'exponents' is a clever way of calculating a polynomial expansion

    else

!       write(stdout,*)'A2',i,dm%tmin_mn,tmin_sd_d(i,:)
!       flush(stdout)

      dm%tmin_sd = sum(tmin_sd_d(i,:) * dm%tmin_mn**exponents)

    end if

    exit

  end if
end do

! write(*,*)'B',i,dm%tmin_sd




!       ! start day-night loop
!       
!       do j = 1,2
! 
!         ! calculate downwelling shortwave radiation, will generally be zero at night
!         
!         ! calculate longwave radiation
!         
!         ! calculate potential evapotranspiration
!         
!         surface radiation budget and potential evapotranspiration
  
        call radpet(pixel(i),solar,albedo,dmet(i,1))
        
        ! write(0,10)m,d,dmet(i,1)

!       end do  ! day-night loop


! ---------------------------------------

real(sp) function Tt(z,phi)

use parametersmod, only : sp,dp

implicit none

! arguments

real(sp), dimension(:), intent(in) :: Tmax  ! vector of smoothed pseudo-daily maximum air temperature (degC)
real(sp),               intent(in) :: z     ! elevation (m)
real(dp),               intent(in) :: phi   ! latitude (degrees)

! parameters                  line 576      | from paper
!                            splash.point.R | (NB the sign of k9 and k10 is inconsistent between paper and code)
real(sp), parameter :: k7  = -0.4710405934  ! -0.5827
real(sp), parameter :: k8  =  1.0473543991  !  1.319
real(sp), parameter :: k9  = -0.0004596581  ! -4.18e-4
real(sp), parameter :: k10 = -0.0110592101  ! -1.140e-2

! 1/(1+exp(-0.4710405934+1.0473543991*as.numeric(tc)-as.numeric(elev)*0.0004596581-abs(as.numeric(lat))*0.0110592101))

! psnow(d) = 1. / (1. + exp(k7 + k8 * Tmax(d) + k9 * z + k10 * abs(phi)))

! local variables

! integer :: d

real(sp), dimension(size(Tmax)) :: psnow  ! probability of snowfall occurrence

! ----
! analytical solution of Sandoval et al. (2024) eqn 19 for Tt at psnow of 0.5

Tt = (-k7 - k9 * z - k10 * abs(phi)) / k8

end function Tt





! from weathergen subroutine
! precipitation occurrence and amount have been moved to a separate subroutine because they are calculated
! for the whole month at once


! 1) Precipitation occurrence

! write(0,*)'flag1',pre,wetf,wetd

! if there is precipitation this month, calculate the precipitation state for today

if (wetf > 0. .and. pre > 0.) then

  ! calculate transitional probabilities for dry to wet and wet to wet days
  ! Relationships from Geng & Auburn, 1986, Weather simulation models based on summaries of long-term data

  if (pday(1)) then !yesterday was raining, use p11

    pwet = p11_1 + p11_2 * wetf

  else if (pday(2)) then ! yesterday was not raining but the day before yesterday was raining, use p101

    pwet = p101_1 + p101_2 * wetf

  else  ! both yesterday and the day before were dry, use p001

    pwet = p001_1 + p001_2 * wetf

  end if

  ! -----
  ! determine the precipitation state of the current day using the Markov chain approach

  u = ranur(rndst)

  if (u <= pwet) then  ! today is a rain day

    pday = eoshift(pday,-1,.true.)

  else  !today is dry

    pday = eoshift(pday,-1,.false.)

  end if

  ! ---------------------------
  ! 2) precipitation amount

  if (pday(1)) then  ! today is a wet day, calculate the rain amount

    ! calculate the mean precipitation on days with 
    
    pbar = pre / wetd

    ! calculate parameters for the distribution function of precipitation amount

    g_scale = g_scale_coeff * pbar
    g_shape = pbar / g_scale

    call gamma_cdf(p_trans,0.,g_scale,g_shape,cdf_thresh)

    call gamma_pdf(p_trans,0.,g_scale,g_shape,pdf_thresh)

    gp_scale = real((1._dp - cdf_thresh) / pdf_thresh)

    i = 1
    do

      !today's precipitation
      
      prec = ran_gamma_gp(rndst,.true.,g_shape,g_scale,p_trans,gp_shape,gp_scale)

      ! enforce positive precipitation that is not more than 5% greater than the monthly total

      if (prec > 0. .and. prec <= 1.05 * pre) exit

      if (i == 1000) then
        write (0,*)'Could not find good precipitation with ', pre, ' mm and ', wetd, ' wet days'
        stop
      else
        i = i + 1
      end if

    end do

  else

    prec = 0.

  end if

else  ! there was no precipitation in this month, so none on this day either

  pday = .false.
  prec = 0.

end if


! deleted from main program 20260108

!       write(0,*)m,d,dmet0(i)%tday,dmet0(i)%tnight,dmet0(i)%prec,dmet0(i)%snow, & 
!                 dmet0(i)%melt,dmet0(i)%swe,dmet0(i)%fsnow,dmet0(i)%asnow,dmet0(i)%Bsw
      
      ! write(0,*)soilw(i)


      ! write(*,'(2i5,3f7.1)')m,d,dmet0(i)%tday,dmet0(i)%tnight,dmet0(i)%prec

      ! write(0,*)m,d,dmet0(i)%tday,dmet0(i)%tnight,dmet0(i)%prec,dmet0(i)%snow,dmet0(i)%melt,dmet0(i)%swe,dmet0(i)%fsnow,dmet0(i)%asnow,dmet0(i)%Bsw
      
      ! debug swe
      
      !call snow(pixel(i),dmet0(i))

      ! if (i == 1 .and. m == 1 .and. d == 1) then
 		! write(0,*) 'DEBUG tday:', dmet0(i)%tday
 		! write(0,*) 'DEBUG prec:', dmet0(i)%prec
 		! write(0,*) 'DEBUG snow:', dmet0(i)%snow
	    ! write(0,*) 'DEBUG melt:', dmet0(i)%melt
  		! write(0,*) 'DEBUG swe :', dmet0(i)%swe
	  ! end if

