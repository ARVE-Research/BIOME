module snowmod

implicit none

public :: Tt

contains

! ---------------------------------------

subroutine snow(pixel,dmet)

! calculate daily snow dynamics

use parametersmod, only : sp,B0
use typesmod,      only : pixeltype,metvars_daily

implicit none

! arguments

type(pixeltype),     intent(in)    :: pixel
type(metvars_daily), intent(inout) :: dmet    ! daily meteorological variables

! parameters

real(sp), parameter :: k11   =   0.443
real(sp), parameter :: k12   =   0.895
real(sp), parameter :: swec  = 140.     ! minimum snow-water equivalent for full snow cover (Sandoval et al., 2024)
real(sp), parameter :: B0snw =   0.85   ! albedo of fresh snow

! local variables

real(sp) :: snowfall_day
real(sp) :: snowfall_night
real(sp) :: Bsnw

! ----
! snowfall

if (dmet%prec > 0.) then

  dmet%prday = dmet%prec * dmet%dayl / 24.
  dmet%prnight = dmet%prec - dmet%prday

  ! day

  snowfall_day = dmet%prday * (1. - frain(dmet%tday,pixel%Tt))
  
  ! night
  
  snowfall_night = dmet%prnight * (1. - frain(dmet%tnight,pixel%Tt))

  ! --

  dmet%snow = snowfall_day + snowfall_night
    
  dmet%swe = dmet%swe + dmet%snow

else 

  dmet%snow = 0.

end if
 
dmet%rain = max(dmet%prec - dmet%snow,0.)

if (dmet%rain < 0.) then
  write(0,*)'error in rain',dmet%rain,dmet%prec,dmet%snow
  stop
end if

! snowmelt - should separate into day and night parts

! write(0,*)'melt',dmet%HNpos,dmet%swe

dmet%melt = snowmelt(pixel%P,dmet%tday,dmet%swe,dmet%HNpos)

! snowpack accounting

dmet%swe = max(dmet%swe - dmet%melt,0.)

! write(0,*)'swe',dmet%HNpos,dmet%swe

! snow cover fraction

dmet%fsnow = dmet%swe / (swec + dmet%swe)  ! eqn. 27

! snowpack age and albedo

if (dmet%swe > 0.) then
  if (dmet%snow >= 3.) then
    dmet%asnow = 0
  else
    dmet%asnow = dmet%asnow + 1
  end if
  
!   if (dmet%asnow < 0 .or. dmet%asnow > 90) then 
!     write(0,*)dmet%asnow
!   end if
  
  dmet%asnow = min(dmet%asnow,90)

  Bsnw = (B0snw - k11) + k11 * exp(-k12 * real(dmet%asnow))  ! eqn. 26
  
  dmet%Bsw = B0 * (1. - dmet%fsnow) + (dmet%fsnow * Bsnw)   ! eqn. 25

else

  dmet%asnow = 0
  
  dmet%Bsw = B0

end if

if (dmet%Bsw < 0. .or. dmet%Bsw > 1.) then
  write(0,*)'albedo error',B0,Bsnw,dmet%fsnow
end if

! write(0,*)'snow',dmet%tday,dmet%tnight,dmet%prec,dmet%snow,dmet%melt,dmet%swe,dmet%fsnow,dmet%asnow

end subroutine snow

! ---------------------------------------

real(sp) function Tt(z,phi)

! function to assess the threshold temperature for snowfall based on a 50% probability
! using an analytical solution of Sandoval et al. (2024) eqn 19 for Tt at psnow of 0.5

use parametersmod, only : sp,dp

implicit none

! arguments

real(sp), intent(in) :: z     ! elevation (m)
real(dp), intent(in) :: phi   ! latitude (degrees)

! parameters                  line 576      | from paper Table 1
!                            splash.point.R | (NB the signs of k9 and k10 are inconsistent between paper and code)
real(sp), parameter :: k7  = 0.4710405934   ! -0.5827
real(sp), parameter :: k8  = 1.0473543991   !  1.319
real(sp), parameter :: k9  = 0.0004596581   ! -4.18e-4
real(sp), parameter :: k10 = 0.0110592101   ! -1.140e-2

real(sp) :: alat

! ----

alat = real(abs(phi))

Tt = (k7 + k9 * z + k10 * alat) / k8

end function Tt

! ---------------------------------------

real(sp) function frain(Tair,Tt)

! Notes: Sandoval et al. (2024) based on Kienzle (2008) suggest that the threshold temperature (Ttm) and 
! temperature range (Trm) for snowfall should vary seasonally, with Ttm smallest in winter and largest in summer
! and Trm being largest in spring and autumn, and smaller in summer and winter.
! As a model formulation based on month of the year as proposed by Kienzle cannot not be applicable outside of the 
! geographical domain for which it was developed (western Canada), we do not use the Sandoval et al. (2008) formulation here,
! but fix Ttm and Trm to single annual values. Ultimately, one could re-express Kienzle's functions in terms of 
! relationships between Ttm and the seasonal cycle of temperature and Trm and the derivative of the temperature seasonal cycle curve.

use parametersmod, only : sp

! arguments

real(sp), intent(in) :: Tair  ! air temperature (degC)
real(sp), intent(in) :: Tt    ! threshold daily maximum temperature for snowfall (degC)

! local variables

real(sp) :: Ttm  ! annual maximum air temperature when the probability of snowfall occurrence meets or exceeds 0.5
real(sp) :: Trm  ! temperature range over which snowfall occurs

real(sp) :: T    ! term inside the polynomial function

! ----
! as noted above the values of Ttm and Trm could vary seasonally but are fixed here

Ttm = Tt
Trm = 13.

T = (Tair - Ttm) / (1.4 * Trm)

frain = 5. * T**3 + 6.76 * T**2 + 3.19 * T + 0.5

frain = max(min(frain,1.),0.)

end function frain

! ---------------------------------------

real(sp) function snowmelt(P,temp,SWE,HNpos)

! Sandoval et al. (2024) eqn. 21, output units mm d-1

use parametersmod, only : sp
use physicsmod,    only : Econ

implicit none

! arguments

real(sp), intent(in) :: P     ! mean air pressure (Pa)
real(sp), intent(in) :: temp  ! temperature (degC)
real(sp), intent(in) :: SWE   ! snow-water equivalent (mm)
real(sp), intent(in) :: HNpos ! daytime accumulated net radiation (J m-2 d-1)

! parameters

real(sp), parameter :: pw = 999.8395  ! density of water at 0 degC (kg m-3) (Kell, 1975)
real(sp), parameter :: Lf =   3.34e5  ! latent heat of fusion of water (J kg-1) 

real(sp), parameter :: pwLf = pw * Lf ! product of the above two terms

! functions used

! real(sp) :: Econ ! energy-to-water equivalent conversion factor (m3 J-1)

! local variable

real(sp) :: psm   ! potential snowmelt
real(sp) :: Eswe  ! evaporating snowmelt
real(sp) :: HApos ! net radiation remaining after snowmelt

! ----

psm = 1000. * HNpos / pwLf

! write(0,*)'snowmelt',HNpos,psm

if (psm <= SWE) then

  snowmelt = psm

else

  snowmelt = SWE

  HApos = HNpos - SWE * pwLf / 1000.

  Eswe = min(snowmelt,HApos / Econ(P,temp) * 1000.)

end if

snowmelt = snowmelt - Eswe

end function snowmelt

! ---------------------------------------

end module snowmod
