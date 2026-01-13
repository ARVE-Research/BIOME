module snowmod

use parametersmod, only : sp, pi

implicit none

public :: Tt

contains

! ---------------------------------------

subroutine snow(pixel,dmet)


! calculate daily snow dynamics

use parametersmod, only : B0
use typesmod,      only : pixeltype,metvars_daily

implicit none

! arguments

type(pixeltype),     intent(in)    :: pixel
type(metvars_daily), intent(inout) :: dmet    ! daily meteorological variables

! parameters

real(sp), parameter :: k11   =   0.443
real(sp), parameter :: k12   =   0.895
! real(sp), parameter :: swec  = 140.     ! minimum snow-water equivalent for full snow cover (Sandoval et al., 2024)
real(sp), parameter :: B0snw =   0.85   ! albedo of fresh snow

! local variables

real(sp) :: snowfall_day
real(sp) :: snowfall_night
real(sp) :: Bsnw
real(sp) :: swe_old       ! track if accumulation or melt

! ----
! Store old SWE to determine if accumulation or melt occurred
swe_old = dmet%swe     

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

if (dmet%swe > 0. .and. dmet%tday >= 0.) then
  call snowmelt(pixel%P,dmet%tday,dmet%HNpos,dmet%swe,dmet%melt)  ! adjusts SWE
else 
  dmet%melt = 0.
end if

! ---------
! snow cover fraction

! dmet%fsnow = dmet%swe / (swec + dmet%swe)  ! eqn. 27

! snow cover fraction - Swenson & Lawrence (2012)

if (dmet%swe <= 0.) then
  ! No snow - reset everything
  dmet%swe_max = 0.
  dmet%fsnow = 0.
else if (dmet%swe > swe_old) then
  ! ACCUMULATION - assume full coverage
  dmet%fsnow = 1.0
  dmet%swe_max = dmet%swe
!   write(0,*)'ACCUM: swe=',dmet%swe,' fsnow=',dmet%fsnow,' swe_max=',dmet%swe_max  ! ADD THIS
else
  ! MELT - use depletion curve (eq. 4)
!   write(0,*)'MELT: swe_old=',swe_old,' swe=',dmet%swe,' swe_max=',dmet%swe_max,' Nmelt=',pixel%Nmelt  ! ADD THIS
  call calc_snow_cover_fraction(pixel, dmet)
!   write(0,*)'MELT result: fsnow=',dmet%fsnow  ! ADD THIS
end if

! ---------

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

! write(0,*)'snow',dmet%tday,dmet%tnight,dmet%prec,dmet%snow,dmet%melt,dmet%swe,dmet%fsnow,dmet%asnow,dmet%Bsw

end subroutine snow

! ---------------------------------------
! UNDER CONSTRUCTION ---- NEW SUBROUTINE FOR SNOW COVER FRACTION
subroutine calc_snow_cover_fraction(pixel, dmet)

! Calculate snow cover fraction during melt using Swenson & Lawrence (2012) eq. 4

use parametersmod, only : sp, pi
use typesmod,      only : pixeltype, metvars_daily

implicit none

type(pixeltype),     intent(in)    :: pixel
type(metvars_daily), intent(inout) :: dmet

real(sp) :: rel_swe     ! relative SWE (W/Wmax)

! ----

if (dmet%swe_max > 0.) then
  rel_swe = dmet%swe / dmet%swe_max
  rel_swe = max(0., min(1., rel_swe))  ! constrain to [0,1]
  
  ! Inverse cosine depletion curve (eq. 4)
  dmet%fsnow = 1. - (1./pi * acos(2. * rel_swe - 1.))**pixel%Nmelt
  
  dmet%fsnow = max(0., min(1., dmet%fsnow))
else
  dmet%fsnow = 0.
end if

end subroutine calc_snow_cover_fraction
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

subroutine snowmelt(P,temp,HNpos,SWE,melt)

! Sandoval et al. (2024) eqn. 21, output units mm d-1

use parametersmod, only : sp
use physicsmod,    only : Econ

implicit none

! arguments

real(sp), intent(in)    :: P      ! mean air pressure (Pa)
real(sp), intent(in)    :: temp   ! air temperature (degC)
real(sp), intent(in)    :: HNpos  ! 24 hour integrated positive net radiation (J m-2 d-1)
real(sp), intent(inout) :: SWE    ! snow-water equivalent (mm)
real(sp), intent(out)   :: melt   ! snowmelt adjusted for sublimation (if any) (mm)

! parameters

real(sp), parameter :: pw = 999.8395  ! density of water at 0 degC (kg m-3) (Kell, 1975)
real(sp), parameter :: Lf =   3.3355e5  ! latent heat of fusion of water (J kg-1) 

real(sp), parameter :: pwLf = pw * Lf ! product of the above two terms

! functions used

! real(sp) :: Econ ! energy-to-water equivalent conversion factor (m3 J-1)

! local variable

real(sp) :: psm   ! potential snowmelt
real(sp) :: Eswe  ! evaporating snowmelt
real(sp) :: HApos ! net radiation remaining after snowmelt

! ----

psm = HNpos / pwLf * 1000.  ! units (J m-2 d-1) / (kg m-3 J kg-1) = m * 1000 = mm

if (psm <= SWE) then

  melt = psm
  
  Eswe = 0.
  
else

  melt = SWE
  
  HApos = HNpos - SWE * pwLf / 1000.

  Eswe = min(melt,HApos / Econ(P,temp) * 1000.)

end if

swe = swe - melt

melt = melt - Eswe  

end subroutine snowmelt

! ---------------------------------------

end module snowmod
