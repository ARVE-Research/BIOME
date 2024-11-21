module typesmod

! repository of derived types

use parametersmod, only : sp,dp,nmos
use randomdistmod, only : randomstate

implicit none

public :: gridinfotype
public :: coordstype
public :: terraintype
public :: climatetype
public :: soiltype
public :: metvars_in
public :: metvars_out

! ---

type gridinfotype
  logical  :: isproj
  real(dp) :: xmin
  real(dp) :: xmax
  real(dp) :: ymin
  real(dp) :: ymax
  integer  :: srtx
  integer  :: cntx
  integer  :: srty
  integer  :: cnty
end type gridinfotype

! ---

type coordstype
  real(dp) :: xcoord
  real(dp) :: ycoord
  real(dp) :: geolon
  real(dp) :: geolat
end type coordstype

! ---

type pixeltype
  integer :: x
  integer :: y
  real(dp) :: lon
  real(dp) :: lat
  real(sp) :: Ratm ! relative atmospheric pressure
  real(sp) :: tcm  ! temperature of the coldest month
  real(sp) :: Pjj  ! precipitation equitability index
  ! logical :: valid
  ! integer, dimension(8) :: neighbors
end type pixeltype

! ---

type terraintype
  real(sp) :: elv
  real(sp) :: slope
  real(sp) :: aspect
  real(sp) :: cti
  real(sp) :: landf
  real(sp) :: waterf
  real(sp) :: icef
end type terraintype

! ---

type climatetype
  real(sp) :: tmp
  real(sp) :: dtr
  real(sp) :: pre
  real(sp) :: cld
  real(sp) :: wnd
  real(sp) :: wet
end type climatetype

! ---

type soiltype
  real(sp) :: Tsat
  real(sp) :: Ksat
  real(sp) :: whc
end type soiltype

! ---

type calendartype
  character(100)          :: calname  ! name of the calendar (optional)
  integer                 :: yrbp     ! year BP (1950 CE)
  integer                 :: ndyr     ! number of days in the year
  real(dp)                :: veqday   ! day of year of the vernal equinox, predefined
  real(dp), dimension(12) :: ndmr     ! number of days per month (real)
  integer,  dimension(12) :: ndmi     ! number of days per month (integer)
end type calendartype

type calendarstype
  type(calendartype) :: noleap
  type(calendartype) :: leapyr
end type calendarstype  

! ---

type orbitpars
  real(dp) :: ecc   ! eccentricity parameter (unitless)
  real(dp) :: pre   ! precession parameter (unitless)
  real(dp) :: perh  ! longitude of perihelion (degrees, heliocentric)
  real(dp) :: xob   ! obliquity (tilt) (degrees)
end type orbitpars

! ---

type solarpars
  real(sp) :: rad0     ! mean daily top-of-the-atmosphere insolation (W m-2)
  real(sp) :: dayl     ! day length (h)
  real(sp) :: delta    ! solar declination (rad)
end type solarpars

! ---

type airmasspars
  real(sp) :: Ratm    ! relative atmospheric pressure 1=sea level
  real(sp) :: mbar    ! daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(sp) :: mo      ! air mass at cosine zenith angle maximum
  real(sp) :: mc      ! air mass at cosine zenith angle medium
  real(sp) :: ml      ! air mass at cosine zenith angle bottom quarter range point
end type airmasspars

! ---

type metvars_in  ! structure for weather generator input

  ! monthly means

  real(sp) :: prec  ! total precipitation amount (mm)
  real(sp) :: wetd  ! number of days in the month with precipitation
  real(sp) :: wetf  ! fraction of days in the month with precipitation

  ! means-preserving interpolated estimate of daily values
  
  real(sp) :: tmin  ! minumum temperture (C)
  real(sp) :: tmax  ! maximum temperture (C)
  real(sp) :: cldf  ! cloud fraction (0=clear sky, 1=overcast) (fraction)
  real(sp) :: wind  ! wind speed (m/s)
  
  ! variables that carry over from one day to the next

  logical, dimension(2)  :: pday   ! precipitation status: true if the day was a rain day
  type(randomstate)      :: rndst  ! state of the random number generator
  real(sp), dimension(4) :: resid  ! previous day's weather residuals

end type metvars_in

! ---

type metvars_out  ! structure for weather generator output (daily)

  ! basic output

  real(sp) :: dayl  ! daylength (h)

  real(sp) :: prec  ! total precipitation (mm)
  real(sp) :: tmin  ! minimum temperature (degC)
  real(sp) :: tmax  ! maximum temperature (degC)
  real(sp) :: cldf  ! mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind  ! wind speed (m s-1)

  real(sp) :: tdew    ! estimated predawn dewpoint
  real(sp) :: tday    ! mean daytime temperature
  real(sp) :: tnight  ! mean nighttime temperature
  real(sp) :: wday    ! mean daytime windspeed
  real(sp) :: wnight  ! mean nighttime windspeed

  real(sp) :: rad0      ! top-of-atmosphere insolation (W m-2)
  real(sp) :: rdirect   ! surface downwelling radiation, direct beam component (kJ m-2 d-1)
  real(sp) :: rdiffuse  ! surface downwelling radiation, diffuse beam component (kJ m-2 d-1)
  real(sp) :: pet       ! daytime potential evapotranspiration (mm d-1)
  
  ! diagnostic output

!   real(sp) :: wind_bias
!   real(sp) :: wind_intercept_bias
!   real(sp) :: tmin_mn
!   real(sp) :: tmin_sd
!   real(sp) :: wind_mn
!   real(sp) :: wind_sd

end type metvars_out

! ---

end module typesmod