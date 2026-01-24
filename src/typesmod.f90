module typesmod

! repository of derived types

use parametersmod, only : i2,sp,dp,nmos
use randomdistmod, only : randomstate

implicit none

public :: gridinfotype
public :: coordstype
public :: terraintype
public :: monclimatetype
public :: dayclimatetype
public :: soiltype
public :: metvars_in
public :: metvars_daily

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

  ! basic coordinates
  integer :: x       ! x position on 2D grid
  integer :: y       ! y position on 2D grid
  real(dp) :: lon    ! longitude (deg)
  real(dp) :: lat    ! latitude (deg)
  real(sp) :: phi    ! latitude (rad)
  ! integer, dimension(8) :: neighbors

  ! terrain characteristics
  real(sp) :: landf      ! part of gridcell that is ice-free land (fraction)
  real(sp) :: elv        ! elevation above sea level (m)
  real(sp) :: slope      ! terrain slope (m m-1)
  real(sp) :: aspect     ! slope azimuth (north = 0) (degrees)
  real(sp) :: cti        ! compound topographic index = log(upstream area / tan(local slope)) (index)
  real(sp) :: hand       ! height above nearest drainage (m)
  real(sp) :: elev_stdev ! standard deviation of elevation (m)
  real(sp) :: thickness  ! soil and regolith thickness (m)

  real(sp) :: srad       ! slope inclination (rad)
  real(sp) :: gamma      ! slope orientation (rad), 0 = S, values increasing clockwise

  ! meteorological variables
  real(sp) :: P      ! mean atmospheric pressure (Pa)
  real(sp) :: Ratm   ! relative atmospheric pressure
  real(sp) :: tcm    ! temperature of the coldest month
  real(sp) :: Pann   ! total annual precipitation
  real(sp) :: Pjj    ! precipitation equitability index
  real(sp) :: Tt     ! snow probability temperature (degC)
  real(sp) :: twm    ! temperature of the warmest month (degC)
  real(sp) :: gdd5   ! growing degree days on a 5-degree base (degC)
  real(sp) :: gdd0   ! growing degree days on a 0-degree base (degC)
  real(sp) :: awm    ! alpha of the warmest month (fraction)
  real(sp) :: aalpha ! mean annual alpha (fraction)
  real(sp) :: Nmelt  ! depletion curve shape parameter

  real(sp), dimension(40) :: dprec

  ! other
  integer :: wm
  integer(i2) :: biome
  
end type pixeltype

! ---

type terraintype
  real(sp) :: landf      ! part of gridcell that is ice-free land (fraction)
  real(sp) :: elv        ! elevation above sea level (m)
  real(sp) :: slope      ! terrain slope (m m-1)
  real(sp) :: aspect     ! slope azimuth (north = 0) (degrees)
  real(sp) :: cti        ! compound topographic index = log(upstream area / tan(local slope)) (index)
  real(sp) :: hand       ! height above nearest drainage (m)
  real(sp) :: elev_stdev ! standard deviation of elevation (m)
  real(sp) :: thickness  ! soil and regolith thickness (m)
end type terraintype

! ---

type monclimatetype
  real(sp) :: tmp
  real(sp) :: dtr
  real(sp) :: pre
  real(sp) :: cld
  real(sp) :: wnd
  real(sp) :: wet
end type monclimatetype

! ---

type dayclimatetype
  real(sp) :: tmin
  real(sp) :: tmax
  real(sp) :: cld
  real(sp) :: wnd
end type dayclimatetype

! ---

type soilcoordstype
  real(sp) :: zpos
  real(sp) :: dz
  real(sp), dimension(2) :: bnds
end type soilcoordstype

type soiltype
  real(sp) :: dz     ! actual layer thickness (cm)
  real(sp) :: Tsat
  real(sp) :: Ksat
  real(sp) :: whc
end type soiltype

! ---

type soilwatertype
  real(sp) :: whc    ! column-integrated water holding capacity (mm)
  real(sp) :: w      ! instantaneous soil water content (mm)
  real(sp) :: theta  ! volumetric water content (fraction)
  real(sp) :: psi    ! soil matric potential (kPa)
end type soilwatertype

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
  real(sp) :: rad0   ! mean daily top-of-the-atmosphere insolation (W m-2)
  real(sp) :: dayl   ! day length (h)
  real(sp) :: phi    ! latitude (rad)
  real(sp) :: delta  ! solar declination (rad)
end type solarpars

! ---

type airmasspars
  real(sp) :: mbar    ! daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(sp) :: mo      ! air mass at cosine zenith angle maximum
  real(sp) :: mc      ! air mass at cosine zenith angle medium
  real(sp) :: ml      ! air mass at cosine zenith angle bottom quarter range point
end type airmasspars

! ---

type metvars_in  ! structure for weather generator input

  ! monthly values

  real(sp) :: prec  ! total precipitation amount (mm)
  real(sp) :: wetf  ! fraction of days in the month with precipitation

  ! means-preserving interpolated estimate of daily values
  
  real(sp) :: tmin  ! minumum temperture (C)
  real(sp) :: tmax  ! maximum temperture (C)
  real(sp) :: cldf  ! cloud fraction (0=clear sky, 1=overcast) (fraction)
  real(sp) :: wind  ! wind speed (m/s)

  ! variables that carry over from one day to the next

  logical, dimension(2)  :: pday   ! precipitation occurrence on the preceding two days
  type(randomstate)      :: rndst  ! state of the random number generator
  real(sp), dimension(4) :: resid  ! previous day's weather residuals

end type metvars_in

! ---

type metvars_daily  ! structure for weather generator output (daily)

  ! basic output

  real(sp) :: delta     ! solar declination (rad)
  real(sp) :: dayl      ! daylength (h)

  real(sp) :: prec      ! total precipitation (mm)
  real(sp) :: tmin      ! minimum temperature (degC)
  real(sp) :: tmax      ! maximum temperature (degC)
  real(sp) :: cldf      ! mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind      ! wind speed (m s-1)

  real(sp) :: tdew      ! mean daily (24-hr) dewpoint temperature (degC)
  real(sp) :: tday      ! mean daytime temperature
  real(sp) :: tnight    ! mean nighttime temperature
  real(sp) :: wday      ! mean daytime windspeed
  real(sp) :: wnight    ! mean nighttime windspeed
  real(sp) :: prday     ! total daytime precipitation
  real(sp) :: prnight   ! total nighttime precipitation

  real(sp) :: rad0      ! top-of-atmosphere insolation (W m-2)
  real(sp) :: rdirect   ! downwelling shortwave radiation, direct beam component (kJ m-2 d-1)
  real(sp) :: rdiffuse  ! downwelling shortwave radiation, diffuse beam component (kJ m-2 d-1)
  real(sp) :: swrad     ! downwelling shortwave radiation, daily mean (W m-2)
  real(sp) :: lwday     ! downwelling longwave radiation, daytime mean (W m-2)
  real(sp) :: lwnight   ! downwelling longwave radiation, nighttime mean (W m-2)
  real(sp) :: dpet      ! daytime potential evapotranspiration (mm d-1)
  
  real(sp) :: Bsw       ! shortwave albedo
  ! real(sp) :: Bsnw      ! snow albedo
  
  real(sp) :: HNpos     ! daytime accumulated net radiation (MJ m-2 d-1)
  
  real(sp) :: aet       ! actual evapotranspiration (mm)
  real(sp) :: alpha     ! ratio of aet to pet
  real(sp) :: soilw     ! soil water content (mm)
  real(sp) :: relsat    ! relative saturation w/whc (0-1)
  
  real(sp) :: rain
  real(sp) :: snow
  real(sp) :: melt

  real(sp) :: swe       ! snow water equivalent (mm)
  real(sp) :: swe_max   ! maximum SWE this season (mm)
  real(sp) :: psnow     ! snow density ()
  real(sp) :: hsnow     ! snow depth (cm)
  real(sp) :: fsnow     ! snow areal cover (fraction)
  integer  :: asnow     ! snowpack age (days)
  
  ! diagnostic output, uncomment this section as needed

  ! real(sp) :: wind_bias
  ! real(sp) :: wind_intercept_bias
  ! real(sp) :: tmin_mn
  ! real(sp) :: tmin_sd
  ! real(sp) :: wind_mn
  ! real(sp) :: wind_sd

end type metvars_daily

type metvars_monthly

  real(sp) :: mpet
  real(sp) :: swe    ! snow water eqiv (mm)
  real(sp) :: snow   ! snowfall today (mm)
  real(sp) :: melt   ! snowmelt today (mm)
  real(sp) :: fsnow  ! snow fraction (0-1)
  real(sp) :: Bsw    ! shortwave albedo (0-1)
  real(sp) :: alpha
  real(sp) :: direct
  real(sp) :: diffuse

end type metvars_monthly

! ---

end module typesmod