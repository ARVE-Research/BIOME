module typesmod

! repository of derived types

use parametersmod, only : sp,dp
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
  ! logical :: valid
  integer :: x
  integer :: y
  ! integer, dimension(8) :: neighbors
end type pixeltype

! ---

type terraintype
  real(sp) :: elv
  real(sp) :: slope
  real(sp) :: aspect
  real(sp) :: cti
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

  logical, dimension(2)  :: pday   ! precipitation status: true if the day was a rain day
  type(randomstate)      :: rndst  ! state of the random number generator
  real(sp), dimension(4) :: resid  ! previous day's weather residuals

end type metvars_in

! ---

type metvars_out  ! structure for weather generator output (daily)

  real(sp) :: prec  ! total precipitation (mm)
  real(sp) :: tmin  ! minimum temperature (degC)
  real(sp) :: tmax  ! maximum temperature (degC)
  real(sp) :: cldf  ! mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind  ! wind speed (m s-1)

  logical,  dimension(2) :: pday   ! precipitation state
  type(randomstate)      :: rndst  ! state of the random number generator, 15 elements
  real(sp), dimension(4) :: resid  ! previous day's weather residuals
  real(sp)               :: wind_bias
  real(sp)               :: wind_intercept_bias
  real(sp)               :: tmin_mn
  real(sp)               :: tmin_sd
  real(sp)               :: wind_mn
  real(sp)               :: wind_sd
  real(sp), dimension(4) :: unorm

end type metvars_out

! ---

end module typesmod