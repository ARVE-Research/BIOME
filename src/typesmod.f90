module typesmod

! repository of derived types

use parametersmod, only : sp,dp

implicit none

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

end module typesmod