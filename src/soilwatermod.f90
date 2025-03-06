module soilwatermod

implicit none

contains

! ----------------------------------------------------------------------------------------------------------------

subroutine calcwhc(terrain,soilcoords,soil,soilw)

use parametersmod, only : sp
use typesmod,      only : soilcoordstype,soiltype,soilwatertype,terraintype
use utilitymod,    only : pos

implicit none

! arguments

type(terraintype),                  intent(in) :: terrain
type(soilcoordstype), dimension(:), intent(in) :: soilcoords
type(soiltype),       dimension(:), intent(inout) :: soil
type(soilwatertype),                intent(out)   :: soilw

! local variables

real(sp) :: depthT
real(sp), dimension(size(soilcoords)) :: lb

integer :: lbot
integer :: nl
integer :: l

real(sp) :: excess

! ----

nl = size(soilcoords)

soil%dz = soilcoords%dz

depthT = terrain%thickness * 100.  ! convert m to cm

lb = soilcoords%bnds(2)

! find layer boundary closest to total soil depth

lbot = pos(soilcoords%bnds(2),depthT)

if (depthT < maxval(lb) .and. lb(lbot) < depthT) lbot = lbot + 1

! zero out thickness below soil depth (for later: prescribe this to bedrock properties)

if (lbot < nl) soil(lbot+1:)%dz = 0.  

if (depthT < lb(lbot)) then
  excess = depthT - sum(soilcoords(1:lbot-1)%dz)
  soil(lbot)%dz = excess
end if

! calculate the column-integrated soil water holding capacity

soilw%whc = sum(soil%whc * 10. * soil%dz)  ! need to convert cm to mm

! diagnostic output
! if (depthT < 200.) then
! 
!   write(0,*)'soil thickness:',depthT
!   ! write(0,*)soilcoords%dz
!   
!   do l = 1,nl
!     write(0,*)l,soilcoords(l)%bnds(2),soilcoords(l)%dz,soil(l)%dz,soil(l)%whc
!   end do
! 
!   write(0,*)'total whc:     ',soilw%whc
! 
! end if

end subroutine calcwhc

! ----------------------------------------------------------------------------------------------------------------

subroutine soilwater(dmet,soilw)

use parametersmod, only : sp
use typesmod

implicit none

! arguments

type(metvars_daily), intent(inout) :: dmet
type(soilwatertype), intent(inout) :: soilw

!parameters

real(sp), parameter :: etmax = 1.  ! maximum equilibrium evapotranspiration rate (mm h-1)

!local variables

real(sp) :: Emd     ! Emax: maximum evapotranspiration rate from saturated soils(mm d-1)
real(sp) :: supply  ! water supply (mm) 
real(sp) :: demand  ! water demand (mm)

! -----

dmet%rain = dmet%prec
dmet%melt = 0.

! ----

Emd = etmax * dmet%dayl

if (soilw%whc > 0. .and. soilw%w > 0.) then
  supply = Emd * soilw%w / soilw%whc  ! min(Emd,soilw%w)  ! alternative formulation
else
  supply = 0.
end if

demand = dmet%dpet

dmet%aet = min(supply,demand)

! zero out low values to avoid underflow

if (soilw%w < 1.e-5) soilw%w = 0.

soilw%w = min(soilw%w + (dmet%rain + dmet%melt - dmet%aet),soilw%whc)

if (demand > 0.) then
  dmet%alpha = dmet%aet / dmet%dpet
else
  dmet%alpha = 1.
end if

! write(0,*)'aet',dmet%dayl,supply,demand,dmet%aet,dmet%alpha

end subroutine soilwater

! ----------------------------------------------------------------------------------------------------------------

end module soilwatermod