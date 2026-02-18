module tsoutputmod

implicit none

public :: writedailymetvars

contains

! ----------------------------------------------------------------------------------------------------------------

subroutine writedailymetvars(m,d,dmet)

use parametersmod, only : dmetfile_unit
use typesmod, only : metvars_daily

integer, intent(in)             :: m      ! month number
integer, intent(in)             :: d      ! day number
type(metvars_daily), intent(in) :: dmet

  write(dmetfile_unit, *) m, d,       &
    dmet%tday,   &   !  3 - daytime temperature (C)
    dmet%tnight, &   !  4 - nighttime temperature (C)
    dmet%prec,   &   !  5 - precipitation (mm)
    dmet%snow,   &   !  6 - snowfall (mm)
    dmet%melt,   &   !  7 - snowmelt (mm)
    dmet%swe,    &   !  8 - snow water equivalent (mm)
    dmet%fsnow,  &   !  9 - fractional snow cover (0-1)
    dmet%asnow,  &   ! 10 - snow albedo
    dmet%Bsw,    &   ! 11 - shortwave radiation albedo
    dmet%alpha,  &   ! 12 - AET/PET ratio (0-1)
    dmet%tdew,   &   ! 13 - dewpoint temperature (C)
    dmet%aet,    &   ! 14 - actual evapotranspiration (mm)
    dmet%soilw,  &   ! 15 - soil water content (mm)
    dmet%relsat, &   ! 16 - relative saturation w/whc (0-1)
    dmet%dpet,   &   ! 17 - daily PET (mm)
    dmet%swrad,  &   ! 18 - total surface shortwave radiation (W m-2)
    dmet%lw_rad, &   ! 19 - net longwave, Sandoval method (W m-2)
    dmet%sunf,   &   ! 20 - sunshine fraction (0-1)
    dmet%tday,   &   ! 21 - daytime temperature (C)
    dmet%lwday,  &   ! 22 - daytime longwave radiation (W m-2)
    dmet%tnight, &   ! 23 - nighttime temperature (C)
    dmet%lwnight,&   ! 24 - nighttime longwave radiation (W m-2)
    dmet%hour_sw, &   ! 25 - hour angle of sunset (hours)
    dmet%hour_net, &   ! 26 - hour angle of net radiation crossover (hours)
    dmet%HNpos,  &   ! 27 - positive (daytime) net radiation (J m-2 d-1)
    dmet%HNneg,  &   ! 28 - negative (nighttime) net radiation (J m-2 d-1)
    dmet%tmin,   &   ! 29 - minimum temperature (C)
    dmet%tmax        ! 30 - maximum temperature (C)

end subroutine writedailymetvars

! ----------------------------------------------------------------------------------------------------------------

end module tsoutputmod
