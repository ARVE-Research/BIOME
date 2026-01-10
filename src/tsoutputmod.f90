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

write(dmetfile_unit,*)m,d,dmet%tday,dmet%tnight,dmet%prec,dmet%snow,dmet%melt,dmet%swe,dmet%fsnow,dmet%asnow,dmet%Bsw

! write(dmetfile_unit,*)m,d,dmet

end subroutine writedailymetvars

! ----------------------------------------------------------------------------------------------------------------

end module tsoutputmod
