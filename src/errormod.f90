module errormod

implicit none

integer :: ncstat

contains

! -------------------------

subroutine netcdf_err(ncstat)

use netcdf,        only : nf90_strerror
use parametersmod, only : stdout

implicit none

integer, intent(in) :: ncstat

write(stdout,'(a,i5,a,a)')' NetCDF error ',ncstat,' encountered: ',trim(nf90_strerror(ncstat))
stop

end subroutine netcdf_err

! -------------------------

end module errormod
