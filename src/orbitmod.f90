module orbitmod

implicit none

! module subprograms

public  :: getorbitpars

contains

! --------------------------------------------------------------------------------------------

subroutine getorbitpars(yrbp,orbit)

! estimate the orbital parameters using linear interpolation from the lookup table 

use parametersmod, only : dp,pir => pir_dp
use typesmod,      only : orbitpars
use utilitymod,    only : pos,spline,cubspline,angleinterp
use netcdf
use errormod,      only : ncstat,netcdf_err

implicit none

integer,         intent(in)  :: yrbp
type(orbitpars), intent(out) :: orbit

integer :: ncid
integer :: dimid
integer :: varid

integer :: nrst
integer :: i0
integer :: i1

real(dp) :: t0
real(dp) :: t1

integer :: tlen

real(dp) :: yrj2

real(dp), allocatable, dimension(:) :: time

integer, parameter :: spw = 11  ! time window over which to evaluate the spline

real(dp) :: wgt

real(dp), dimension(2)   :: perh
real(dp), dimension(spw) :: ecc
real(dp), dimension(spw) :: pre
real(dp), dimension(spw) :: xob

real(dp), dimension(spw) :: splineobj

character(50), parameter :: orbitparsfile = 'share/La2004.nc'

! --------------------------------
! convert yrbp into yr since J2000

yrj2 = (-50. - real(yrbp))

! find the bracketing period

ncstat = nf90_open(orbitparsfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ncid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(time(tlen))

ncstat = nf90_inq_varid(ncid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,time)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! find the time nearest to the input time in the orbital parameters file

nrst = pos(time,yrj2)

! define a window +/-5000 yr around the nearest point (this does not need to be more precise)

i0 = nrst - 5
i1 = nrst + 5

t0 = time(i0)
t1 = time(i1)

! ----
! read the precalculated Laskar 2004 orbital parameters at 1000-year intevals

ncstat = nf90_inq_varid(ncid,'ecc',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ecc,start=[i0],count=[spw])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ncid,'pre',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,pre,start=[i0],count=[spw])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ncid,'xob',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,xob,start=[i0],count=[spw])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! create the spline model and evaluate the spline for each variable
! the boundary conditions are set to zero (it won't matter with such a broad time window)

call spline(time(i0:i1),ecc,0._dp,0._dp,splineobj)

orbit%ecc = cubspline(time(i0:i1),ecc,splineobj,yrj2)

call spline(time(i0:i1),pre,0._dp,0._dp,splineobj)

orbit%pre = cubspline(time(i0:i1),pre,splineobj,yrj2)

call spline(time(i0:i1),xob,0._dp,0._dp,splineobj)

orbit%xob = cubspline(time(i0:i1),xob,splineobj,yrj2)

! ----
! longitude of precession changes linearly and because of the discontinuity should not be interpolated with the spline
! but does require interpolation between angles, so we use a special function (angleinterp) for that

i0 = nrst

t0 = time(i0)

if (yrj2 > t0) then

  i1 = i0 + 1

else

  i1 = i0
  i0 = i0 - 1

end if

ncstat = nf90_inq_varid(ncid,'perh',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,perh,start=[i0],count=[2])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

t0 = time(i0)
t1 = time(i1)

wgt = 1. - (t1 - yrj2) / (t1 - t0)

orbit%perh = angleinterp(perh(1),perh(2),wgt)

end subroutine getorbitpars

! ------------------------------------

end module orbitmod
