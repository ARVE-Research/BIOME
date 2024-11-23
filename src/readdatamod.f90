module readdatamod

use parametersmod, only : sp

implicit none

public  :: readcoords
public  :: readterrain
public  :: readclimate
public  :: readsoil
private :: getvar_i2
private :: getvar_sp

contains

! -----------------------------------------------------

subroutine readcoords(coordsfile,gridinfo,coords)

use netcdf
use errormod,  only : ncstat,netcdf_err
use typesmod,  only : gridinfotype,coordstype

implicit none

! arguments

character(*),                     intent(in)  :: coordsfile
type(gridinfotype),               intent(in)  :: gridinfo
type(coordstype), dimension(:,:), intent(out) :: coords

! local variables

integer :: ncid
integer :: varid

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

! -----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

write(0,*)'readcoords',srtx,cntx,srty,cnty

ncstat = nf90_open(coordsfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

if (gridinfo%isproj) then

  ncstat = nf90_inq_varid(ncid,'x',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords(:,1)%xcoord,start=[srtx],count=[cntx])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ncid,'y',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords(1,:)%ycoord,start=[srty],count=[cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ncid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords%geolon,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ncid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords%geolat,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else

  ncstat = nf90_inq_varid(ncid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords(:,1)%xcoord,start=[srtx],count=[cntx])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  coords%geolon = coords%xcoord

  ncstat = nf90_inq_varid(ncid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords(1,:)%ycoord,start=[srty],count=[cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  coords%geolat = coords%ycoord

end if
  
ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readcoords

! -----------------------------------------------------

subroutine readterrain(terrainfile,gridinfo,terrain)

use netcdf
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype,terraintype
use parametersmod, only : i2,stdout,stderr

implicit none

! arguments

character(*),                      intent(in)  :: terrainfile
type(gridinfotype),                intent(in)  :: gridinfo
type(terraintype), dimension(:,:), intent(out) :: terrain

! local variables

integer :: ncid
integer :: varid

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

real(sp), dimension(2) :: actual_range

integer(i2) :: missing_value

! -----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

ncstat = nf90_open(terrainfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ncid,'elv',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,terrain%elv,start=[srtx,srty],count=[cntx,cnty])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

actual_range(1) = minval(terrain%elv,mask = terrain%elv /= missing_value)
actual_range(2) = maxval(terrain%elv,mask = terrain%elv /= missing_value)

write(stderr,*)'reading elevation',actual_range
  
ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readterrain

! -----------------------------------------------------

subroutine readclimate(climatefile,gridinfo,climate)

use netcdf
use errormod,  only : ncstat,netcdf_err
use typesmod,  only : gridinfotype,climatetype

implicit none

! arguments

character(*),                        intent(in)  :: climatefile
type(gridinfotype),                  intent(in)  :: gridinfo
type(climatetype), dimension(:,:,:), intent(out) :: climate

! local variables

integer :: ncid

! -----

ncstat = nf90_open(climatefile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

call getvar_i2(ncid,'tmp',gridinfo,climate%tmp)
call getvar_i2(ncid,'dtr',gridinfo,climate%dtr)
call getvar_i2(ncid,'pre',gridinfo,climate%pre)
call getvar_i2(ncid,'wet',gridinfo,climate%wet)
call getvar_i2(ncid,'cld',gridinfo,climate%cld)
call getvar_i2(ncid,'wnd',gridinfo,climate%wnd)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readclimate

! -----------------------------------------------------

subroutine readsoil(soilfile,gridinfo,soil)

use netcdf
use errormod,  only : ncstat,netcdf_err
use typesmod,  only : gridinfotype,soiltype

implicit none

! arguments

character(*),                     intent(in)  :: soilfile
type(gridinfotype),               intent(in)  :: gridinfo
type(soiltype), dimension(:,:,:), intent(out) :: soil

! local variables

integer :: ncid

! -----

ncstat = nf90_open(soilfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(0,*)soilfile

call getvar_sp(ncid,'Tsat',gridinfo,soil%Tsat)
call getvar_sp(ncid,'Ksat',gridinfo,soil%Ksat)
call getvar_sp(ncid,'whc',gridinfo,soil%whc)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readsoil

! -----------------------------------------------------

subroutine getvar_i2(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i2,sp,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype

implicit none

integer,                    intent(in)  :: ncid
character(*),               intent(in)  :: varname
type(gridinfotype),         intent(in)  :: gridinfo
real(sp), dimension(:,:,:), intent(out) :: ovar

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

integer :: zlen

integer :: varid

integer(i2), allocatable, dimension(:,:,:) :: ivar

real(sp) :: scale_factor
real(sp) :: add_offset

integer(i2) :: missing_value

real(sp), dimension(2) :: actual_range

! ----

zlen = size(ovar,dim=3)

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

allocate(ivar(cntx,cnty,zlen))

ovar = rmissing

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ivar,start=[srtx,srty,1],count=[cntx,cnty,zlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'add_offset',add_offset)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (ivar /= missing_value) ovar = real(ivar) * scale_factor + add_offset

actual_range(1) = minval(ovar,mask = ivar /= missing_value)
actual_range(2) = maxval(ovar,mask = ivar /= missing_value)

write(stderr,*)'reading ',trim(varname),actual_range

end subroutine getvar_i2

! -----------------------------------------------------

subroutine getvar_sp(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i2,sp,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype

implicit none

integer,                    intent(in)  :: ncid
character(*),               intent(in)  :: varname
type(gridinfotype),         intent(in)  :: gridinfo
real(sp), dimension(:,:,:), intent(out) :: ovar

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

integer :: varid

real(sp) :: missing_value

real(sp), dimension(2) :: actual_range

! ----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

write(0,*)'getvar sp',gridinfo%cntx,gridinfo%cnty
write(0,*)size(ovar,dim=1),size(ovar,dim=2),size(ovar,dim=3)

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ovar,start=[srtx,srty,1],count=[cntx,cnty,6])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! actual_range(1) = minval(ovar,mask = ovar /= missing_value)
! actual_range(2) = maxval(ovar,mask = ovar /= missing_value)

write(stderr,*)'reading ',trim(varname),actual_range

where (ovar == missing_value) ovar = rmissing

end subroutine getvar_sp

! -----------------------------------------------------

end module readdatamod

