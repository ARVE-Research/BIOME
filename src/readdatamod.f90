module readdatamod

use parametersmod, only : sp

implicit none

public  :: readcoords
public  :: readterrain
public  :: readclimate
public  :: readsoil
private :: getvar_i2
private :: getvar_sp

interface getvar_sp
  module procedure getvar_sp_2d
  module procedure getvar_sp_3d
end interface getvar_sp

interface getvar_i2
  module procedure getvar_i2_2d
  module procedure getvar_i2_3d
end interface getvar_i2

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

integer :: x
integer :: y

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

  do y = 1,cnty
    coords(:,y)%geolon = coords(:,1)%xcoord
  end do

  ncstat = nf90_inq_varid(ncid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_get_var(ncid,varid,coords(1,:)%ycoord,start=[srty],count=[cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  do x = 1,cntx
    coords(x,:)%geolat = coords(1,:)%ycoord
  end do

end if

! do y = 1,cnty
!   do x = 1,cntx
!     write(0,*)x,y,coords(x,y)%geolon,coords(x,y)%geolat
!   end do
! end do
! 
! stop

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

ncstat = nf90_open(terrainfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

call getvar_i2(ncid,'landf',gridinfo,terrain%landf)
call getvar_i2(ncid,'elv',gridinfo,terrain%elv)
call getvar_i2(ncid,'slope',gridinfo,terrain%slope)
call getvar_i2(ncid,'aspect',gridinfo,terrain%aspect)
call getvar_i2(ncid,'cti',gridinfo,terrain%cti)
call getvar_i2(ncid,'hand',gridinfo,terrain%hand)
call getvar_i2(ncid,'elev_stdev',gridinfo,terrain%elev_stdev)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readterrain

! -----------------------------------------------------

subroutine readclimate(climatefile,gridinfo,climate)

use netcdf
use errormod,  only : ncstat,netcdf_err
use typesmod,  only : gridinfotype,monclimatetype

implicit none

! arguments

character(*),                           intent(in)  :: climatefile
type(gridinfotype),                     intent(in)  :: gridinfo
type(monclimatetype), dimension(:,:,:), intent(out) :: climate

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

subroutine readsoil(soilfile,gridinfo,terrain,soilcoords,soilinput)

use netcdf
use errormod,  only : ncstat,netcdf_err
use typesmod,  only : gridinfotype,soilinputtype,soilcoordstype,terraintype

implicit none

! arguments

character(*),                           intent(in)    :: soilfile
type(gridinfotype),                     intent(in)    :: gridinfo
type(soilcoordstype), dimension(:),     intent(out)   :: soilcoords
type(terraintype),    dimension(:,:),   intent(inout) :: terrain
type(soilinputtype),  dimension(:,:,:), intent(out)   :: soilinput

! local variables

integer :: ncid
integer :: varid

integer :: l
integer :: nl

real(sp), dimension(2,size(soilcoords)) :: layer_bnds

! -----

nl = size(soilcoords)

ncstat = nf90_open(soilfile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ---

ncstat = nf90_inq_varid(ncid,'depth',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,soilcoords%zpos)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ncid,'dz',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,soilcoords%dz)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ncid,'layer_bnds',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,layer_bnds)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do l = 1,nl
  soilcoords(l)%layer_bnds = layer_bnds(:,l)
end do

! ---

call getvar_i1(ncid,'USDA',gridinfo,terrain%USDA)
call getvar_sp(ncid,'thickness',gridinfo,terrain%thickness)
call getvar_i2(ncid,'sand',gridinfo,soilinput%sand)
call getvar_i2(ncid,'clay',gridinfo,soilinput%clay)
call getvar_i2(ncid,'cfvo',gridinfo,soilinput%cfvo)
call getvar_i2(ncid,'soc',gridinfo,soilinput%soc)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine readsoil

! -----------------------------------------------------

subroutine getvar_i2_3d(ncid,varname,gridinfo,ovar)

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
if (ncstat == nf90_enotatt) then
  scale_factor = 1.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

ncstat = nf90_get_att(ncid,varid,'add_offset',add_offset)
if (ncstat == nf90_enotatt) then
  add_offset = 0.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

where (ivar /= missing_value) ovar = real(ivar) * scale_factor + add_offset

actual_range(1) = minval(ovar,mask = ivar /= missing_value)
actual_range(2) = maxval(ovar,mask = ivar /= missing_value)

write(stderr,*)'reading ',trim(varname),actual_range

end subroutine getvar_i2_3d

! -----------------------------------------------------

subroutine getvar_i2_2d(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i2,sp,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype

implicit none

integer,                  intent(in)  :: ncid
character(*),             intent(in)  :: varname
type(gridinfotype),       intent(in)  :: gridinfo
real(sp), dimension(:,:), intent(out) :: ovar

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

integer :: varid

integer(i2), allocatable, dimension(:,:) :: ivar

real(sp) :: scale_factor
real(sp) :: add_offset

integer(i2) :: missing_value

real(sp), dimension(2) :: actual_range

! ----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

allocate(ivar(cntx,cnty))

ovar = rmissing

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ivar,start=[srtx,srty],count=[cntx,cnty])
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

write(stderr,*)'read ',trim(varname),actual_range

end subroutine getvar_i2_2d

! -----------------------------------------------------

subroutine getvar_i1(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i1,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype

implicit none

integer,                     intent(in)  :: ncid
character(*),                intent(in)  :: varname
type(gridinfotype),          intent(in)  :: gridinfo
integer(i1), dimension(:,:), intent(out) :: ovar

integer :: srtx
integer :: cntx
integer :: srty
integer :: cnty

integer :: varid

real(sp) :: scale_factor
real(sp) :: add_offset

integer(i1) :: missing_value

integer(i1), dimension(2) :: actual_range

! ----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ovar,start=[srtx,srty],count=[cntx,cnty])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

actual_range(1) = minval(ovar,mask = ovar /= missing_value)
actual_range(2) = maxval(ovar,mask = ovar /= missing_value)

write(stderr,*)'read ',trim(varname),actual_range

end subroutine getvar_i1

! -----------------------------------------------------

subroutine getvar_sp_3d(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i2,sp,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype
use ieee_arithmetic

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

integer :: nl

! ----

srtx = gridinfo%srtx
cntx = gridinfo%cntx
srty = gridinfo%srty
cnty = gridinfo%cnty

nl = size(ovar,dim=3)

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ovar,start=[srtx,srty,1],count=[cntx,cnty,nl])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (ieee_is_nan(ovar)) ovar = missing_value

actual_range(1) = minval(ovar,mask = ovar /= missing_value)
actual_range(2) = maxval(ovar,mask = ovar /= missing_value)

write(stderr,*)'reading ',trim(varname),actual_range

where (ovar == missing_value) ovar = rmissing

end subroutine getvar_sp_3d

! -----------------------------------------------------

subroutine getvar_sp_2d(ncid,varname,gridinfo,ovar)

use netcdf
use parametersmod, only : i2,sp,stdout,stderr,rmissing
use errormod,      only : ncstat,netcdf_err
use typesmod,      only : gridinfotype
use ieee_arithmetic

implicit none

integer,                  intent(in)  :: ncid
character(*),             intent(in)  :: varname
type(gridinfotype),       intent(in)  :: gridinfo
real(sp), dimension(:,:), intent(out) :: ovar

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

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ovar,start=[srtx,srty],count=[cntx,cnty])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (ieee_is_nan(ovar)) ovar = missing_value

actual_range(1) = minval(ovar,mask = ovar /= missing_value)
actual_range(2) = maxval(ovar,mask = ovar /= missing_value)

write(stderr,*)'reading ',trim(varname),actual_range

where (ovar == missing_value) ovar = rmissing

end subroutine getvar_sp_2d

! -----------------------------------------------------

end module readdatamod

