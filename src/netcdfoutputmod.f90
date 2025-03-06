module netcdfoutputmod

implicit none

contains

! -----------------------------------------------------

subroutine genoutputfile(jobfile,outfile,gridinfo,coords,ofid)

use netcdf
use errormod,      only : ncstat,netcdf_err
use parametersmod, only : i2,dp,imissing,rmissing,nmos,npft
use typesmod,      only : gridinfotype,coordstype

implicit none

character(*),                     intent(in)  :: jobfile
character(*),                     intent(in)  :: outfile
type(gridinfotype),               intent(in)  :: gridinfo
type(coordstype), dimension(:,:), intent(in)  :: coords
integer,                          intent(out) :: ofid

integer, dimension(4) :: dimids
integer, dimension(4) :: chunks

character(8)  :: today
character(10) :: now

integer :: dimid
integer :: varid

integer :: m

integer, dimension(nmos), parameter :: month = [(m,m=1,12)]
integer, dimension(npft), parameter :: PFT   = [(m,m=1,13)]

real(dp), dimension(2) :: xrange = [0.,0.]
real(dp), dimension(2) :: yrange = [0.,0.]

! ----------------------

call date_and_time(today,now)

ncstat = nf90_create(outfile,nf90_hdf5,ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','BIOME1 netCDF output file')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'timestamp',today//' '//now(1:4))
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Jobfile',jobfile)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','CF-1.11')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----

if (gridinfo%isproj) then

  ! projected grid: supply separate variables for geodetic longitude and latitude
  
  ! ----
  ! x coordinate
  
  ncstat = nf90_def_dim(ofid,'x',gridinfo%cntx,dimid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_def_var(ofid,'x',nf90_double,dimid,varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','easting')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','meters')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  dimids(1) = dimid
  
  ! ----
  ! y coordinate
  
  ncstat = nf90_def_dim(ofid,'y',gridinfo%cnty,dimid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_def_var(ofid,'y',nf90_double,dimid,varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','northing')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','meters')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  dimids(2) = dimid

  ! ----
  ! geodetic longitude

  ncstat = nf90_def_var(ofid,'lon',nf90_double,dimids(1:2),varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


  ! ----
  ! geodetic latitude

  ncstat = nf90_def_var(ofid,'lat',nf90_double,dimids(1:2),varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
else

  ! lon-lat grid

  ! ----
  ! x coordinate
  
  ncstat = nf90_def_dim(ofid,'lon',gridinfo%cntx,dimid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  dimids(1) = dimid
  
  ! ----
  ! y coordinate
  
  ncstat = nf90_def_dim(ofid,'lat',gridinfo%cnty,dimid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  dimids(2) = dimid

end if

! ----
! PFT coordinate

ncstat = nf90_def_dim(ofid,'PFT',npft,dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'PFT',nf90_int,dimid,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','plant functional type')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','type')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',[1,13])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

dimids(3) = dimid

! ----
! coordinate month

ncstat = nf90_def_dim(ofid,'month',nmos,dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'month',nf90_int,dimid,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','month')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','month')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',[1,12])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

dimids(4) = dimid

! ----

chunks(1) = min(gridinfo%cntx,200)
chunks(2) = min(gridinfo%cnty,200)
chunks(3) = 1
chunks(4) = 1

! ----
! monthly mean direct-beam surface radiation

ncstat = nf90_def_var(ofid,'rdirect',nf90_float,[dimids(1),dimids(2),dimids(4)],varid,chunksizes=[chunks(1),chunks(2),chunks(4)],deflate_level=1,shuffle=.false.)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','direct-beam surface shortwave radiation')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','MJ m-2 d-1')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! monthly mean diffuse surface radiation

ncstat = nf90_def_var(ofid,'rdiffuse',nf90_float,[dimids(1),dimids(2),dimids(4)],varid,chunksizes=[chunks(1),chunks(2),chunks(4)],deflate_level=1,shuffle=.false.)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','diffuse-beam surface shortwave radiation')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','MJ m-2 d-1')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! monthly PET

ncstat = nf90_def_var(ofid,'mpet',nf90_float,[dimids(1),dimids(2),dimids(4)],varid,chunksizes=[chunks(1),chunks(2),chunks(4)],deflate_level=1,shuffle=.false.)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','potential evapotranspiration')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','mm')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! monthly alpha

ncstat = nf90_def_var(ofid,'alpha',nf90_float,[dimids(1),dimids(2),dimids(4)],varid,chunksizes=[chunks(1),chunks(2),chunks(4)],deflate_level=1,shuffle=.false.)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','ratio of AET to PET')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','fraction')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',rmissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! biome (category)

ncstat = nf90_def_var(ofid,'biome',nf90_short,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.false.)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','biome')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','biome')
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',imissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',imissing)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----

ncstat = nf90_enddef(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! write the coordinate variables

if (gridinfo%isproj) then

  ncstat = nf90_inq_varid(ofid,'x',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_put_var(ofid,varid,coords(:,1)%xcoord)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ofid,'y',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,coords(1,:)%ycoord)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ofid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,coords%geolon)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ofid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,coords%geolat)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else

  ncstat = nf90_inq_varid(ofid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,coords(:,1)%geolon)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(ofid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,coords(1,:)%geolat)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end if

ncstat = nf90_inq_varid(ofid,'month',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,month)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'PFT',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,pft)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine genoutputfile

! -----------------------------------------------------

subroutine writereal3d(ofid,gridinfo,pixel,varname,ovar)

use parametersmod, only : sp,rmissing
use typesmod,      only : gridinfotype,pixeltype,metvars_monthly
use netcdf
use errormod,      only : ncstat,netcdf_err

implicit none

! arguments

integer,                       intent(in) :: ofid
type(gridinfotype),            intent(in) :: gridinfo
type(pixeltype), dimension(:), intent(in) :: pixel
character(*),                  intent(in) :: varname
real(sp), dimension(:,:),      intent(in) :: ovar

! local variables

integer :: cntx
integer :: cnty

integer :: n
integer :: i
integer :: x
integer :: y

integer :: varlen
integer :: varid

real(sp), allocatable, dimension(:,:,:) :: mvar

! ----

varlen = size(ovar,dim=2)

cntx = gridinfo%cntx
cnty = gridinfo%cnty

allocate(mvar(cntx,cnty,varlen))

mvar = rmissing

n = size(pixel)

do i = 1,n

  x = pixel(i)%x
  y = pixel(i)%y
  
  mvar(x,y,:) = ovar(i,:)

end do

ncstat = nf90_inq_varid(ofid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,mvar)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine writereal3d

! -----------------------------------------------------

subroutine writereal2d(ofid,gridinfo,pixel,varname,ovar)

use parametersmod, only : sp,i2,imissing
use typesmod,      only : gridinfotype,pixeltype,metvars_monthly
use netcdf
use errormod,      only : ncstat,netcdf_err

implicit none

! arguments

integer,                       intent(in) :: ofid
type(gridinfotype),            intent(in) :: gridinfo
type(pixeltype), dimension(:), intent(in) :: pixel
character(*),                  intent(in) :: varname
real(sp), dimension(:),        intent(in) :: ovar

! local variables

integer :: cntx
integer :: cnty

integer :: n
integer :: i
integer :: x
integer :: y

integer :: varid

integer(sp), allocatable, dimension(:,:) :: mvar

! ----

cntx = gridinfo%cntx
cnty = gridinfo%cnty

allocate(mvar(cntx,cnty))

mvar = imissing

n = size(pixel)

do i = 1,n

  x = pixel(i)%x
  y = pixel(i)%y
  
  mvar(x,y) = ovar(i)

end do

ncstat = nf90_inq_varid(ofid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,mvar)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine writereal2d

! -----------------------------------------------------

subroutine writeinteger2d(ofid,gridinfo,pixel,varname,ovar)

use parametersmod, only : sp,i2,imissing
use typesmod,      only : gridinfotype,pixeltype,metvars_monthly
use netcdf
use errormod,      only : ncstat,netcdf_err

implicit none

! arguments

integer,                       intent(in) :: ofid
type(gridinfotype),            intent(in) :: gridinfo
type(pixeltype), dimension(:), intent(in) :: pixel
character(*),                  intent(in) :: varname
integer(i2), dimension(:),      intent(in) :: ovar

! local variables

integer :: cntx
integer :: cnty

integer :: n
integer :: i
integer :: x
integer :: y

integer :: varid

integer(sp), allocatable, dimension(:,:) :: mvar

! ----

cntx = gridinfo%cntx
cnty = gridinfo%cnty

allocate(mvar(cntx,cnty))

mvar = imissing

n = size(pixel)

do i = 1,n

  x = pixel(i)%x
  y = pixel(i)%y
  
  mvar(x,y) = ovar(i)

end do

ncstat = nf90_inq_varid(ofid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,mvar)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine writeinteger2d

! -----------------------------------------------------

subroutine closeoutput(ofid)

use netcdf
use errormod, only : ncstat,netcdf_err

implicit none

integer, intent(in) :: ofid

! ---

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine closeoutput

! -----------------------------------------------------

end module netcdfoutputmod