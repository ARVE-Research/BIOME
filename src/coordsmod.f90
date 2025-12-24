module coordsmod

use parametersmod, only : dp

implicit none

public :: parsecoords
public :: calcpixels

contains

! -----------------------------------------------------

subroutine parsecoords(coordstring,gridinfo)

! subroutine to parse a coordinate string in characters separated by "/" into 4 separate real values

use typesmod, only : gridinfotype

implicit none

! arguments

character(*),       intent(in)  :: coordstring
type(gridinfotype), intent(out) :: gridinfo

! local variables

real(dp), dimension(4) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

! -----

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

! if only two values are provided (single point), assibn 

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

gridinfo%xmin = val(1)
gridinfo%xmax = val(2)
gridinfo%ymin = val(3)
gridinfo%ymax = val(4)

end subroutine parsecoords

! -----------------------------------------------------

subroutine calcpixels(infile,gridinfo)

! given infile and coordinate boundaries, calculate srtx, cntx, srty, and cnty to extract data from a window

use netcdf
use errormod, only   : ncstat,netcdf_err
use utilitymod, only : pos
use typesmod,   only : gridinfotype

implicit none

! arguments

character(*),       intent(in)    :: infile
type(gridinfotype), intent(inout) :: gridinfo

! local variables

integer :: ncid
integer :: dimid
integer :: varid

integer :: xlen
integer :: ylen

integer,  dimension(2) :: xpos
integer,  dimension(2) :: ypos

real(dp), allocatable, dimension(:) :: xcoords
real(dp), allocatable, dimension(:) :: ycoords

character(3) :: xcoordname
character(3) :: ycoordname

! -----------------------------------------------------

if (gridinfo%isproj) then
  xcoordname = 'x'
  ycoordname = 'y'
else
  xcoordname = 'lon'
  ycoordname = 'lat'
end if

! open input file

ncstat = nf90_open(infile,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ---
! retrieve x and y coordinate vectors
  
ncstat = nf90_inq_dimid(ncid,xcoordname,dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(xcoords(xlen))

ncstat = nf90_inq_varid(ncid,xcoordname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,xcoords)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ---

ncstat = nf90_inq_dimid(ncid,ycoordname,dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(ycoords(ylen))

ncstat = nf90_inq_varid(ncid,ycoordname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ycoords)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ---

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! -----------------------------------------------------
! find nearest pixel to minimum coordinates (single point or window)
 
xpos = [pos(xcoords,gridinfo%xmin),pos(xcoords,gridinfo%xmax)]

ypos = [pos(ycoords,gridinfo%ymin),pos(ycoords,gridinfo%ymax)]

! -----------------------------------------------------
! calculate the indices of the starting point for the desired pixel, and the number of pixels to be retrieved 

gridinfo%srtx = xpos(1)
gridinfo%cntx = max(xpos(2) - xpos(1),1)

gridinfo%srty = ypos(1)
gridinfo%cnty = max(ypos(2) - ypos(1),1)

! write(0,*)'input data xlen: ',xlen
! write(0,*)'input data ylen: ',ylen
! write(0,*)gridinfo%srtx,gridinfo%srty
! write(0,*)gridinfo%cntx,gridinfo%cnty

end subroutine calcpixels

! -----------------------------------------------------

end module coordsmod
