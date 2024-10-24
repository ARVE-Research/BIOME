program biome1

use parametersmod
use coordsmod
use readdatamod
use typesmod
use netcdfoutputmod
use newsplinemod
use calendarmod, only : nd_365,present_mon_noleap,nm

implicit none

character(200) :: jobfile
character(200) :: terrainfile
character(200) :: climatefile
character(200) :: soilfile
character(200) :: outfile

character(60) :: coordstring

type(gridinfotype) :: gridinfo

type(coordstype),  allocatable, dimension(:,:)   :: coords
type(terraintype), allocatable, dimension(:,:)   :: terrain
type(climatetype), allocatable, dimension(:,:,:) :: climate
type(soiltype),    allocatable, dimension(:,:,:) :: soil
type(climatetype), allocatable, dimension(:,:)   :: daily

integer :: x
integer :: y
integer :: m
integer :: i

integer :: cntx
integer :: cnty

integer :: ofid

integer, dimension(nm) :: ndm

integer :: valid
integer :: nd

namelist /joboptions/ gridinfo,terrainfile,climatefile,soilfile

! ---------------------------------
! read the job file with input data

call getarg(1,jobfile)

open(10,file=jobfile)

read(10,nml=joboptions)

! ---------------------------------
! get the coordinates for the run

call getarg(2,coordstring)

call parsecoords(coordstring,gridinfo)

call calcpixels(climatefile,gridinfo)

! write(0,*)gridinfo

cntx = gridinfo%cntx
cnty = gridinfo%cnty

allocate(coords(cntx,cnty))
allocate(terrain(cntx,cnty))
allocate(climate(cntx,cnty,12))
allocate(soil(cntx,cnty,6))

! ---------------------------------
! read the input data (coordinate variables, terrain, and land fraction)

call readcoords(climatefile,gridinfo,coords)

call readterrain(climatefile,gridinfo,terrain)

! read the input data (monthly climate and soil)

call readclimate(climatefile,gridinfo,climate)

call readsoil(soilfile,gridinfo,soil)

! ---------------------------------
! generate the output file

call getarg(3,outfile)

call genoutputfile(jobfile,outfile,gridinfo,ofid)

! ---------------------------------
! check for valid pixels

valid = count(soil(:,:,1)%whc > 0.)

write(0,'(a,i0,a)')' there are ',valid,' valid pixels'

nd  = nd_365
ndm = int(present_mon_noleap)

! at this point should make a quick check for total memory requirement and take appropriate action if the amount is too large

allocate(daily(valid,nd))

write(0,'(a,i0,a)')' daily smoothed metvars size: ',sizeof(daily) / 1048576,' MB'

i = 1

do y = 1,cnty
  do x = 1,cntx
  
    if (soil(x,y,1)%whc <= 0.) cycle
    
    call newspline(climate(x,y,:)%tmp,ndm,[climate(x,y,12)%tmp,climate(x,y,1)%tmp],daily(i,:)%tmp)
    call newspline(climate(x,y,:)%dtr,ndm,[climate(x,y,12)%dtr,climate(x,y,1)%dtr],daily(i,:)%dtr)
    call newspline(climate(x,y,:)%cld,ndm,[climate(x,y,12)%cld,climate(x,y,1)%cld],daily(i,:)%cld,llim=0.,ulim=100.)
    call newspline(climate(x,y,:)%wnd,ndm,[climate(x,y,12)%wnd,climate(x,y,1)%wnd],daily(i,:)%wnd,llim=0.)
    
    i = i + 1

  end do
end do

! daily loop

  ! grid loop


! ---------------------------------
! write model output

end program biome1