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
type(climatetype), allocatable, dimension(:)     :: daily

integer :: x
integer :: y

integer :: cntx
integer :: cnty

integer :: ofid

integer, dimension(nm) :: ndm

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

write(0,*)gridinfo

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
! allocate daily vectors for meteorological variables

allocate(daily(nd_365))

ndm = int(present_mon_noleap)

! grid loop

do y = 1,cnty
  do x = 1,cntx
  
    if (soil(x,y,1)%whc < 0.) cycle

    ! interpolate climate to pseudo-daily using means-preserving algorithm
  
    ! call newspline(climate(x,y,:)%tmp,ndm,daily%tmp,prec=1)
  
  
    ! use weather generator to generate daily weather
    
    ! ---------------------------------
    ! calculate daily surface insolation and potential evapotranspiration (PET)

  end do
end do

! end grid loop

! ---------------------------------
! write model output

end program biome1