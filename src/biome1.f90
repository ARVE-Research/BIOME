program biome1

use parametersmod
use coordsmod
use readdatamod
use typesmod  ! going to use all of the types, but should specify
use netcdfoutputmod
use newsplinemod
use calendarmod, only : nd_365,present_mon_noleap,nm
use randomdistmod, only : ran_seed
use weathergenmod, only : weathergen

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


type(pixeltype),   allocatable, dimension(:) :: pixel
type(metvars_in),  allocatable, dimension(:) :: met_in
type(metvars_out), allocatable, dimension(:) :: met_out

integer :: ncells

integer :: x
integer :: y
integer :: m
integer :: i
integer :: j

integer :: cntx
integer :: cnty

integer :: ofid

integer, dimension(nm) :: ndm

! logical, allocatable, dimension(:,:) :: valid

integer :: d
integer :: nd

integer :: memreq
integer :: maxmem = 20000  ! default amount of memory allowed for input data


namelist /joboptions/ gridinfo,terrainfile,climatefile,soilfile,maxmem

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

ncells = count(soil(:,:,1)%whc > 0. .and. climate(:,:,1)%pre >= 0.)

write(0,'(a,i0,a)')' there are ',ncells,' valid pixels'

allocate(pixel(ncells))

i = 1

do y = 1,cnty
  do x = 1,cntx

    if (soil(x,y,1)%whc < 0. .or. climate(x,y,1)%pre < 0.) cycle
    
    pixel(i)%x = x
    pixel(i)%y = y
    
    i = i + 1

  end do
end do

nd  = nd_365
ndm = int(present_mon_noleap)

! at this point should make a quick check for total memory requirement and take appropriate action if the amount is too large

allocate(daily(1,nd))

memreq = ncells * sizeof(daily) / 1048576

write(0,'(i0,a,i0,a)')ncells,' valid cells takes ',memreq,' MB memory'

deallocate(daily)

if (memreq > maxmem) then
  write(0,'(a,i0,a)')' ERROR daily data memory requirement of ',memreq,' MB exceeds user set maximum, aborting.'
  stop
end if

allocate(daily(ncells,nd))

! interpolate selected monthly meteorological variables to means-preserving smooth daily estimates

do i = 1,ncells

  x = pixel(i)%x
  y = pixel(i)%y
  
  call newspline(climate(x,y,:)%tmp,ndm,[climate(x,y,12)%tmp,climate(x,y,1)%tmp],daily(i,:)%tmp)
  call newspline(climate(x,y,:)%dtr,ndm,[climate(x,y,12)%dtr,climate(x,y,1)%dtr],daily(i,:)%dtr)
  call newspline(climate(x,y,:)%cld,ndm,[climate(x,y,12)%cld,climate(x,y,1)%cld],daily(i,:)%cld,llim=0.,ulim=100.)
  call newspline(climate(x,y,:)%wnd,ndm,[climate(x,y,12)%wnd,climate(x,y,1)%wnd],daily(i,:)%wnd,llim=0.)
    
end do

allocate(met_in(ncells))
allocate(met_out(ncells))

! initalize the random number generator

do i = 1,ncells
  call ran_seed(-10,met_in(i)%rndst)
end do

! initialize prior precipitation and the weather residuals

met_in%pday(1) = .false.
met_in%pday(2) = .false.

do i = 1,4
  met_in%resid(i) = 0.
end do

! start daily loop

write(0,*)'start daily loop'

j = 1

do m = 1,12
  do d = 1,ndm(m)
  
    ! loop over valid cells

    do i = 1,ncells
    
      x = pixel(i)%x
      y = pixel(i)%y
    
      met_in(i)%prec = climate(x,y,m)%pre
      met_in(i)%wetf = climate(x,y,m)%wet
      met_in(i)%wetd = climate(x,y,m)%wet * real(present_mon_noleap(m))
      
      met_in(i)%tmin = daily(i,j)%tmp - 0.5 * daily(i,j)%dtr
      met_in(i)%tmax = daily(i,j)%tmp + 0.5 * daily(i,j)%dtr
      met_in(i)%cldf = daily(i,j)%cld
      met_in(i)%wind = daily(i,j)%wnd

      ! generate daily meteorology for all valid cells for one day
      
      write(0,*)m,d,j,i,'weathergen'
  
      call weathergen(met_in(i),met_out(i))

    end do ! cells
    
    j = j + 1
    
  end do   ! days in the month
end do     ! month

! ---------------------------------
! write model output

end program biome1