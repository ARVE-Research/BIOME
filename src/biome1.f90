program biome1

use parametersmod
use coordsmod
use readdatamod
use typesmod  ! going to use all of the types, but should specify
use netcdfoutputmod
use newsplinemod
use orbitmod, only : getorbitpars
use calendarmod, only : initcalendar
use randomdistmod, only : ran_seed
use weathergenmod, only : weathergen
use calendarmod

implicit none

character(200) :: jobfile
character(200) :: terrainfile
character(200) :: climatefile
character(200) :: soilfile
character(200) :: outfile

character(60) :: coordstring

type(gridinfotype) :: gridinfo
type(orbitpars)    :: orbit

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

integer, dimension(nmos) :: ndm

type(calendartype) :: noleap
type(calendartype) :: leapyr

! logical, allocatable, dimension(:,:) :: valid

integer :: d
integer :: nd

integer(i8) :: memreq
integer(i8) :: maxmem = 20000_i8  ! default amount of memory allowed for input data

integer :: yrbp

integer,  dimension(nmos) :: imonlen
real(dp), dimension(nmos) :: rmonlen
real(dp), dimension(nmos) :: rmonbeg
real(dp), dimension(nmos) :: rmonmid
real(dp), dimension(nmos) :: rmonend

namelist /joboptions/ gridinfo,terrainfile,climatefile,soilfile,maxmem

! ---------------------------------
! read the job file with input data

call getarg(1,jobfile)

open(10,file=jobfile)

read(10,nml=joboptions)

! ---------------------------------
! calculate orbital parameters and month lengths for the simulation year

yrbp = 6000

call initcalendar(yrbp,orbit,noleap,leapyr)

stop

! ---------------------------------
! get the coordinates for the run

call getarg(2,coordstring)

call parsecoords(coordstring,gridinfo)

call calcpixels(climatefile,gridinfo)

cntx = gridinfo%cntx
cnty = gridinfo%cnty

write(0,*)'allocate rectangular arrays',cntx,cnty

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

write(0,*)'go genoutput'

call genoutputfile(jobfile,outfile,gridinfo,ofid)

! ---------------------------------
! check for valid pixels

ncells = count(soil(:,:,1)%whc /= rmissing .and. climate(:,:,1)%tmp /= rmissing)

write(0,'(a,i0,a)')' there are ',ncells,' valid gridcells'

allocate(pixel(ncells))

i = 1

do y = 1,cnty
  do x = 1,cntx

    if (soil(x,y,1)%whc == rmissing .or. climate(x,y,1)%pre == rmissing) cycle
    
    pixel(i)%x = x
    pixel(i)%y = y
    
    i = i + 1

  end do
end do

nd  = nd_365
ndm = int(present_mon_noleap)

! at this point should make a quick check for total memory requirement and take appropriate action if the amount is too large

allocate(daily(1,nd))

memreq = ncells * sizeof(daily) / 1048576_i8

write(0,'(a,i0,a,i0,a)')' ',ncells,' valid gridcells takes ',memreq,' MB memory'

deallocate(daily)

if (memreq > maxmem) then
  write(0,'(a,i0,a)')' ERROR daily data memory requirement of ',memreq,' MB exceeds user set maximum, aborting.'
  stop
end if

allocate(daily(ncells,nd))

! interpolate selected monthly meteorological variables to means-preserving smooth daily estimates

write(0,*)'go newspline'

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

write(0,*)'go random seed'

do i = 1,ncells
  call ran_seed(-104576,met_in(i)%rndst)
end do

! initialize prior precipitation and the weather residuals

met_in%pday(1) = .false.
met_in%pday(2) = .false.

do i = 1,4
  met_in%resid(i) = 0.
end do

! ---------------------------------
! initialization daily loop

write(0,*)'start initialization daily loop'

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
      met_in(i)%cldf = daily(i,j)%cld * 0.01
      met_in(i)%wind = daily(i,j)%wnd

      ! generate daily meteorology for all valid cells for one day
  
      call weathergen(met_in(i),met_out(i))

      write(*,*)'2022',m,d,j,met_out(i)%prec,met_out(i)%tmin,met_out(i)%tmax,met_out(i)%cldf,met_out(i)%wind

    end do ! cells
    
    j = j + 1
    
  end do   ! days in the month
end do     ! month

! ---------------------------------
! computation daily loop

write(0,*)'start computation daily loop'

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
      met_in(i)%cldf = daily(i,j)%cld * 0.01
      met_in(i)%wind = daily(i,j)%wnd

      ! generate daily meteorology for all valid cells for one day
  
      call weathergen(met_in(i),met_out(i))

      write(*,*)'2023',m,d,j,met_out(i)%prec,met_out(i)%tmin,met_out(i)%tmax,met_out(i)%cldf,met_out(i)%wind
      
      ! calculate daylength and insolation
      
      ! estimate integrated daytime and nighttime temperatures
      
      ! estimate humidity
      
      ! estimate potential evapotranspiration for day and nighttime separately

    end do ! cells
    
    j = j + 1
    
  end do   ! days in the month
end do     ! month

! ---------------------------------
! write model output

end program biome1