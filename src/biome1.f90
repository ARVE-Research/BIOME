program biome1

use parametersmod
use coordsmod
use readdatamod
use typesmod         ! going to use all of the types, but should specify
use netcdfoutputmod, only : genoutputfile,writereal3d,closeoutput
use utilitymod,      only : bp2ce,leapyear,overprint
use newsplinemod
use orbitmod,        only : getorbitpars
use calendarmod,     only : initcalendar
use randomdistmod,   only : ran_seed
use weathergenmod,   only : weathergen
use calendarmod
use insolationmod,   only : truelon,insol
use diurnaltempmod,  only : diurnaltemp
use radiationmod,    only : initairmass,elev_corr,radpet,Pjj
use soilwatermod,    only : calcwhc,soilwater

implicit none

character(200) :: jobfile
character(200) :: terrainfile
character(200) :: climatefile
character(200) :: soilfile
character(200) :: outfile

character(60) :: coordstring

type(gridinfotype) :: gridinfo
type(orbitpars)    :: orbit
type(solarpars)    :: solar
type(calendartype) :: cal

type(soilcoordstype),  allocatable, dimension(:)     :: soilcoords
type(pixeltype),       allocatable, dimension(:)     :: pixel
type(metvars_in),      allocatable, dimension(:)     :: met_in
type(coordstype),      allocatable, dimension(:,:)   :: coords
type(terraintype),     allocatable, dimension(:,:)   :: terrain
type(climatetype),     allocatable, dimension(:,:)   :: daily
type(soilwatertype),   allocatable, dimension(:)     :: soilw
type(metvars_daily),   allocatable, dimension(:,:)   :: dmet
type(metvars_monthly), allocatable, dimension(:,:)   :: mmet
type(climatetype),     allocatable, dimension(:,:,:) :: climate
type(soiltype),        allocatable, dimension(:,:,:) :: soil


integer :: ncells

integer :: x
integer :: y
integer :: m
integer :: i
integer :: j

integer :: doy

integer :: cntx
integer :: cnty

integer :: ofid

integer :: ndy
integer, dimension(nmos) :: ndm

integer :: d

integer(i8) :: memreq
integer(i8) :: maxmem = 20000_i8  ! default amount of memory allowed for input data

integer :: yrbp

real(dp) :: slon  ! true solar longitude, for insolation calculations (radians)
real(dp) :: latr  ! geodesic latitude (radians)
real(dp) :: latd  ! geodesic latitude (degrees)

real(sp) :: albedo

character(40) :: status_line

namelist /joboptions/ gridinfo,terrainfile,climatefile,soilfile,maxmem

! ---------------------------------
! read the job file with input data

call getarg(1,jobfile)

open(10,file=jobfile)

read(10,nml=joboptions)

! ---------------------------------
! calculate orbital parameters and month lengths for the simulation year

yrbp = -73

call initcalendar(yrbp,orbit,cal)

ndy = cal%ndyr
ndm = cal%ndmi

write(0,*)'orbit: ',orbit
write(0,*)'calendar: ',ndy,ndm

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
allocate(climate(cntx,cnty,nmos))
allocate(soilcoords(6))
allocate(soil(cntx,cnty,6))

! ---------------------------------
! read the input data (coordinate variables, terrain, and land fraction)

call readcoords(climatefile,gridinfo,coords)

call readterrain(climatefile,gridinfo,terrain)

! read the input data (monthly climate and soil)

call readclimate(climatefile,gridinfo,climate)

call readsoil(soilfile,gridinfo,terrain,soilcoords,soil)

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

allocate(mmet(ncells,12))  ! monthly output

i = 1

do y = 1,cnty
  do x = 1,cntx

    if (soil(x,y,1)%whc == rmissing .or. climate(x,y,1)%pre == rmissing) cycle
    
    pixel(i)%x = x
    pixel(i)%y = y
    
    pixel(i)%lon = coords(x,y)%geolon
    pixel(i)%lat = coords(x,y)%geolat
    
    ! temperature of the coldest month

    pixel(i)%tcm = minval(climate(x,y,:)%tmp)
    
    ! precipitation equitability index, used in airmass calculations for surface radiation
    
    pixel(i)%Pann = sum(climate(x,y,:)%pre)

    pixel(i)%Pjj = Pjj(climate(x,y,:)%tmp,climate(x,y,:)%pre) 
    
    ! relative atmospheric pressure
    
    pixel(i)%Ratm = elev_corr(terrain(x,y)%elv)
    
    i = i + 1

  end do
end do

! at this point should make a quick check for total memory requirement and take appropriate action if the amount is too large

allocate(daily(1,ndy))

memreq = ncells * sizeof(daily) / 1048576_i8

write(0,'(a,i0,a,i0,a)')' ',ncells,' valid gridcells takes ',memreq,' MB memory'

deallocate(daily)

if (memreq > maxmem) then
  write(0,'(a,i0,a)')' ERROR daily data memory requirement of ',memreq,' MB exceeds user set maximum, aborting.'
  stop
end if

allocate(daily(ncells,ndy))

allocate(soilw(ncells))

! interpolate selected monthly meteorological variables to means-preserving smooth daily estimates

write(0,*)'calculate smoothed meteorology'

do i = 1,ncells

  x = pixel(i)%x
  y = pixel(i)%y
  
  call newspline(climate(x,y,:)%tmp,ndm,[climate(x,y,nmos)%tmp,climate(x,y,1)%tmp],daily(i,:)%tmp)
  call newspline(climate(x,y,:)%dtr,ndm,[climate(x,y,nmos)%dtr,climate(x,y,1)%dtr],daily(i,:)%dtr)
  call newspline(climate(x,y,:)%cld,ndm,[climate(x,y,nmos)%cld,climate(x,y,1)%cld],daily(i,:)%cld,llim=0.,ulim=100.)
  call newspline(climate(x,y,:)%wnd,ndm,[climate(x,y,nmos)%wnd,climate(x,y,1)%wnd],daily(i,:)%wnd,llim=0.)
  
  ! calculate the column-integrated soil water holding capacity
  
  call calcwhc(terrain(x,y),soilcoords,soil(x,y,:),soilw(i))

end do

allocate(met_in(ncells))
allocate(dmet(ncells,2))  ! for the current and next day

! initalize the random number generator

write(0,*)'seed random number generator'

do i = 1,ncells
  call ran_seed(-104576,met_in(i)%rndst)
end do

! initialize airmass parameters

call initairmass()

! initialize prior precipitation and the weather residuals

met_in%pday(1) = .false.
met_in%pday(2) = .false.

do i = 1,4
  met_in%resid(i) = 0.
end do

dmet%tmin = 0.
dmet%tmax = 0.
dmet%prec = 0.
dmet%cldf = 0.
dmet%wind = 0.
dmet%rad0 = 0.
dmet%dayl = 0.
dmet%tday = 0.
dmet%tnight= 0.

soilw%w = soilw%whc

! write(0,*)'initial soilw%w',soilw%w

! ---------------------------------
! initialization daily loop

! 10 format(2i4,f8.3,f6.1,2f7.1,f7.3,f7.1,3f7.1,2f7.1,f7.1,2f10.1,f7.1)

write(0,*)'start initialization daily loop'

doy = 1  ! day of year counter

do m = 1,nmos
  do d = 1,ndm(m)
  
    write(status_line,'(a,i0,a,i0)')' working on ',m,' ',d
    call overprint(status_line)
    
    ! calculate true solar longitude (for daylength)
    ! only need to do this once per day, not for each cell
    
    slon = truelon(orbit,cal,doy)
  
    ! loop over valid cells

    do i = 1,ncells

      x = pixel(i)%x
      y = pixel(i)%y
      
      latd = coords(x,y)%geolat
      latr = latd * pir
      
      ! variable initializations

      met_in(i)%prec = climate(x,y,m)%pre
      met_in(i)%wetf = climate(x,y,m)%wet
      met_in(i)%wetd = climate(x,y,m)%wet * cal%ndmr(m)

      met_in(i)%tmin = daily(i,doy)%tmp - 0.5 * daily(i,doy)%dtr
      met_in(i)%tmax = daily(i,doy)%tmp + 0.5 * daily(i,doy)%dtr
      met_in(i)%cldf = daily(i,doy)%cld * 0.01
      met_in(i)%wind = daily(i,doy)%wnd

      ! generate daily meteorology for all valid cells for one day

      call weathergen(met_in(i),dmet(i,2))
      
      ! write(0,*)'WG:',dmet(i,2)%tmin,dmet(i,2)%tmax

      ! calculate top-of-the-atmosphere insolation and daylength
      
      call insol(slon,orbit,latr,solar)
      
      dmet(i,2)%rad0 = solar%rad0
      dmet(i,2)%dayl = solar%dayl

      ! dmet(:,1) = current day values, dmet(:,2) = next day day values
      
      dmet(i,:) = eoshift(dmet(i,:),-1,dmet(i,2))
      
      ! calculate integrated day- and night-time temperature

      call diurnaltemp(dmet(i,:))
      
      ! estimate dewpoint temperature 
      
      albedo = 0.17  ! shortwave albedo for vegetated surfaces. should vary in space and time; placeholder for now
      
      ! surface radiation budget and potential evapotranspiration
  
      call radpet(pixel(i),solar,albedo,dmet(i,1))
      
      ! soil water balance, including actual evapotranspiration and alpha
      
      call soilwater(dmet(i,1),soilw(i))

    end do ! cells

    doy = doy + 1
    
  end do   ! days in the month
end do     ! month

! ---------------------------------
! computation daily loop

write(0,*)
write(0,*)'start computation daily loop'

mmet%mpet = 0.

doy = 1

do m = 1,nmos
  do d = 1,ndm(m)

    write(status_line,'(a,i0,a,i0)')' working on ',m,' ',d
    call overprint(status_line)

    ! calculate true solar longitude for daylength, 
    ! only need to do this once per day, not for each cell
    
    slon = truelon(orbit,cal,doy)
  
    ! loop over valid cells

    do i = 1,ncells
    
      x = pixel(i)%x
      y = pixel(i)%y
    
      met_in(i)%prec = climate(x,y,m)%pre
      met_in(i)%wetf = climate(x,y,m)%wet
      met_in(i)%wetd = climate(x,y,m)%wet * real(present_mon_noleap(m))
      
      met_in(i)%tmin = daily(i,doy)%tmp - 0.5 * daily(i,doy)%dtr
      met_in(i)%tmax = daily(i,doy)%tmp + 0.5 * daily(i,doy)%dtr
      met_in(i)%cldf = daily(i,doy)%cld * 0.01
      met_in(i)%wind = daily(i,doy)%wnd

      ! generate daily meteorology for all valid cells for one day
  
      call weathergen(met_in(i),dmet(i,1))
      
      ! calculate daylength and insolation
      
      latr = coords(x,y)%geolat * pir
      
      call insol(slon,orbit,latr,solar)
      
      dmet(i,2)%rad0 = solar%rad0
      dmet(i,2)%dayl = solar%dayl
      
      ! dmet(:,1) = current day values, dmet(:,2) = next day day values
      
      dmet(i,:) = eoshift(dmet(i,:),-1,dmet(i,2))
      
      ! calculate integrated day- and night-time temperature

      call diurnaltemp(dmet(i,:))
      
      albedo = 0.17  ! shortwave albedo for vegetated surfaces. should vary in space and time; placeholder for now

      ! calculate surface solar radiation and potential evapotranspiration

      call radpet(pixel(i),solar,albedo,dmet(i,1))
     
      ! soil water balance, including actual evapotranspiration and alpha

      call soilwater(dmet(i,1),soilw(i))
      
      ! monthly summaries

      mmet(i,m)%mpet  = mmet(i,m)%mpet  + dmet(i,1)%dpet
      mmet(i,m)%alpha = mmet(i,m)%alpha + dmet(i,1)%alpha

    end do ! cells
    
    doy = doy + 1
    
  end do   ! days in the month
end do     ! month

stop
write(0,*)

! ---------------------------------
! write model output

call writereal3d(ofid,gridinfo,pixel,'mpet',mmet%mpet)
call writereal3d(ofid,gridinfo,pixel,'alpha',mmet%alpha)

! ---------------------------------

call closeoutput(ofid)

! ---------------------------------
! finish

end program biome1
