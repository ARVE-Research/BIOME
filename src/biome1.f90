program biome1

use parametersmod
use coordsmod
use readdatamod
use typesmod         ! going to use all of the types, but should specify
use netcdfoutputmod
use utilitymod,      only : bp2ce,leapyear,overprint,imaxloc
use newsplinemod
use orbitmod,        only : getorbitpars
use calendarmod,     only : initcalendar
use randomdistmod,   only : ran_seed
use weathergenmod,   only : weathergen
use calendarmod
use insolationmod,   only : truelon,insol
use diurnaltempmod,  only : diurnaltemp
use physicsmod,      only : stdP
use airmassmod,      only : initairmass,elev_corr,Pjj
use radiationmod,    only : radpet
use soilwatermod,    only : calcwhc,soilwater
use snowmod,         only : Tt,snow
use calcbiomemod,    only : calcbiome

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
type(dayclimatetype),  allocatable, dimension(:,:)   :: daily
type(soilwatertype),   allocatable, dimension(:)     :: soilw
type(metvars_daily),   allocatable, dimension(:)     :: dmet0  ! current day meteorology
type(metvars_daily),   allocatable, dimension(:)     :: dmet1  ! next day meteorology
type(metvars_monthly), allocatable, dimension(:,:)   :: mmet
type(monclimatetype),  allocatable, dimension(:,:,:) :: climate
type(soiltype),        allocatable, dimension(:,:,:) :: soil

real(sp), allocatable, dimension(:) :: tmin
real(sp), allocatable, dimension(:) :: tmax

integer :: ncells

integer :: x
integer :: y
integer :: m
integer :: i
integer :: j

integer :: doy
integer :: d1

integer :: cntx
integer :: cnty

integer :: ofid

integer :: ndy
integer, dimension(nmos) :: ndm

integer :: wm
integer :: d

integer(i8) :: memreq
integer(i8) :: maxmem = 20000_i8  ! default amount of memory allowed for input data

integer :: yrbp

real(dp) :: slon  ! true solar longitude, for insolation calculations (radians)
real(dp) :: latr  ! geodesic latitude (radians)
real(dp) :: latd  ! geodesic latitude (degrees)

real(sp) :: albedo

real(sp) :: dpd_day
real(sp) :: dpd_night

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

call genoutputfile(jobfile,outfile,gridinfo,coords,ofid)

! ---------------------------------
! check for valid pixels

ncells = count(soil(:,:,1)%whc /= rmissing .and. climate(:,:,1)%tmp /= rmissing .and. terrain(:,:)%elv /= rmissing)

write(0,'(a,i0,a)')' there are ',ncells,' valid gridcells'

allocate(pixel(ncells))

allocate(mmet(ncells,12))  ! monthly output

i = 1

do y = 1,cnty
  do x = 1,cntx

    if (soil(x,y,1)%whc == rmissing .or. climate(x,y,1)%pre == rmissing .or. terrain(x,y)%elv == rmissing) cycle
    
    pixel(i)%x = x
    pixel(i)%y = y
    
    pixel(i)%lon = coords(x,y)%geolon
    pixel(i)%lat = coords(x,y)%geolat

    pixel(i)%elv = terrain(x,y)%elv
    
    pixel(i)%slope  = 0. ! terrain(x,y)%slope
    pixel(i)%aspect = 0. ! terrain(x,y)%aspect
    
    ! mean temperature of the coldest month

    pixel(i)%tcm = minval(climate(x,y,:)%tmp)
    pixel(i)%twm = maxval(climate(x,y,:)%tmp)
    pixel(i)%wm  = imaxloc(climate(x,y,:)%tmp)

    ! precipitation equitability index, used in airmass calculations for surface radiation
    
    pixel(i)%Pann = sum(climate(x,y,:)%pre)

    pixel(i)%Pjj = Pjj(climate(x,y,:)%tmp,climate(x,y,:)%pre) 
    
    ! relative atmospheric pressure
    
    pixel(i)%Ratm = elev_corr(pixel(i)%elv)
    pixel(i)%P    = stdP(pixel(i)%elv)
    
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

allocate(tmin(nmos))
allocate(tmax(nmos))

allocate(soilw(ncells))

! ----------------------------------------------------------------------------------------------------------------
! interpolate selected monthly meteorological variables to means-preserving smooth daily estimates

write(0,*)'calculate smoothed meteorology'

do i = 1,ncells

  x = pixel(i)%x
  y = pixel(i)%y
  
  tmin = climate(x,y,:)%tmp - 0.5 * climate(x,y,:)%dtr
  tmax = climate(x,y,:)%tmp + 0.5 * climate(x,y,:)%dtr
  
  call newspline(tmin,ndm,[tmin(nmos),tmin(1)],daily(i,:)%tmin)
  call newspline(tmax,ndm,[tmax(nmos),tmax(1)],daily(i,:)%tmax)
    
  call newspline(climate(x,y,:)%cld,ndm,[climate(x,y,nmos)%cld,climate(x,y,1)%cld],daily(i,:)%cld,llim=0.,ulim=100.)
  call newspline(climate(x,y,:)%wnd,ndm,[climate(x,y,nmos)%wnd,climate(x,y,1)%wnd],daily(i,:)%wnd,llim=0.)
  
  ! calculate the column-integrated soil water holding capacity
  
  call calcwhc(terrain(x,y),soilcoords,soil(x,y,:),soilw(i))
  
  ! calculate the snow probability temperature Tt
  
  pixel(i)%Tt = Tt(pixel(i)%elv,pixel(i)%lat)
  
end do

deallocate(tmin)
deallocate(tmax)

allocate(met_in(ncells))
allocate(dmet0(ncells))    ! for the current day
allocate(dmet1(ncells))   ! for the next day

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

dmet0%tmin = 0.
dmet0%tmax = 0.
dmet0%prec = 0.
dmet0%cldf = 0.
dmet0%wind = 0.
dmet0%rad0 = 0.
dmet0%dayl = 0.
dmet0%tday = 0.
dmet0%tnight = 0.
dmet0%swe = 0.
dmet0%asnow = 0

soilw%w = soilw%whc

! write(0,*)'initial soilw%w',soilw%w

! ----------------------------------------------------------------------------------------------------------------
! initialize first day meteorology
! --- notes ---
! model operation starts at sunrise on the first day of the year
! daytime temperature reflects this day
! nighttime temperature is temperature going into the following day
! on the first day of operation, need to get the current and following day meteorology
! after first day, move next day into current
! -------------

write(0,*)'start initialization day one'

doy = 1  ! day of year counter

m = 1
d = 1
    
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

  met_in(i)%tmin = daily(i,doy)%tmin
  met_in(i)%tmax = daily(i,doy)%tmax
  met_in(i)%cldf = daily(i,doy)%cld * 0.01
  met_in(i)%wind = daily(i,doy)%wnd
  
  if (met_in(i)%tmin > met_in(i)%tmax) then
    write(0,*)'unphysical tmin tmax'
    write(0,*)met_in(i)%tmin,met_in(i)%tmax
  end if

  ! initialize daily meteorology for the current day, NB this is the only time this is called for the current day

  call weathergen(met_in(i),dmet0(i))
  
  ! write(*,'(i5,2f8.2)')doy,dmet0(i)%tmin,dmet0(i)%tmax
  
  ! calculate top-of-the-atmosphere insolation and daylength
  
  call insol(slon,orbit,latr,solar)
  
  pixel(i)%phi   = solar%phi
  dmet0(i)%delta = solar%delta
  dmet0(i)%rad0  = solar%rad0
  dmet0(i)%dayl  = solar%dayl
  
  dmet0(i)%Bsw   = B0  ! initialize to background albedo

end do

! ----------------------------------------------------------------------------------------------------------------
! initialization year loop

write(0,*)'start initialization daily loop'

doy = 1  ! day of year counter

do m = 1,nmos
  do d = 1,ndm(m)
  
    d1 = doy + 1
    if (d1 > ndy) d1 = 1
  
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

      met_in(i)%tmin = daily(i,d1)%tmin
      met_in(i)%tmax = daily(i,d1)%tmax
      met_in(i)%cldf = daily(i,d1)%cld * 0.01
      met_in(i)%wind = daily(i,d1)%wnd
      
      if (met_in(i)%tmin > met_in(i)%tmax) then
        write(0,*)'unphysical tmin tmax'
        write(0,*)met_in(i)%tmin,met_in(i)%tmax
      end if

      ! generate daily meteorology for the next day (to have tmin for nighttime temperature)

      call weathergen(met_in(i),dmet1(i))
      
      ! calculate top-of-the-atmosphere insolation and daylength
      
      call insol(slon,orbit,latr,solar)
      
      dmet0(i)%delta = solar%delta
      dmet0(i)%rad0  = solar%rad0
      dmet0(i)%dayl  = solar%dayl
            
      ! calculate integrated day- and night-time temperature
      ! nighttime is goes into the next day, need to know tmin of the next day
      ! in situations of polar day, day temperature covers 23 hrs
      ! in situations of polar night, day temperature covers 1 hr

      call diurnaltemp(dmet0(i),dmet1(i))
      
      ! surface radiation budget and potential evapotranspiration
  
      call radpet(pixel(i),dmet0(i))
      
      ! snow dynamics

      call snow(pixel(i),dmet0(i))
      
      ! soil water balance, including actual evapotranspiration and alpha
      
      call soilwater(dmet0(i),soilw(i))
      
!       write(0,*)m,d,dmet0(i)%tday,dmet0(i)%tnight,dmet0(i)%prec,dmet0(i)%snow,dmet0(i)%melt,dmet0(i)%swe,dmet0(i)%fsnow,dmet0(i)%asnow,dmet0(i)%Bsw
      
      ! write(0,*)soilw(i)

      ! store today's meteorology for tomorrow

      dmet1(i)%swe   = dmet0(i)%swe
      dmet1(i)%asnow = dmet0(i)%asnow
      dmet1(i)%Bsw   = dmet0(i)%Bsw

      dmet0(i) = dmet1(i)

    end do ! cells

    doy = doy + 1
    
    if (doy > ndy) doy = 1
    
  end do   ! days in the month
end do     ! month

! ----------------------------------------------------------------------------------------------------------------
! computation daily loop

write(0,*)
write(0,*)'start computation daily loop'

mmet%mpet  = 0.
mmet%alpha = 0.
mmet%direct  = 0.
mmet%diffuse = 0.

pixel%gdd0 = 0.
pixel%gdd5 = 0.

doy = 1

do m = 1,nmos
  do d = 1,ndm(m)
  
    d1 = doy + 1
    if (d1 > ndy) d1 = 1

    write(status_line,'(a,i0,a,i0)')' working on ',m,' ',d
    call overprint(status_line)

    ! calculate true solar longitude (position in earth's orbit around the sun) for daylength, 
    ! only need to do this once per day, not for each cell
    
    slon = truelon(orbit,cal,doy)
  
    ! loop over valid cells

    do i = 1,ncells
    
      x = pixel(i)%x
      y = pixel(i)%y
    
      met_in(i)%prec = climate(x,y,m)%pre
      met_in(i)%wetf = climate(x,y,m)%wet
      met_in(i)%wetd = climate(x,y,m)%wet * real(present_mon_noleap(m))
      
      met_in(i)%tmin = daily(i,d1)%tmin
      met_in(i)%tmax = daily(i,d1)%tmax
      met_in(i)%cldf = daily(i,d1)%cld * 0.01
      met_in(i)%wind = daily(i,d1)%wnd

      ! generate daily meteorology for next day
  
      call weathergen(met_in(i),dmet1(i))
      
      ! calculate daylength and insolation
      
      latr = coords(x,y)%geolat * pir
      
      call insol(slon,orbit,latr,solar)
      
      dmet0(i)%delta = solar%delta
      dmet0(i)%rad0  = solar%rad0
      dmet0(i)%dayl  = solar%dayl
      
      ! calculate integrated day- and night-time temperature

      call diurnaltemp(dmet0(i),dmet1(i))

      ! surface radiation budget and potential evapotranspiration
  
      call radpet(pixel(i),dmet0(i))
      
      ! snow dynamics
            
      call snow(pixel(i),dmet0(i))

      ! soil water balance, including actual evapotranspiration and alpha
      
      call soilwater(dmet0(i),soilw(i))

!       write(0,*)m,d,dmet0(i)%tday,dmet0(i)%tnight,dmet0(i)%prec,dmet0(i)%snow,dmet0(i)%melt,dmet0(i)%swe,dmet0(i)%fsnow,dmet0(i)%asnow,dmet0(i)%Bsw
      
      ! monthly summaries

      mmet(i,m)%direct  = mmet(i,m)%direct  + dmet0(i)%rdirect  / real(ndm(m))
      mmet(i,m)%diffuse = mmet(i,m)%diffuse + dmet0(i)%rdiffuse / real(ndm(m))

      mmet(i,m)%mpet  = mmet(i,m)%mpet  + dmet0(i)%dpet
      mmet(i,m)%alpha = mmet(i,m)%alpha + dmet0(i)%alpha / real(ndm(m))

      pixel(i)%gdd5 = pixel(i)%gdd5 + max(dmet0(i)%tday - 5.,0.)
      pixel(i)%gdd0 = pixel(i)%gdd0 + max(dmet0(i)%tday,0.)

      ! store today's meteorology for tomorrow

      dmet1(i)%swe   = dmet0(i)%swe
      dmet1(i)%asnow = dmet0(i)%asnow
      dmet1(i)%Bsw   = dmet0(i)%Bsw
      
      dmet0(i) = dmet1(i)

    end do ! cells
    
    doy = doy + 1
    
    if (doy > ndy) doy = 1
    
  end do   ! days in the month
end do     ! month

write(0,*)

! ----
! calculate annual alpha and biome

do i = 1,ncells

  pixel(i)%aalpha = sum(mmet(i,:)%alpha) / 12.
  
  wm = pixel(i)%wm
  
  pixel(i)%awm = mmet(i,wm)%alpha

  call calcbiome(pixel(i))

end do

! ---------------------------------
! write model output

write(0,*)'writing'

call writereal3d(ofid,gridinfo,pixel,'rdirect',mmet%direct)
call writereal3d(ofid,gridinfo,pixel,'rdiffuse',mmet%diffuse)

call writereal3d(ofid,gridinfo,pixel,'mpet',mmet%mpet)
call writereal3d(ofid,gridinfo,pixel,'alpha',mmet%alpha)

call writereal2d(ofid,gridinfo,pixel,'tcm',pixel%tcm)
call writereal2d(ofid,gridinfo,pixel,'twm',pixel%twm)
call writereal2d(ofid,gridinfo,pixel,'GDD0',pixel%gdd0)
call writereal2d(ofid,gridinfo,pixel,'GDD5',pixel%gdd5)
call writereal2d(ofid,gridinfo,pixel,'awm',pixel%awm)
call writereal2d(ofid,gridinfo,pixel,'aalpha',pixel%aalpha)

call writeinteger2d(ofid,gridinfo,pixel,'biome',pixel%biome)

! ---------------------------------

call closeoutput(ofid)

! ---------------------------------
! finish

end program biome1
