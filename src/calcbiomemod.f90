module calcbiomemod

implicit none

contains

! -----------------------------------------------------

subroutine calcbiome(pixel)

use parametersmod, only : sp,i2
use typesmod,      only : pixeltype

implicit none

! This subroutine calculates the BIOME

! Choosing the plant types that are present
! The planttype array has 13 positions which the following plant types occupying the following positions:
!  1 Tropical Evergreen
!  2 Tropical Raingreen
!  3 Warm-Temperate Evergreen
!  4 Temperate Summergreen
!  5 Cool-Temperate Conifer
!  6 Boreal Evergreen Conifer
!  7 Boreal Summergreen
!  8 Sclerophyll/succulent
!  9 Warm grass/shrub
! 10 Cool grass/shrub
! 11 Cold grass/shrub
! 12 Hot desert shrub
! 13 Cold desert shrub

! --- list of biome number/names ---
!  1 tropical rain forest
!  2 tropical seasonal forest
!  3 tropical dry forests/savanna
!  4 warm mixed forest
!  5 temperate deciduous forest
!  6 cool mixed forest
!  7 cool conifer forest
!  8 cold evergreen forest
!  9 cold mixed forest
! 10 cold deciduous forest
! 11 xerophytic woods and scrub
! 12 warm grass/shrub
! 13 cool grass/shrub
! 14 tundra
! 15 hot desert
! 16 semidesert
! 17 polar desert

! argument

type(pixeltype), intent(inout) :: pixel

! local variables

real(sp) :: alpha
real(sp) :: tcm
real(sp) :: twm
real(sp) :: awm
real(sp) :: acm
real(sp) :: GDD
real(sp) :: GDD0
integer  :: biome
integer  :: plantcase

logical, dimension(13) :: planttype

! ---------------------------------

alpha = pixel%aalpha

tcm = pixel%tcm
twm = pixel%twm
awm = pixel%awm
acm = pixel%acm
 
GDD = pixel%gdd5

GDD0 = pixel%gdd0

planttype = .false.


! plant types: 1) tropical evergreen, 2) tropical raingreen
! EDIT: calc independently so both can be true simultaneously
if (tcm >= 15.5) then
  if (alpha >= 0.80) then
    planttype(1) = .true.
  end if
  if (alpha >= 0.45 .and. alpha <= 0.95) then
    planttype(2) = .true.
  end if
end if

! plant type: 3) warm-temperate evergreen
! 0.65-->.33
! Include climates with very wet winters and dry summers
if (tcm >= 5) then
   if (alpha >= 0.33 .and. acm >= 0.6) then
!   if (alpha >= 0.33 .and. acm >= .98 .and. awm >= 0.008) then
    planttype(3) = .true.
  end if
end if

!Planttypes #4 and #5 I may want to lower their alpha values to 0.30, so that they can be true in pixel cells destined to have decidious forests. 

! plant type: 4) temperate summergreen
! 0.65-->.33
if (tcm >= -15 .and. tcm <=15.5) then
  if (GDD >= 1200) then
  if (alpha >= 0.33) then
    planttype(4) = .true.
  end if
  end if
end if

! plant type: 5) cool-temp conifer
! 0.65-->.33
if (tcm >= -19 .and. tcm <= 5) then
  if (GDD >= 900) then
  if (alpha >= 0.33) then
    planttype(5) = .true.
  end if
  end if
end if

! plant type: 6) boreal evergreen conifer
! 0.75-->.38
if (tcm >= -35 .and. tcm <= -2) then
  if (GDD >= 350) then
  if (alpha >= 0.38) then
    planttype(6) = .true.
  end if
  end if
end if 

! plant type: 7) boreal summergreen
! 0.65-->.21
if (tcm <= 5) then 
  if (GDD >= 350) then
  if (alpha >= 0.33) then
  planttype(7) = .true.
  end if
  end if
end if

! plant type: 8) sclerophyll/succulent
!0.28-->.14-->.09 ??????????????????????????????????????
if (tcm >= 5) then
  if (alpha >= 0.09) then
  planttype(8) = .true.
  end if
end if

! plant type: 9) warm grass/shrub
! 0.18-->.06
if (twm >= 22) then
  if (alpha >= 0.06) then
  planttype(9) = .true.
  end if
end if

! plant type: 10) cool grass/shrub
! 0.33-->.11-->.08
if (GDD >= 500) then
  if (alpha >= 0.08) then
  planttype(10) = .true.
  end if
end if

! plant type: 11) cold grass/shrub
! 0.33-->.11-->.08
if (GDD0 >= 100) then
  if (alpha >= 0.08) then
  planttype(11) = .true. 
  end if
end if

! plant type: 12) hot desert shrub
if (twm >= 22) then
  planttype(12) = .true.
end if

! plant type: 13) cold desert shrub
if (GDD0 >= 100) then
  planttype(13) = .true.
end if

! BIOME Determining using present plant types. 

! Here is my thinking. The select case only works using scalars. Since the planttype is an array here is what I will perform. 
! I will set each .true./.false. combination of the planttype array equal to a integer value in the integer scalar "plantcase". Then I will perform a select case(plantcase).

! write(0,*)alpha,tcm,twm,GDD,GDD0
! write(0,*)planttype
! read(*,*)

! dominance class 1

if (planttype(1)) then
  if (.not. planttype(2)) then
    pixel%biome = 1
  else
    pixel%biome = 2
  end if
  return
end if

if (planttype(2) .and. .not. planttype(1)) then 
  pixel%biome = 3
  return
end if

! dominance class 2

if (planttype(3)) then
  pixel%biome = 4
  return
end if 

! dominance class 3; planttypes #4-#7

if (planttype(4) .and. planttype(5) .and. planttype(7)) then 
  if (.not. planttype(6)) then
  pixel%biome = 5
  else
 pixel%biome = 6
  end if
  return
end if 

if (planttype(5) .and. planttype(6) .and. planttype(7) .and. .not. planttype(4)) then 
  pixel%biome = 7
  return
end if

if (planttype(6) .and. planttype(7) .and. .not. any(planttype(4:5))) then 
  pixel%biome = 8
  return
end if

if (planttype(5) .and. planttype(7) .and. .not. planttype(4) .and. .not. planttype(6)) then
  pixel%biome = 9
  return
end if 

if (planttype(7) .and. .not. any(planttype(4:6))) then 
  pixel%biome = 10
  return
end if


! I am swapping the order of dominance class #4 and #5 because biome 11 is being inferred in areas where biome #12 should be present. 

! dominance class 4; planttype #8

if (planttype(8)) then
  pixel%biome = 11
  return
end if

 ! domiance class 5; planttype #9

if (planttype(9)) then 
  pixel%biome = 12
  return
end if

! dominance class 6; planttypes #10 #11

! if (planttype(11)) then
!   if (.not. planttype(10)) then 
!     pixel%biome = 14
!   else 
!     pixel%biome = 13
!   end if
!   return
! end if

! dominance class 6; planttypes #10 #11
! EDIT handle cases where only cool grass/shrub is present

if (planttype(10) .or. planttype(11)) then
  if (planttype(11) .and. .not. planttype(10)) then 
    pixel%biome = 14
  else 
    pixel%biome = 13
  end if
  return
end if

! dominance class 7; planttype #12

if (planttype(12)) then
  pixel%biome = 15
  return
end if

! dominance class 8; planttype #13

if (planttype(13)) then
  pixel%biome = 16
  return
end if

! Biome 17 for polar desert when no plant types present
pixel%biome = 17

end subroutine calcbiome

end module calcbiomemod