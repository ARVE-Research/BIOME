module calcbiomemod

implicit none

contains

!---------------------------------------------------------------------------------------

subroutine calcbiome(pixel)

use parametersmod, only : sp,i2
use typesmod,      only : pixeltype

implicit none

! This subroutine calculates the BIOME

! Choosing the plant types that are present
! The planttype array has 12 positions which the following plant types occupying the following positions:
!  1 Warm-Temperate Evergreen
!  2 Coastal evergreen (fog-belt forest)
!  3 Temperate Summergreen
!  4 Cool-Temperate Conifer
!  5 Boreal Evergreen Conifer
!  6 Boreal Summergreen
!  7 Sclerophyll/succulent (xerophytic, residual)
!  8 Chaparral sclerophyll (warm-winter coastal)
!  9 Valley savanna grass/oak (Central Valley)
! 10 Warm grass/shrub
! 11 Cool grass/shrub
! 12 Hot desert shrub

! --- list of biome number/names ---
!  1 mixed oak and pine forest
!  2 coastal forest
!  3 yellow pine forest
!  4 montane forest
!  5 upper montane forest
!  6 tundra
!  7 pinyon-juniper woodland
!  8 subalpine woodland
!  9 arid scrub
! 10 coastal scrub
! 11 oak savanna
! 12 sagebrush steppe
! 13 desert steppe
! 14 hot desert

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

logical, dimension(12) :: planttype

!---------------------------------------------------------------------------------------

alpha = pixel%aalpha

tcm = pixel%tcm
twm = pixel%twm
awm = pixel%awm
acm = pixel%acm
 
GDD = pixel%gdd5

GDD0 = pixel%gdd0

planttype = .false.

!---------------------------------------------------------------------------------------

! plant type: 1) warm-temperate evergreen
!Prentice Alpha Value = 0.65
! Include climates with very wet winters and dry summers
if (tcm >= 5) then
   if (alpha >= 0.45 .or. (alpha >= 0.4 .and. twm <= 22)) then
    planttype(1) = .true.
  end if
end if

! plant type: 2) coastal evergreen (fog-belt forest)
! Redwood, Douglas fir, tanoak — high year-round moisture with mild winters
! Subset of warm-temperate evergreen climate space but wetter
if (tcm > 8.0 .and. twm < 23.0) then
  if (alpha > 0.55) then
    planttype(2) = .true.
  end if
end if

! plant type: 3) temperate summergreen
!Prentice Alpha Value = 0.65
if (tcm >= -15 .and. tcm <=15.5) then
  if (GDD >= 1200) then
  if (alpha >= 0.33) then
    planttype(3) = .true.
  end if
  end if
end if

! plant type: 4) cool-temp conifer
!Prentice Alpha Value = 0.65
if (tcm >= -19 .and. tcm <= 5) then
  if (GDD >= 900) then
  if (alpha >= 0.25) then
    planttype(4) = .true.
  end if
  end if
end if

! plant type: 5) boreal evergreen conifer
!Prentice Alpha Value = 0.75
if (tcm >= -35 .and. tcm <= -2) then
  if (GDD >= 350) then
  if (alpha >= 0.38) then
    planttype(5) = .true.
  end if
  end if
end if 

! plant type: 6) boreal summergreen
!Prentice Alpha Value = 0.65
if (tcm <= 5) then 
  if (GDD >= 350) then
  if (alpha >= 0.25) then
  planttype(6) = .true.
  end if
  end if
end if

! plant type: 7) sclerophyll/succulent (xerophytic residual)
!Prentice Alpha Value = 0.28
if (tcm >= 5) then
  if (alpha >= 0.12) then
  planttype(7) = .true.
  end if
end if

! plant type: 8) chaparral sclerophyll (warm-winter coastal)
! New biome, not in Prentice paper
! Warm winters (tcm > 12) separate from inland scrub; alpha > 0.20
if (tcm > 12.0) then
  if (alpha > 0.20) then
    planttype(8) = .true.
  end if
end if

! plant type: 9) valley savanna grass/oak (Central Valley)
! New biome, not in Prentice paper
! Moderate moisture (alpha 0.28-0.48) with cooler winters (tcm < 9)
if (tcm < 11.0) then
  if (alpha >= 0.28 .and. alpha <= 0.48) then
    planttype(9) = .true.
  end if
end if

! plant type: 10) warm grass/shrub
!Prentice Alpha Value = 0.18
if (twm >= 22) then
  if (alpha >= 0.09) then
  planttype(10) = .true.
  end if
end if

! plant type: 11) cool grass/shrub
!Prentice Alpha Value = 0.33
! Now also covers old cold grass/shrub and tundra pixels
if (GDD >= 500 .or. GDD0 >= 100) then
  if (alpha >= 0.08) then
  planttype(11) = .true.
  end if
end if

! plant type: 12) hot desert shrub
! Now also covers old cold desert shrub pixels
if (twm >= 22 .or. GDD0 >= 100) then
  planttype(12) = .true.
end if

!---------------------------------------------------------------------------------------
! BIOME Determining using present plant types. 
!---------------------------------------------------------------------------------------

! dominance class 1; planttypes #2, #1

! Coastal forest checked first — more specific, wetter type
if (planttype(2)) then
  pixel%biome = 2   ! coastal forest
  return
end if

if (planttype(1)) then
  pixel%biome = 1   ! mixed oak and pine forest
  return
end if 

! dominance class 2; planttypes #3-#6

if (planttype(3) .and. planttype(4) .and. planttype(6)) then 
  if (.not. planttype(5)) then
  pixel%biome = 3   ! yellow pine forest
  else
 pixel%biome = 4    ! montane forest
  end if
  return
end if 

if (planttype(4) .and. planttype(5) .and. planttype(6) .and. .not. planttype(3)) then 
  pixel%biome = 5   ! upper montane forest
  return
end if

if (planttype(5) .and. planttype(6) .and. .not. any(planttype(3:4))) then 
  pixel%biome = 6   ! tundra
  return
end if

if (planttype(4) .and. planttype(6) .and. .not. planttype(3) .and. .not. planttype(5)) then
  pixel%biome = 7   ! pinyon-juniper woodland
  return
end if 

if (planttype(6) .and. .not. any(planttype(3:5))) then 
  pixel%biome = 8   ! subalpine woodland
  return
end if

! dominance class 3; planttypes #7, #8, #9
! Check specific sclerophyll subtypes before falling back to generic arid scrub

if (planttype(8) .or. planttype(9) .or. planttype(7)) then

  ! Coastal scrub dominates where warm-winter coastal sclerophyll is present
  if (planttype(8)) then
    pixel%biome = 10  ! coastal scrub
    return
  end if

  ! Valley savanna where moderate moisture + cool winters on valley floor
  if (planttype(9)) then
    pixel%biome = 11  ! oak savanna
    return
  end if

  ! Residual arid scrub — dry interior
  pixel%biome = 9    ! arid scrub
  return

end if

 ! dominance class 4; planttype #10

if (planttype(10)) then 
  pixel%biome = 13   ! desert steppe
  return
end if

! dominance class 5; planttype #11

if (planttype(11)) then
  pixel%biome = 12   ! sagebrush steppe
  return
end if

! dominance class 6; planttype #12

if (planttype(12)) then
  pixel%biome = 14   ! hot desert
  return
end if

! Fallback: hot desert when no plant types present
pixel%biome = 14

!---------------------------------------------------------------------------------------

end subroutine calcbiome

end module calcbiomemod