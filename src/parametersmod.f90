module parametersmod

use iso_fortran_env

implicit none

integer, parameter :: i1 = int8    ! 1 byte integer
integer, parameter :: i2 = int16   ! 2 byte integer
integer, parameter :: i4 = int32   ! 4 byte integer
integer, parameter :: i8 = int64   ! 8 byte integer
integer, parameter :: sp = real32  ! 4 byte real
integer, parameter :: dp = real64  ! 8 byte real

integer, parameter :: stdin  = input_unit
integer, parameter :: stdout = output_unit
integer, parameter :: stderr = error_unit

real(sp), parameter :: rmissing = -9999.

real(sp), parameter :: hsp = huge(1._sp)    ! largest positive 4-byte real

real(sp), parameter :: Tfreeze = 273.15 ! freezing temperature of freshwater (K)

real(dp), parameter :: pi  = 3.14159265358979323846_dp !26433 83279 50288 41971 69399 37510 (unitless)
real(dp), parameter :: pir = pi / 180._dp

end module parametersmod
