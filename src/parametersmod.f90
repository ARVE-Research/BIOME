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

end module parametersmod
