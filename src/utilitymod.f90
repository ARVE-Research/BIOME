module utilitymod

implicit none

public :: pos
public :: matsol

interface pos
  module procedure pos_sp,pos_dp
end interface pos

contains

! ------------------------

integer function pos_dp(vect,val)

! finds the index in the vector "vect" that has a value nearest to the input scalar "val"

use parametersmod, only : dp

implicit none

! arguments

real(dp), dimension(:), intent(in) :: vect
real(dp),               intent(in) :: val

! ---

pos_dp = minloc(abs(vect - val),dim=1)

end function pos_dp

! ------------------------

integer function pos_sp(vect,val)

! finds the index in the vector "vect" that has a value nearest to the input scalar "val"

use parametersmod, only : sp

implicit none

! arguments

real(sp), dimension(:), intent(in) :: vect
real(sp),               intent(in) :: val

! ---

pos_sp = minloc(abs(vect - val),dim=1)

end function pos_sp

! ------------------------

subroutine matsol(mat,sol)

! Provides matrix solution to X matrix in a A * X = B system using LU decomposition
! Code adapted from Press et al. (1996) 
! Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing 2nd Edition (P.1016-1017)
! Coded by Leo Lai ca. 2021

use parametersmod, only : sp

implicit none

! arguments

real(sp), dimension(:,:), intent(inout) :: mat
real(sp), dimension(:),   intent(inout) :: sol

! parameter

real(sp), parameter :: tiny_sp = 1.0e-38_sp

! local variables

integer,  dimension(size(sol)) :: indx
real(sp), dimension(size(sol)) :: mv

real(sp), allocatable, dimension(:,:) :: prod
integer,               dimension(1)   :: maxl

integer  :: i
integer  :: n
integer  :: k
integer  :: ll
integer  :: max
real(sp) :: summ

! -------------------

n = size(sol)

mv = 1. / maxval(abs(mat), dim=2)

!---

do i = 1, n

  maxl = maxloc(mv(i:n) * abs(mat(i:n,i)))

  max = (i - 1) + maxl(1)

  indx(i) = max

  !---

  if (mat(i,i) == 0.) mat(i,i) = tiny_sp

  mat(i+1:n, i) = mat(i+1:n, i) / mat(i,i)

  !---

  allocate(prod(i+1:n, i+1:n))

  prod = spread(mat(i+1:n, i), dim=2, ncopies=size(mat(i, i+1:n)))

  prod = prod * spread(mat(i, i+1:n), dim=1, ncopies=size(mat(i+1:n, i)))

  !---

  mat(i+1:n, i+1:n) = mat(i+1:n, i+1:n) - prod

  deallocate(prod)

end do

!---

k = 0

do i = 1, n

  ll = indx(i)
  summ = sol(ll)
  sol(ll) = sol(i)

  if (k /= 0) then

    summ = summ - dot_product(mat(i, k:i-1), sol(k:i-1))

  else if (summ /= 0.) then

    k = i

  end if

  sol(i) = summ

end do

!---

do i = n, 1, -1

  sol(i) = (sol(i) - dot_product(mat(i, i+1:n), sol(i+1:n))) / mat(i,i)

end do

end subroutine matsol

! ------------------------

end module utilitymod
