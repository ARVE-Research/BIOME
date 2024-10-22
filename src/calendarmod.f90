module calendarmod

use parametersmod, only : dp

implicit none

integer, parameter :: nm = 12              ! number of months in the year
integer, parameter :: nd_360 = 360         ! number of days in a 360-day year
integer, parameter :: daysinmonth360 = 30  ! number of days in a month of a 360-day year
integer, parameter :: nd_365 = 365         ! number of days in a 365-day "noleaps" year
integer, parameter :: nd_366 = 366         ! number of days in a 366-day leap year

! other calendar-related variables

real(dp) :: veqday_360    = 80._dp           ! fixed vernal equinox day, 360-day year
real(dp) :: veqday_365    = 80.5_dp          ! fixed vernal equinox day, 365-day year
real(dp) :: veqday_366    = 81.5_dp          ! fixed vernal equinox day, 366-day year
real(dp) :: ssday_360     = 170.5_dp         ! fixed (northern) summer solstice day, 360-day year
real(dp) :: ssday_365     = 173._dp          ! fixed (northern) summer solstice day, 365-day year
real(dp) :: ssday_366     = 173.5_dp         ! fixed (northern) summer solstice day, 366-day year
real(dp) :: tropical_year = 365.24219876_dp  ! length of a tropical year (days)
real(dp) :: progreg_year  = 365.2425_dp      ! length of a Gregorian year (days)

integer  :: nd_progreg                       ! number of days in a 365 or 366-day year proleptic_gregorian calendar
real(dp) :: veqday_progreg                   ! vernal equinox day in a 365 or 366-day year proleptic_gregorian calendar
real(dp) :: midMarch                         ! mid-March day
real(dp) :: perihelion                       ! perihelion day (in a proleptic Gregorian calendar)
real(dp) :: aphelion                         ! aphelion day (in a proleptic Gregorian calendar)

! month-length definitions

real(dp) :: present_mon_360(nm) = 30._dp     ! present-day month lengths in 360-day year

real(dp) :: present_beg_360(nm) = &          ! present-day month beginning day in 360-day year
    [  0._dp, 30._dp, 60._dp, 90._dp, 120._dp, 150._dp, 180._dp, 210._dp, 240._dp, 270._dp, 300._dp, 330._dp ]

real(dp) :: present_mid_360(nm) = &          ! present-day month middle day in 360-day year
    [ 15._dp, 45._dp, 75._dp, 105._dp, 135._dp, 165._dp, 195._dp, 225._dp, 255._dp, 285._dp, 315._dp, 345._dp ]

real(dp) :: present_end_360(nm) = &          ! present-day month ending day in 360-day year
    [ 30._dp, 60._dp, 90._dp, 120._dp, 150._dp, 180._dp, 210._dp, 240._dp, 270._dp, 300._dp, 330._dp, 360._dp ]   

real(dp) :: present_mon_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
    [ 31._dp, 28._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp ]

real(dp) :: present_beg_365(nm) = &          ! present-day month beginning day in 365-day (noleap) year
    [  0._dp, 31._dp, 59._dp, 90._dp, 120._dp, 151._dp, 181._dp, 212._dp, 243._dp, 273._dp, 304._dp, 334._dp ]

real(dp) :: present_mid_365(nm) = &          ! present-day month middle day in 365-day (noleap) year
    [ 15.5_dp, 45._dp, 74.5_dp, 105._dp, 135.5_dp, 166._dp, 196.5_dp, 227.5_dp, 258._dp, 288.5_dp, 319._dp, 349.5_dp ]

real(dp) :: present_end_365(nm) = &          ! present-day month ending day in 365-day (noleap) year
    [ 31._dp, 59._dp, 90._dp, 120._dp, 151._dp, 181._dp, 212._dp, 243._dp, 273._dp, 304._dp, 334._dp, 365._dp ]    

real(dp) :: present_mon_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
    [ 31._dp, 29._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp ]

real(dp) :: present_beg_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
    [  0._dp, 31._dp, 60._dp, 91._dp, 121._dp, 152._dp, 182._dp, 213._dp, 244._dp, 274._dp, 305._dp, 335._dp ]

real(dp) :: present_mid_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
    [ 15.5_dp, 45.5_dp, 75.5_dp, 106._dp, 136.5_dp, 167._dp, 197.5_dp, 228.5_dp, 259._dp, 289.5_dp, 320._dp, 350.5_dp ]

real(dp) :: present_end_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
    [ 31._dp, 60._dp, 91._dp, 121._dp, 152._dp, 182._dp, 213._dp, 244._dp, 274._dp, 305._dp, 335._dp, 366._dp ]   

real(dp) :: present_mon_365_trop(nm) = &     ! present-day month lengths in a tropical year (note Feb.)
    [ 31._dp, 28.24219876_dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp ]

real(dp) :: present_beg_365_trop(nm) = &     ! present-day month beginning day in a tropical year
    [  0.0000_dp, 31.0000_dp, 59.2422_dp, 90.2422_dp, 120.2422_dp, 151.2422_dp, 181.2422_dp, 212.2422_dp, 243.2422_dp, &
        273.2422_dp, 304.2422_dp, 334.2422_dp ]

real(dp) :: present_mid_365_trop(nm) = &     ! present-day month middle day in a tropical year
    [ 15.5000_dp, 45.1211_dp, 74.7422_dp, 105.2422_dp, 135.7422_dp, 166.2422_dp, 196.7422_dp, 227.7422_dp, 258.2422_dp, &
        288.7422_dp, 319.2422_dp, 349.7422_dp ]

real(dp) :: present_end_365_trop(nm) = &     ! present-day month ending day in a tropical year
    [ 31.0000_dp, 59.2422_dp, 90.2422_dp, 120.2422_dp, 151.2422_dp, 181.2422_dp, 212.2422_dp, 243.2422_dp, 273.2422_dp, &
        304.2422_dp, 334.2422_dp, 365.2422_dp ]    

real(dp) :: present_mon_365_progreg(nm) = &  ! present-day month lengths in a Gregorian year (note Feb.)
    [ 31._dp, 28.2425_dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp, 31._dp, 30._dp, 31._dp, 30._dp, 31._dp ]

real(dp) :: present_beg_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
    [  0.0000_dp, 31.0000_dp, 59.2425_dp, 90.2425_dp, 120.2425_dp, 151.2425_dp, 181.242_dp, 212.2425_dp, 243.2425_dp, &
        273.2425_dp, 304.2425_dp, 334.2425_dp ]

real(dp) :: present_mid_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
    [ 15.5000_dp, 45.1213_dp, 74.7425_dp, 105.2425_dp, 135.7425_dp, 166.2425_dp, 196.7425_dp, 227.7425_dp, 258.2425_dp, &
    288.7425_dp, 319.2425_dp, 349.7425_dp ]

real(dp) :: present_end_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
    [ 31.0000_dp, 59.2425_dp, 90.2425_dp, 120.2425_dp, 151.2425_dp, 181.2425_dp, 212.2425_dp, 243.2425_dp, 273.2425_dp, &
    304.2425_dp, 334.2425_dp, 365.2425_dp ]

end module calendarmod
