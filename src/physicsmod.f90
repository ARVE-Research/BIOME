module physicsmod

implicit none

public :: pw
public :: gamma
public :: lvap
public :: Econ
public :: Esat
public :: desdT

contains

! ----------------------------------------------------------------------------------------------------------------

real(sp) function pw(T)

! function to calculate the temperature-dependent density of liquid water, based on:
! Kell, G. S. (1975). Density, Thermal Expansivity, and Compressibility of Liquid Water from 0° to 150°C: 
! Correlations and Tables for Atmospheric Pressure and Saturation Reviewed and Expressed on 1968 Temperature Scale. 
! Journal of Chemical and Engineering Data, 20(1), 97-105. 

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: T ! water temperature (degC)

! parameters

real(sp), parameter :: a = 999.83952
real(sp), parameter :: b =  16.945176
real(sp), parameter :: c =   7.9870401e-3
real(sp), parameter :: d =  46.170461e-6
real(sp), parameter :: e = 105.56302e-9
real(sp), parameter :: f = 280.54253e-12
real(sp), parameter :: g =  16.879850e-3

! ----

pw = (a + b * T - c * T**2 - d * T**3 + e * T**4 + f * T**5) / (1. + g * T)  ! Kell (1975) eqn. 16

end function pw

! ----------------------------------------------------------------------------------------------------------------

real(sp) function gamma(Tair)

! function to calculate the temperature-dependent psychrometer constant (Pa K-1)

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: Tair  ! air temperature (degC)

! ----

gamma = 65.05 + Tair * 0.064  ! psychrometer constant

end function gamma

! ----------------------------------------------------------------------------------------------------------------

real(sp) function lvap(Tair)

! function to calculate the temperature-dependent latent heat of vaporization of water (kJ kg-1), based on:
! Henderson-Sellers, B. (1984). A new formula for latent heat of vaporization of water as a function of temperature.
! Quarterly Journal of the Royal Meteorological Society, 110(466), 1186-1190. doi:10.1002/qj.49711046626

use parametersmod, only : sp,tfreeze

implicit none

! argument

real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variable

real(sp) :: Tk  ! air temperature (degC)

! ----

Tk = Tair + Tfreeze

lvap = 0.001 * 1.91846e6 * (Tk / (Tk - 33.91))**2  ! Henderson-Sellers (1984), eqn. 1

end function lvap

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Econ(Tair)

! energy-to-water conversion factor, Sandoval et al. (2024) eqn 51

use parametersmod, only : sp,tfreeze

implicit none

! argument

real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variables

real(sp) :: ss  ! rate of change of saturation vapor pressure with temperature (Pa K-1)

! ----

ss = desdT(Tair + tfreeze)

Econ = ss / (lvap(Tair) * pw(Tair) * (ss + 0.24 * gamma(Tair)))

end function Econ

! ----------------------------------------------------------------------------------------------------------------

real(sp) function esat(temp)
  
! Function to calculate saturation vapor pressure in water and ice
! From CLM formulation, table 5.2, after Flatau et al. 1992
! calculated in double precision because temperatures close to zero cause underflow

use parametersmod, only : sp,dp,tfreeze

implicit none

! argument

real(sp), intent(in) :: temp  ! temperature in K

! parameters

! coefficients for liquid water

real(dp), dimension(9), parameter :: al = [ 6.11213476     ,  &
                                            4.44007856e-1  ,  &
                                            1.43064234e-2  ,  &
                                            2.64461437e-4  ,  &
                                            3.05903558e-6  ,  &
                                            1.96237241e-8  ,  &
                                            8.92344772e-11 ,  &
                                           -3.73208410e-13 ,  &
                                            2.09339997e-16 ]

! coefficients for ice

real(dp), dimension(9), parameter :: ai = [ 6.11123516     ,  &
                                            5.03109514e-1  ,  &
                                            1.88369801e-2  ,  &
                                            4.20547422e-4  ,  &
                                            6.14396778e-6  ,  &
                                            6.02780717e-8  ,  &
                                            3.87940929e-10 ,  &
                                            1.49436277e-12 ,  &
                                            2.62655803e-15 ]

integer, dimension(9), parameter :: p = [0,1,2,3,4,5,6,7,8]

! local variables

real(dp) :: esatdp
real(dp) :: T  

real(dp), dimension(9) :: a

! ----

if (temp <= tfreeze) then   ! these coefficients are for temperature values in Celcius
  a = ai
else
  a = al
end if

T = temp - tfreeze

esatdp = sum(a * T**p)

esat = real(100._dp * esatdp)
   
end function esat

! ----------------------------------------------------------------------------------------------------------------

real(sp) function desdT(temp)

! Function to calculate the first derivative of saturation vapor pressure in water and ice vs. temperature
! From CLM formulation, table 5.3, after Flatau et al. 1992

use parametersmod, only : sp,dp,tfreeze

implicit none

! argument

real(sp), intent(in) :: temp ! temperature in K

! parameters

! coefficients for liquid water

real(sp), dimension(9), parameter :: bl = [ 4.44017302e-1  ,  &
                                            2.86064092e-2  ,  &
                                            7.94683137e-4  ,  &
                                            1.21211669e-5  ,  &
                                            1.03354611e-7  ,  &
                                            4.04125005e-10 ,  &
                                           -7.88037859e-13 ,  &
                                           -1.14596802e-14 ,  &
                                            3.81294516e-17 ]

! coefficients for ice

real(sp), dimension(9), parameter :: bi = [ 5.03277922e-1  ,  &
                                            3.77289173e-2  ,  &
                                            1.26801703e-3  ,  &
                                            2.49468427e-5  ,  &
                                            3.13703411e-7  ,  &
                                            2.57180651e-9  ,  &
                                            1.32268878e-11 ,  &
                                            3.94116744e-14 ,  &
                                            4.98070196e-17 ]

integer, dimension(9), parameter :: p = [0,1,2,3,4,5,6,7,8]

! local variables

real(dp), dimension(9) :: b ! coefficients

real(dp) :: T
real(dp) :: desdTdp

! ----

if (temp <= tfreeze) then
  b = bi
else
  b = bl
end if

T = temp - tfreeze  ! these coefficients are for temperature values in Celcius

desdTdp = sum(b * T**p)
  
desdT = real(100._dp * desdTdp)

end function desdT

! ----------------------------------------------------------------------------------------------------------------

end module physicsmod
