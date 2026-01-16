module constants

real(8), parameter :: PI = 4.d0*datan(1.d0)
real(8), parameter :: PIover2 = 0.5d0*PI
real(8), parameter :: twoPI = 2.d0*PI

real(8), parameter :: deg2rad = PI/180.d0
real(8), parameter :: rad2deg = 180.d0/PI

complex(8), parameter :: czero = dcmplx(0.d0, 0.d0)
complex(8), parameter :: cone  = dcmplx(0.d0, 1.d0)


end module constants
