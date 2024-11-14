module calc_azimuth

use constants, only: PI, twoPI, PIover2, deg2rad, rad2deg

implicit none

real, parameter :: Earth_raduis = 6371.e3

real(8), parameter :: Earth_ellipsoid = 1.d0/298.257d0


contains


!subroutine backaz(evla, evlo, stla, stlo, az, baz)
!
!implicit none
!
!real(8), intent(in) :: evla, evlo, stla, stlo
!
!real(8), intent(out) :: az, baz
!
!
!real(8) coeff, x, a, b, C, cosc, sinc
!
!
!coeff = (1.d0 - Earth_ellipsoid)*(1.d0 - Earth_ellipsoid)
!a = PIover2 - atan(coeff*tan(stla*deg2rad))
!b = PIover2 - atan(coeff*tan(evla*deg2rad))
!C = (stlo - evlo)*deg2rad
!
!cosc = cos(a)*cos(b) + sin(a)*sin(b)*cos(C)
!sinc = sqrt(1.d0 - cosc*cosc)
!
!x = sin(a)*sin(C)/sinc
!if (abs(x-1.d0) < 2.5e-16) then
!   az = 90.d0*sign(1.d0,x)
!else
!   az = asin(x)*rad2deg
!end if
!x = -sin(b)*sin(C)/sinc
!if (abs(x-1.d0) < 2.5e-16) then
!   baz = 90.d0*sign(1.d0,x)
!else
!   baz = asin(x)*rad2deg
!end if
!
!
!if (((stlo - evlo) >= 0) .and. ((stla - evla) > 0)) then
!    ! first quadrant
!    baz = 180.d0 - baz
!else if (((stlo - evlo) < 0) .and. ((stla - evla) > 0)) then
!    ! second quadrant
!    az = 360.d0 + az
!    baz = 180.d0 - baz
!else if (((stlo - evlo) <= 0) .and. ((stla - evla) < 0)) then
!    ! third quadrant
!    az = 180.d0 - az
!    baz = -baz
!else if (((stlo - evlo) > 0) .and. ((stla - evla) < 0)) then
!    ! fourth quadrant
!    az = 180.d0 - az
!    baz = 360.d0 + baz
!end if
!
!
!end subroutine backaz




subroutine backaz(evla, evlo, stla, stlo, baz, az, delta)
!
! Subroutine to calculate the Great Circle Arc distance
!    between two sets of geographic coordinates
!
! Given: stla => Latitude of first point (+N, -S) in degrees
!	      stlo => Longitude of first point (+E, -W) in degrees
!	      evla => Latitude of second point
!	      evlo => Longitude of second point
!
! Returns:  delta => Great Circle Arc distance in degrees
!	         az    => Azimuth from pt. 1 to pt. 2 in degrees
!	         baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
!
! If you are calculating station-epicenter pairs, pt. 1 is the station
!
! Equations take from Bullen, pages 154, 155
!
! T. Owens, September 19, 1991
!           Sept. 25 -- fixed az and baz calculations
!           Dec. 2006, changed for fortran95
!           May, 2007 -- added predel to get around OSX acos round-off NaN issue
!
!

implicit none

real(8), intent(in) :: evla, evlo, stla, stlo

real(8), intent(out) :: baz
real(8), optional, intent(out) :: az, delta


real(8) scola, slon, ecola, elon
real(8) a, b, c, d, e, aa, bb, cc, dd, ee
real(8) g, gg, h, hh, k, kk, coeff, predel
real(8) rhs1, rhs2, sph, rad, del, daz, dbaz


!
! scola and ecola are the geocentric colatitudes
! as defined by Richter (pg. 318)
!
! Earth Flattening of 1/298.257 take from Bott (pg. 3)
!


coeff = (1.d0 - Earth_ellipsoid)*(1.d0 - Earth_ellipsoid)
scola = PIover2 - atan(coeff*tan(stla*deg2rad))
ecola = PIover2 - atan(coeff*tan(evla*deg2rad))
slon = stlo*deg2rad
elon = evlo*deg2rad
!
!  a - e are as defined by Bullen (pg. 154, Sec 10.2)
!     These are defined for the pt. 1
!
a = sin(scola)*cos(slon)
b = sin(scola)*sin(slon)
c = cos(scola)
d = sin(slon)
e = -cos(slon)
g = -c*e
h = c*d
k = -sin(scola)
!
!  aa - ee are the same as a - e, except for pt. 2
!
aa = sin(ecola)*cos(elon)
bb = sin(ecola)*sin(elon)
cc = cos(ecola)
dd = sin(elon)
ee = -cos(elon)
gg = -cc*ee
hh = cc*dd
kk = -sin(ecola)
!
!  Bullen, Sec 10.2, eqn. 4
!
predel = a*aa + b*bb + c*cc
if(abs(predel + 1.d0) < 1.d-6) then
   predel = -1.d0
endif
if(abs(predel - 1.d0) < 1.d-6) then
   predel = 1.d0
endif
del = acos(predel)
if (present(delta)) then
   delta = del*rad2deg
end if
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 1 is unprimed, so this is technically the baz
!
!  Calculate baz this way to avoid quadrant problems
!
rhs1 = (aa-d)*(aa-d) + (bb-e)*(bb-e) + cc*cc - 2.d0
rhs2 = (aa-g)*(aa-g) + (bb-h)*(bb-h) + (cc-k)*(cc-k) - 2.d0
dbaz = atan2(rhs1,rhs2)
if(dbaz < 0.d0) dbaz = dbaz + twoPI
baz = dbaz*rad2deg
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 2 is unprimed, so this is technically the az
!
if (present(az)) then
   rhs1 = (a-dd)*(a-dd) + (b-ee)*(b-ee) + c*c - 2.d0
   rhs2 = (a-gg)*(a-gg) + (b-hh)*(b-hh) + (c-kk)*(c-kk) - 2.d0
   daz = atan2(rhs1,rhs2)
   if(daz < 0.d0) daz = daz+ twoPI
   az = daz*rad2deg
end if
!
!   Make sure 0.0 is always 0.0, not 360.
!
if (abs(baz-360.d0) < 1.d-5) baz = 0.d0
if (present(az)) then
   if (abs(az-360.d0)  < 1.d-5) az = 0.d0
end if


return


end subroutine backaz



subroutine geogr2cart(lat, lon, x, y, z)
! (r, theta, phi) -> (x, y, z)

!use constants, only: Earth_raduis, deg2rad

implicit none

real(8), intent(in) :: lat, lon

real(8), intent(out) :: x, y, z


real(8) colatitude, longitude


colatitude = (90.d0 - lat)*deg2rad
longitude = lon*deg2rad

x = Earth_raduis*sin(colatitude)*cos(longitude)
y = Earth_raduis*sin(colatitude)*sin(longitude)
z = Earth_raduis*cos(colatitude)


end subroutine geogr2cart


subroutine cart2geogr(x, y, z, lat, lon)
! (x, y, z) -> (lat, lon)

!use constants, only: Earth_raduis, rad2deg

implicit none

real(8), intent(in) :: x, y, z

real(8), intent(out) :: lat, lon


real(8) r, tinyval


r = hypot(x, y)
tinyval = 1.d-6

if (r < tinyval) then
    if (z > 0.0) then
        lat = 90.d0
        lon = 0.d0
    else
        lat = -90.d0
        lon = 0.d0
    end if
    return
end if

lat = atan2(z, r)*rad2deg
lon = atan2(y, x)*rad2deg

if (abs(lat) < 1.e-9) then
   lat = 0.0
end if
if (abs(lon) < 1.e-9) then
   lon = 0.0
end if


end subroutine cart2geogr


end module calc_azimuth
