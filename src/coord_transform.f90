subroutine compute_backazimuth_grids(nx, coordx, evla, evlo, bazs)

use calc_azimuth
use constants, only: deg2rad
use FDFK_par, only: lat_org, lon_org, az_org, &
                    rotate_matrix, inv_rotate_matrix

implicit none

integer, intent(in) :: nx

real(8), intent(in) :: evla, evlo

real, intent(in) :: coordx(1:nx)

real(8), intent(out) :: bazs(1:nx)


integer ix

real(8) x, y, z
real(8) xo, yo, zo

real(8) stla, stlo, baz

real(8) XYZ(3,1)
real(8) xyzr(3,1)



call define_rotation_matrix(lat_org, lon_org, az_org, rotate_matrix, inv_rotate_matrix)


! Geographic coordinates to cartesian coordinates at center of chunk
call geogr2cart(lat_org, lon_org, xo, yo, zo)



! grids (x, y, z) to (theta, phi, r)
do ix = 1, nx, 1

   XYZ(1,1) = coordx(ix)
   XYZ(2,1) = 0.d0
   XYZ(3,1) = 0.d0

   ! from local cartesian coordinates to spherical cartesian coordinates
   xyzr = matmul(inv_rotate_matrix, XYZ)

   x = xyzr(1,1) + xo
   y = xyzr(2,1) + yo
   z = xyzr(3,1) + zo
   ! from spherical cartesian coordinates to geographic coordinates
   call cart2geogr(x, y, z, stla, stlo)

   call backaz(evla, evlo, stla, stlo, baz)

   bazs(ix) = baz*deg2rad

end do


end subroutine compute_backazimuth_grids
