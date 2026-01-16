subroutine define_rotation_matrix(lat_org, lon_org, az_org, rotate_matrix, inv_rotate_matrix)

use constants, only: PIover2, deg2rad

implicit none

real(8), intent(in) :: lat_org, lon_org, az_org

real(8), intent(out) :: rotate_matrix(3,3), inv_rotate_matrix(3,3)


real(8) lat, lon, az, theta
real(8) sin_az, cos_az, factor1, factor2
real(8) sin_theta, cos_theta, sin_phi, cos_phi

!real(8), dimension(3,3) :: R, R1, R2, R3


az  = dble(az_org)*deg2rad
lat = dble(lat_org)*deg2rad
lon = dble(lon_org)*deg2rad


!! Spherical cartesian coordinate to local cartesian coordinate at station
!! (x, y, z) -> (x', y', z')
!theta = PIover2 + lon
!R1(1,1) =  cos(theta); R1(1,2) = sin(theta); R1(1,3) = 0.d0
!R1(2,1) = -sin(theta); R1(2,2) = cos(theta); R1(2,3) = 0.d0
!R1(3,1) = 0.d0; R1(3,2) = 0.d0; R1(3,3) = 1.d0
!
!
!! (x', y', z') -> (T, R, Z)
!theta = PIover2 - lat
!R2(1,1) = 1.d0; R2(1,2) = 0.d0; R2(1,3) = 0.d0
!R2(2,1) = 0.d0; R2(2,2) =  cos(theta); R2(2,3) = sin(theta)
!R2(3,1) = 0.d0; R2(3,2) = -sin(theta); R2(3,3) = cos(theta)
!
!
!! (T, R, Z) -> (X, Y, Z)
!theta = PIover2 - az
!R3(1,1) = cos(theta); R3(1,2) = -sin(theta); R3(1,3) = 0.d0
!R3(2,1) = sin(theta); R3(2,2) =  cos(theta); R3(2,3) = 0.d0
!R3(3,1) = 0.d0; R3(3,2) = 0.d0; R3(3,3) = 1.d0
!
!
!
!! rotation matrix
!R = matmul(R2, R1)
!rotate_matrix = matmul(R3, R)
!
!
!! inverse rotation matrix
!R = matmul(R2, R3)
!R = matmul(R1, R)
!R(1,2) = -R(1,2)
!R(2,1) = -R(2,1)
!R(2,3) = -R(2,3)
!R(3,2) = -R(3,2)
!inv_rotate_matrix = R


! or
sin_theta = sin(lat); cos_theta = cos(lat)
sin_phi   = sin(lon); cos_phi   = cos(lon)
sin_az    = sin(az) ; cos_az    = cos(az)
factor1 = sin_theta*cos_phi
factor2 = sin_theta*sin_phi
rotate_matrix(1,1) = -sin_az*sin_phi + cos_az*factor1
rotate_matrix(1,2) =  sin_az*cos_phi + cos_az*factor2
rotate_matrix(1,3) = -cos_az*cos_theta
rotate_matrix(2,1) = -cos_az*sin_phi - sin_az*factor1
rotate_matrix(2,2) =  cos_az*cos_phi - sin_az*factor2
rotate_matrix(2,3) =  sin_az*cos_theta
rotate_matrix(3,1) =  cos_theta*cos_phi
rotate_matrix(3,2) =  cos_theta*sin_phi
rotate_matrix(3,3) =  sin_theta

inv_rotate_matrix(1,1) = -sin_az*sin_phi + cos_az*factor1
inv_rotate_matrix(1,2) = -cos_az*sin_phi - sin_az*factor2
inv_rotate_matrix(1,3) = -cos_theta*cos_phi
inv_rotate_matrix(2,1) =  sin_az*cos_phi + cos_az*factor1
inv_rotate_matrix(2,2) =  cos_az*cos_phi - sin_az*factor2
inv_rotate_matrix(2,3) =  cos_theta*cos_phi
inv_rotate_matrix(3,1) = -cos_az*cos_theta
inv_rotate_matrix(3,2) =  sin_az*cos_phi
inv_rotate_matrix(3,3) =  sin_theta


end subroutine define_rotation_matrix

