!submain.for designed by liuyou-shan 2010-3-18
!=================================================================!
module dirac_coeff

! define the coefficients of the discretized Dirac-delta function

! abs(x) <= h
real(8), parameter :: C1_0 =   1.d0
real(8), parameter :: C1_1 =   0.d0
real(8), parameter :: C1_2 =  -5.d0 /  4.d0
real(8), parameter :: C1_3 = -35.d0 / 12.d0
real(8), parameter :: C1_4 =  21.d0 /  4.d0
real(8), parameter :: C1_5 = -25.d0 / 12.d0

! h < abs(x) <= 2*h
real(8), parameter :: C2_0 =   -4.d0
real(8), parameter :: C2_1 =   75.d0 /  4.d0
real(8), parameter :: C2_2 = -245.d0 /  8.d0
real(8), parameter :: C2_3 =  545.d0 / 24.d0
real(8), parameter :: C2_4 =  -63.d0 /  8.d0
real(8), parameter :: C2_5 =   25.d0 / 24.d0

! 2*h < abs(x) < 3*h
real(8), parameter :: C3_0 =   18.d0
real(8), parameter :: C3_1 = -153.d0 /  4.d0
real(8), parameter :: C3_2 =  255.d0 /  8.d0
real(8), parameter :: C3_3 = -313.d0 / 24.d0
real(8), parameter :: C3_4 =   21.d0 /  8.d0
real(8), parameter :: C3_5 =   -5.d0 / 24.d0

end module dirac_coeff
!=================================================================!
subroutine loadseis(nx, nz, npml, izt, x0, z0, dx, dz, rx, rz, ux, uz, displx, displz)

use dirac_coeff

implicit none

integer(4), intent(in) :: nx, nz
integer(4), intent(in) :: npml, izt

real(4), intent(in) :: rx, rz
real(4), intent(in) :: x0, z0, dx, dz

real(4), intent(in) :: ux(1:2, izt:nz+npml, 1-npml:nx+npml)
real(4), intent(in) :: uz(1:2, izt:nz+npml, 1-npml:nx+npml)

real(4), intent(out) :: displx, displz


integer(4) ix, iz, irx, irz

real(8) ba
real(8) x, z, bax, baz
real(8) x2, x3, z2, z3



irx = floor((rx-x0)/dx) + 1
irz = floor((rz-z0)/dz) + 1


displx = 0.0
displz = 0.0
!$omp parallel default(shared) private(ix, iz, x, x2, x3, z, z2, z3, ba, bax, baz) reduction(+: displx, displz)
!$omp do schedule(dynamic)
do ix = max(1-npml,irx-3), min(nx+npml,irx+3), 1

   x = dabs( dble( dble(x0 + (ix-1)*dx) - rx ) ) / dble(dx)
   x2 = x*x
   x3 = x2*x

   if(x <= 1.0)then
      !bax = (1.d0 - x2*(5.d0-21.d0*x2)/4.d0 - x3*(35.d0+25.d0*x2)/12.d0)
      bax = (C1_0 + x2*(C1_2 + x2*C1_4) + x3*(C1_3 + x2*C1_5))
   else if(x <= 2.0)then
      !bax = (-4.d0 + 75.d0/4.d0*x - x2*(245.d0+63.d0*x2)/8.d0 + x3*(545.d0+25.d0*x2)/24.d0)
      bax = (C2_0 + C2_1*x + x2*(C2_2 + x2*C2_4) + x3*(C2_3 + x2*C2_5))
   else if(x <= 3.0)then
      !bax = (18.d0 - 153.d0/4.d0*x + x2*(255.d0+21.d0*x2)/8.d0 - x3*(313.d0+5.d0*x2)/24.d0)
      bax = (C3_0 + C3_1*x + x2*(C3_2 + x2*C3_4) + x3*(C3_3 + x2*C3_5))
   else
      bax = 0.0
   end if

   do iz = max(izt,irz-3), min(nz+npml,irz+3), 1

      z = dabs( dble( dble(z0 + (iz-1)*dz) - rz ) ) / dble(dz)
      z2 = z*z
      z3 = z2*z

      if(z <= 1.0)then
         !baz = (1.d0 - z2*(5.d0-21.d0*z2)/4.d0 - z3*(35.d0+25.d0*z2)/12.d0)
         baz = (C1_0 + z2*(C1_2 + z2*C1_4) + z3*(C1_3 + z2*C1_5))
      else if(z <= 2.0)then
         !baz = (-4.d0 + 75.d0/4.d0*z - z2*(245.d0+63.d0*z2)/8.d0 + z3*(545.d0+25.d0*z2)/24.d0)
         baz = (C2_0 + C2_1*z + z2*(C2_2 + z2*C2_4) + z3*(C2_3 + z2*C2_5))
      else if(z <= 3.0)then
         !baz = (18.d0 - 153.d0/4.d0*z + z2*(255.d0+21.d0*z2)/8.d0 - z3*(313.d0+5.d0*z2)/24.d0)
         baz = (C3_0 + C3_1*z + z2*(C3_2 + z2*C3_4) + z3*(C3_3 + z2*C3_5))
      else
         baz = 0.0
      end if

      ba = bax*baz
      displx = displx + ba*ux(1,iz,ix)
      displz = displz + ba*uz(1,iz,ix)

   end do

end do
!$omp end do
!$omp end parallel


end subroutine loadseis
!=================================================================!
subroutine check_stability(nx, nz, npml, izt, norder, dx, dz, dt, coef2, c)

implicit none

integer(4), intent(in) :: nx, nz
integer(4), intent(in) :: npml, izt
integer(4), intent(in) :: norder

real(8), intent(in) :: dt

real(4), intent(in) :: dx, dz

real(8), intent(in) :: coef2(1:norder, 1:norder)

real, intent(in) :: c(1:3, izt:nz+npml, 1-npml:nx+npml)


real CFL, dt_stable
real Vpmax, h, cx, cz


h = max(dx, dz)
cx = (h/dx)*(h/dx)
cz = (h/dz)*(h/dz)
CFL = sqrt(2.0 / ((cx + cz)*(sum(abs(coef2(1:norder,norder))) + sum(coef2(1:norder,norder)))))
Vpmax = maxval(sqrt((c(2,:,:) + 2.0*c(3,:,:)) / c(1,:,:)))


dt_stable = CFL * h / Vpmax


if ( dt > dt_stable ) then
   write(*,*)
   write(*,"(A)") 'ERROR: The dt is too large to stablize the simulation !'
   write(*,"(A, F15.7, A)") 'The stable dt = ', dt_stable, ' is recommended'
   write(*,*)
   stop
end if

if ( dt < 0.4*dt_stable ) then
   write(*,*)
   write(*,"(A)") 'WARNING: The dt is too small for an efficient modeling!'
   write(*,"(A, F15.7, A)") 'The stable dt = ', dt_stable, ' is recommended'
   write(*,*)
end if
   
end subroutine check_stability
!=================================================================!
!dire denotes the direction of gradient, 1 for x direction; 2 for z direction
subroutine diff1(norder, dire, nx, nz, npml, izt, ix, iz, d, coef, u, tmp)

implicit none

integer, intent(in) :: ix, iz, izt
integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, dire

real, intent(in) :: d

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: tmp


integer k, kx, kz
integer order, node
integer norderx, norderz
integer lx1, lx2, lz1, lz2

real mysum


norderx = min(ix+npml-1, nx+npml-ix)
norderx = min(norderx, norder)
norderz = min(iz-izt, nz+npml-iz)
norderz = min(norderz, norder)

kx = 2-dire
kz = dire-1
node = kx*norderx + kz*norderz

mysum = 0.0

do k = 1, node, 1
   lx2 = ix+kx*k; lx1 = ix-kx*k
   lz2 = iz+kz*k; lz1 = iz-kz*k
   mysum = mysum + coef(k,node)*(u(lz2,lx2) - u(lz1,lx1))
end do

tmp = mysum / d


end subroutine diff1
!=================================================================!
!dire denotes the direction of gradient, 1 for x direction; 2 for z direction
subroutine diffx1(norder, nx, npml, ix, dx, coef, u, tmp)

implicit none

integer, intent(in) :: nx, npml
integer, intent(in) :: norder, ix

real, intent(in) :: dx

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(1-npml:nx+npml)


real, intent(out) :: tmp


integer k, norderx

real mysum


norderx = min(ix+npml-1, nx+npml-ix)
norderx = min(norderx, norder)

mysum = 0.0

do k = 1, norderx, 1
   mysum = mysum + coef(k,norderx)*(u(ix+k) - u(ix-k))
end do

tmp = mysum / dx


end subroutine diffx1
!=================================================================!
subroutine diffzx_reg(norder, nx, nz, npml, izt, dx, dz, coef, mat, u, uzx)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx, dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)

real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)


real, intent(out) :: uzx(izt:nz+npml, 1-npml:nx+npml)


integer ix, iz, k
integer norderx, norderz


real mysum

integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v



!! slow
!allocate(v(izt:nz+npml, 1-npml:nx+npml))
!allocate(mask(izt:nz+npml, 1-npml:nx+npml))
!
!
!! first-order derivative in z direction
!mask = 0
!mask(2:nz,1-norder:nx+norder) = 1
!
!
!v = 0.0
!!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!!$omp do schedule(dynamic)
!do iz = 2, nz, 1
!
!   norderz = min(iz-izt, nz+npml-iz)
!   norderz = min(norderz, norder)
!
!   do ix = 1-norder, nx+norder, 1
!
!      if(0 == mask(iz,ix)) cycle
!
!      mysum  = 0.0
!      do k = 1, norderz, 1
!         mysum = mysum + coef(k,norderz) * (u(iz+k,ix) - u(iz-k,ix))
!      end do
!      v(iz,ix) = mat(iz,ix)*mysum / dz
!
!   end do
!
!end do
!!$omp end do
!!$omp end parallel
!
!
!! first-order derivative in x direction
!mask = 0
!mask(2:nz,1:nx) = 1
!
!
!uzx = 0.0
!!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!!$omp do schedule(dynamic)
!do ix = 1, nx, 1
!
!   norderx = min(ix+npml-1, nx+npml-ix)
!   norderx = min(norderx, norder)
!
!   do iz = 2, nz, 1
!
!      if(0 == mask(iz,ix)) cycle
!
!      mysum = 0.0
!      do k = 1, norderx, 1
!         mysum = mysum + coef(k,norderx) * (v(iz,ix+k) - v(iz,ix-k))
!      end do
!      uzx(iz,ix) = mysum / dx
!
!   end do
!
!end do
!!$omp end do
!!$omp end parallel
!
!
!deallocate(v, mask)



! fast
allocate(v(izt:nz+npml, 1-npml:nx+npml))


! first-order derivative in z direction

v = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!$omp do schedule(dynamic)
do iz = 2, nz, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 1-norder, nx+norder, 1

      mysum  = 0.0
      do k = 1, norderz, 1
         mysum = mysum + coef(k,norderz) * (u(iz+k,ix) - u(iz-k,ix))
      end do
      v(iz,ix) = mat(iz,ix)*mysum / dz

   end do

end do
!$omp end do
!$omp end parallel


! first-order derivative in x direction


uzx = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!$omp do schedule(dynamic)
do ix = 1, nx, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = 2, nz, 1

      mysum = 0.0
      do k = 1, norderx, 1
         mysum = mysum + coef(k,norderx) * (v(iz,ix+k) - v(iz,ix-k))
      end do
      uzx(iz,ix) = mysum / dx

   end do

end do
!$omp end do
!$omp end parallel


deallocate(v)


end subroutine diffzx_reg
!=================================================================!
subroutine diffzx_pml(norder, nx, nz, npml, izt, dx, dz, coef, mat, u, uzx)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx, dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: uzx(izt:nz+npml, 1-npml:nx+npml)


integer ix, iz, k
integer norderx, norderz

real mysum

integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v



allocate(v(izt:nz+npml, 1-npml:nx+npml))
allocate(mask(izt:nz+npml, 1-npml:nx+npml))


! first-order derivative in z direction
mask = 1
mask(1:nz,norder+1:nx-norder) = 0
mask(izt,:) = 0


v = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!$omp do schedule(dynamic)
do iz = izt+1, nz+npml-1, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 2-npml, nx+npml-1, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderz, 1
         mysum = mysum + coef(k,norderz) * (u(iz+k,ix) - u(iz-k,ix))
      end do
      v(iz,ix) = mat(iz,ix)*mysum / dz

   end do

end do
!$omp end do
!$omp end parallel


! first-order derivative in x direction
mask = 1
mask(1:nz,1:nx) = 0
mask(izt,:) = 0

uzx = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!$omp do schedule(dynamic)
do ix = 2-npml, nx+npml-1, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = izt+1, nz+npml-1, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderx, 1
         mysum = mysum + coef(k,norderx) * (v(iz,ix+k) - v(iz,ix-k))
      end do
      uzx(iz,ix) = mysum / dx

   end do

end do
!$omp end do
!$omp end parallel


deallocate(v, mask)


end subroutine diffzx_pml
!=================================================================!
subroutine diffxz_reg(norder, nx, nz, npml, izt, dx, dz, coef, mat, u, uxz)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx, dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: uxz(izt:nz+npml, 1-npml:nx+npml)


integer ix, iz, k
integer norderx, norderz


real mysum

integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v



!! slow
!allocate(v(izt:nz+npml, 1-npml:nx+npml))
!allocate(mask(izt:nz+npml, 1-npml:nx+npml))
!
!
!! first-order derivative in x direction
!mask = 0
!mask(1:nz+norder,1:nx) = 1
!
!
!v = 0.0
!!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!!$omp do schedule(dynamic)
!do ix = 1, nx, 1
!
!   norderx = min(ix+npml-1, nx+npml-ix)
!   norderx = min(norderx, norder)
!
!   do iz = 1, nz+norder, 1
!
!      if(0 == mask(iz,ix)) cycle
!
!      mysum = 0.0
!      do k = 1, norderx, 1
!         mysum = mysum + coef(k,norderx) * (u(iz,ix+k) - u(iz,ix-k))
!      end do
!      v(iz,ix) = mat(iz,ix)*mysum / dx
!
!   end do
!
!end do
!!$omp end do
!!$omp end parallel
!
!
!! first-order derivative in z direction
!mask = 0
!mask(2:nz,1:nx) = 1
!
!
!uxz = 0.0
!!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!!$omp do schedule(dynamic)
!do iz = 2, nz, 1
!
!   norderz = min(iz-izt, nz+npml-iz)
!   norderz = min(norderz, norder)
!
!   do ix = 1, nx, 1
!
!      if(0 == mask(iz,ix)) cycle
!
!      mysum = 0.0
!      do k = 1, norderz, 1
!         mysum = mysum + coef(k,norderz) * (v(iz+k,ix) - v(iz-k,ix))
!      end do
!      uxz(iz,ix) = mysum / dz
!
!   end do
!
!end do
!!$omp end do
!!$omp end parallel
!
!
!deallocate(v, mask)



! slow
allocate(v(izt:nz+npml, 1-npml:nx+npml))


! first-order derivative in x direction

v = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!$omp do schedule(dynamic)
do ix = 1, nx, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = 1, nz+norder, 1

      mysum = 0.0
      do k = 1, norderx, 1
         mysum = mysum + coef(k,norderx) * (u(iz,ix+k) - u(iz,ix-k))
      end do
      v(iz,ix) = mat(iz,ix)*mysum / dx

   end do

end do
!$omp end do
!$omp end parallel


! first-order derivative in z direction

uxz = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!$omp do schedule(dynamic)
do iz = 2, nz, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 1, nx, 1

      mysum = 0.0
      do k = 1, norderz, 1
         mysum = mysum + coef(k,norderz) * (v(iz+k,ix) - v(iz-k,ix))
      end do
      uxz(iz,ix) = mysum / dz

   end do

end do
!$omp end do
!$omp end parallel


deallocate(v)


end subroutine diffxz_reg
!=================================================================!
subroutine diffxz_pml(norder, nx, nz, npml, izt, dx, dz, coef, mat, u, uxz)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx, dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: uxz(izt:nz+npml, 1-npml:nx+npml)


integer ix, iz, k
integer norderx, norderz

real mysum

integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v



allocate(v(izt:nz+npml, 1-npml:nx+npml))
allocate(mask(izt:nz+npml, 1-npml:nx+npml))


! derivative in x direction
mask = 1
mask(1:nz-norder,1:nx) = 0


v = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderx, mysum)
!$omp do schedule(dynamic)
do ix = 2-npml, nx+npml-1, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = izt, nz+npml-1, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderx, 1
         mysum = mysum + coef(k,norderx) * (u(iz,ix+k) - u(iz,ix-k))
      end do
      v(iz,ix) = mat(iz,ix)*mysum / dx

   end do

end do
!$omp end do
!$omp end parallel


! first-order derivative in z direction
mask = 1
mask(1:nz,1:nx) = 0
mask(izt,:) = 0


uxz = 0.0
!$omp parallel default(shared) private(ix, iz, k, norderz, mysum)
!$omp do schedule(dynamic)
do iz = izt+1, nz+npml-1, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 2-npml, nx+npml-1, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderz, 1
         mysum = mysum + coef(k,norderz) * (v(iz+k,ix) - v(iz-k,ix))
      end do
      uxz(iz,ix) = mysum / dz

   end do

end do
!$omp end do
!$omp end parallel


deallocate(v, mask)


end subroutine diffxz_pml
!=================================================================!
subroutine diffx2_reg(norder, nx, nz, npml, izt, dx, coef, mat, u, ux2)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: ux2(izt:nz+npml, 1-npml:nx+npml)


integer i, j, k, ix, iz
integer norderx, norderx1, norderx2


real dx2, mysum
real mysum1, mysum2

integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v


allocate(mask(izt:nz+npml, 1-npml:nx+npml))


! second-order derivative in x direction
mask = 0
mask(1:nz,1:nx) = 1


ux2 = 0.d0
dx2 = 1.d0/(dx*dx)
!$omp parallel default(shared) private(ix, iz, i, j, k, norderx, norderx1, norderx2, mysum, mysum1, mysum2)
!$omp do schedule(dynamic)
do ix = 1, nx, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = izt, nz, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderx, 1

         i = int((2*k-1)/2)

         mysum1 = 0.0
         norderx1 = min(ix+i+npml-1, nx+npml-ix-i)
         norderx1 = min(norderx1, norder)
         do j = 1, norderx1, 1
            mysum1 = mysum1 + coef(j,norderx1)*(u(iz,ix+j+k-1) - u(iz,ix-j+k))
         end do
         
         mysum2 = 0.0
         norderx2 = min(ix-i+npml-1, nx+npml-ix+i)
         norderx2 = min(norderx2, norder)
         do j = 1, norderx2, 1
            mysum2 = mysum2 + coef(j,norderx2)*(u(iz,ix+j-k) - u(iz,ix-j-k+1))
         end do

         mysum = mysum + coef(k,norderx)*(2.0/(1.0/mat(iz,ix+k)+1.0/mat(iz,ix+k-1))*mysum1 - &
                                          2.0/(1.0/mat(iz,ix-k+1)+1.0/mat(iz,ix-k))*mysum2)

      end do
      ux2(iz,ix) = mysum * dx2

   end do

end do
!$omp end do
!$omp end parallel


deallocate(mask)


end subroutine diffx2_reg
!=================================================================!
subroutine diffx2_pml(norder, nx, nz, npml, izt, dx, coef, mat, u, ux2)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dx

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: ux2(izt:nz+npml, 1-npml:nx+npml)


integer i, j, k, ix, iz
integer norderx, norderx1, norderx2


real dx2, mysum
real mysum1, mysum2


integer(1), allocatable, dimension(:,:) :: mask

real, allocatable, dimension(:,:) :: v, ux3


allocate(mask(izt:nz+npml, 1-npml:nx+npml))


! second-order derivative in x direction
mask = 1
mask(1:nz,1:nx) = 0


ux2 = 0.0
dx2 = 1.0/(dx*dx)
!$omp parallel default(shared) private(ix, iz, i, j, k, norderx, norderx1, norderx2, mysum, mysum1, mysum2)
!$omp do schedule(dynamic)
do ix = 2-npml, nx+npml-1, 1

   norderx = min(ix+npml-1, nx+npml-ix)
   norderx = min(norderx, norder)

   do iz = izt, nz+npml-1, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderx, 1

         i = int((2*k-1)/2)

         mysum1 = 0.0
         norderx1 = min(ix+i+npml-1, nx+npml-ix-i)
         norderx1 = min(norderx1, norder)
         do j = 1, norderx1, 1
            mysum1 = mysum1 + coef(j,norderx1)*(u(iz,ix+j+k-1) - u(iz,ix-j+k))
         end do
         
         mysum2 = 0.0
         norderx2 = min(ix-i+npml-1, nx+npml-ix+i)
         norderx2 = min(norderx2, norder)
         do j = 1, norderx2, 1
            mysum2 = mysum2 + coef(j,norderx2)*(u(iz,ix+j-k) - u(iz,ix-j-k+1))
         end do

         mysum = mysum + coef(k,norderx)*(2.0/(1.0/mat(iz,ix+k)+1.0/mat(iz,ix+k-1))*mysum1 - &
                                          2.0/(1.0/mat(iz,ix-k+1)+1.0/mat(iz,ix-k))*mysum2)

      end do
      ux2(iz,ix) = mysum * dx2

   end do

end do
!$omp end do
!$omp end parallel


deallocate(mask)


end subroutine diffx2_pml
!=================================================================!
subroutine diffz2_reg(norder, nx, nz, npml, izt, dz, coef, mat, u, uz2)

implicit none

integer, intent(in) :: norder, izt
integer, intent(in) :: nx, nz, npml

real, intent(in) :: dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: uz2(izt:nz+npml, 1-npml:nx+npml)


integer i, j, k, ix, iz
integer norderz, norderz1, norderz2


real dz2, mysum
real mysum1, mysum2


integer(1), allocatable, dimension(:,:) :: mask


allocate(mask(izt:nz+npml, 1-npml:nx+npml))


! second-order derivative in z direction
mask = 0
mask(2:nz,1:nx) = 1


uz2 = 0.0
dz2 = 1.0/(dz*dz)
!$omp parallel default(shared) private(ix, iz, i, j, k, norderz, norderz1, norderz2, mysum, mysum1, mysum2)
!$omp do schedule(dynamic)
do iz = 2, nz, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 1, nx, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderz, 1

         i = int((2*k-1)/2)

         mysum1 = 0.0
         norderz1 = min(iz+i-izt, nz+npml-iz-i)
         norderz1 = min(norderz1, norder)
         do j = 1, norderz1, 1
            mysum1 = mysum1 + coef(j,norderz1)*(u(iz+j+k-1,ix) - u(iz-j+k,ix))
         end do

         mysum2 = 0.0
         norderz2 = min(iz-i-izt, nz+npml-iz+i)
         norderz2 = min(norderz2, norder)
         do j = 1, norderz2, 1
            mysum2 = mysum2 + coef(j,norderz2)*(u(iz+j-k,ix) - u(iz-j-k+1,ix))
         end do

         mysum = mysum + coef(k,norderz)*(2.0/(1.0/mat(iz+k,ix)+1.0/mat(iz+k-1,ix))*mysum1 - &
                                          2.0/(1.0/mat(iz-k+1,ix)+1.0/mat(iz-k,ix))*mysum2)

      end do
      uz2(iz,ix) = mysum * dz2

   end do

end do
!$omp end do
!$omp end parallel


deallocate(mask)


end subroutine diffz2_reg
!=================================================================!
subroutine diffz2_pml(norder, nx, nz, npml, izt, dz, coef, mat, u, uz2)

implicit none

integer, intent(in) :: nx, nz, npml
integer, intent(in) :: norder, izt

real, intent(in) :: dz

real(8), intent(in) :: coef(1:norder, 1:norder)

real, intent(in) :: u(izt:nz+npml, 1-npml:nx+npml)
real, intent(in) :: mat(izt:nz+npml, 1-npml:nx+npml)

real, intent(out) :: uz2(izt:nz+npml, 1-npml:nx+npml)


integer i, j, k, ix, iz
integer norderz, norderz1, norderz2


real dz2, mysum
real mysum1, mysum2

integer(1), allocatable, dimension(:,:) :: mask


allocate(mask(izt:nz+npml, 1-npml:nx+npml))



! derivative in z direction
mask = 1
mask(1:nz,1:nx) = 0


uz2 = 0.0
dz2 = 1.0/(dz*dz)
!$omp parallel default(shared) private(ix, iz, i, j, k, norderz, norderz1, norderz2, mysum, mysum1, mysum2)
!$omp do schedule(dynamic)
do iz = izt+1, nz+npml-1, 1

   norderz = min(iz-izt, nz+npml-iz)
   norderz = min(norderz, norder)

   do ix = 1-npml, nx+npml, 1

      if(0 == mask(iz,ix)) cycle

      mysum = 0.0
      do k = 1, norderz, 1

         i = int((2*k-1)/2)

         mysum1 = 0.0
         norderz1 = min(iz+i-izt, nz+npml-iz-i)
         norderz1 = min(norderz1, norder)
         do j = 1, norderz1, 1
            mysum1 = mysum1 + coef(j,norderz1)*(u(iz+j+k-1,ix) - u(iz-j+k,ix))
         end do

         mysum2 = 0.0
         norderz2 = min(iz-i-izt, nz+npml-iz+i)
         norderz2 = min(norderz2, norder)
         do j = 1, norderz2, 1
            mysum2 = mysum2 + coef(j,norderz2)*(u(iz+j-k,ix) - u(iz-j-k+1,ix))
         end do

         mysum = mysum + coef(k,norderz)*(2.0/(1.0/mat(iz+k,ix)+1.0/mat(iz+k-1,ix))*mysum1 - &
                                          2.0/(1.0/mat(iz-k+1,ix)+1.0/mat(iz-k,ix))*mysum2)

      end do
      uz2(iz,ix) = mysum * dz2

   end do

end do
!$omp end do
!$omp end parallel


deallocate(mask)


end subroutine diffz2_pml
!=====================================The parameter of Dire denotes the times of gradient==============================!
subroutine fd_coef0(norder, coef)

! first derivative of FD cofficients in staggered grids

implicit none

integer(4), intent(in) :: norder

real(8), intent(out) :: coef(1:norder, 1:norder)

integer(4) i, j, k

real(8) tmp1, tmp2, tmp3


coef = 0.d0

! second-order derivative (regular grid)
do i = 1, norder, 1

   do j = 1, i, 1

      tmp1 = 1.d0
      tmp2 = 1.d0
      tmp3 = 1.d0

      do k = 1, j-1, 1
         tmp1 = tmp1*dble((2*k-1)*(2*k-1))
         tmp2 = tmp2*(dble((2*j-1)*(2*j-1))-dble((2*k-1)*(2*k-1)))
      end do
      do k = j+1, i, 1
         tmp1 = tmp1*dble((2*k-1)*(2*k-1))
         tmp3 = tmp3*(dble((2*k-1)*(2*k-1))-dble((2*j-1)*(2*j-1)))
      end do
      tmp1 = tmp1*(-1)**(j+1)
      tmp2 = tmp2*dble(2*j-1)

      coef(j,i) = tmp1/tmp2/tmp3

   end do

end do


end subroutine fd_coef0
!=================================================================!
subroutine fd_coef1(norder, coef)

! first derivative of FD cofficients in collocated grids

implicit none

integer(4), intent(in) :: norder

real(8), intent(out) :: coef(1:norder, 1:norder)

integer(4) i, j, k

real(8) tmp1, tmp2, tmp3


coef = 0.d0

do i = 1, norder, 1

   do j = 1, i, 1

      tmp1 = 1.d0
      tmp2 = 1.d0
      tmp3 = 1.d0

      do k = 1, j-1, 1
         tmp1 = tmp1*dble(k*k)
         tmp2 = tmp2*(dble(j*j)-dble(k*k))
      end do
      do k = j+1, i, 1
         tmp1 = tmp1*dble(k*k)
         tmp3 = tmp3*(dble(k*k)-dble(j*j))
      end do
      tmp1 = tmp1*(-1)**(j+1)
      tmp2 = tmp2*dble(2.d0*j)

      coef(j,i) = tmp1/tmp2/tmp3

   end do

end do


end subroutine fd_coef1
!=====================================The parameter of Dire denotes the times of gradient==============================!
subroutine fd_coef2(norder, coef)

! second derivative of FD cofficients in collocated grids

implicit none

integer(4), intent(in) :: norder

real(8), intent(out) :: coef(1:norder, 1:norder)

integer(4) i, j, k

real(8) tmp1, tmp2, tmp3


coef = 0.d0

do i = 1, norder, 1

   do j = 1, i, 1

      tmp1 = 1.d0
      tmp2 = 1.d0
      tmp3 = 1.d0

      do k = 1, j-1, 1
         tmp1 = tmp1*dble(k*k)
         tmp2 = tmp2*(dble(j*j)-dble(k*k))
      end do
      do k = j+1, i, 1
         tmp1 = tmp1*dble(k*k)
         tmp3 = tmp3*(dble(k*k)-dble(j*j))
      end do
      tmp1 = tmp1*(-1)**(j+1)
      tmp2 = tmp2*dble(j*j)

      coef(j,i) = tmp1/tmp2/tmp3

   end do

end do


end subroutine fd_coef2
!=================================================================!
subroutine zoeppritz(Vp1, Vs1, rho1, alpha1, beta1, Vp2, Vs2, rho2, alpha2, beta2, Rpp, Rps, Tpp, Tps)

use constants, only: deg2rad

implicit none

real, intent(in) :: Vp1, Vp2
real, intent(in) :: Vs1, Vs2
real, intent(in) :: rho1, rho2

real(8), intent(in) :: beta1, beta2
real(8), intent(in) :: alpha1, alpha2


real(8), intent(out) :: Rpp, Rps, Tpp, Tps


real(8) rho_ratio

real(8) A(1:4, 1:4)
real(8) B(1:4, 1:1)
real(8) E(1:4, 1:4)
real(8) X(1:4, 1:1)


!write(*,*) rho1, Vp1, Vs1
!write(*,*) rho2, Vp2, Vs2

rho_ratio = rho2/rho1

A(1,1) =  sin(alpha1*deg2rad)
A(1,2) =  cos(beta1*deg2rad)
A(1,3) = -sin(alpha2*deg2rad)
A(1,4) =  cos(beta2*deg2rad)

A(2,1) =  cos(alpha1*deg2rad)
A(2,2) = -sin(beta1*deg2rad)
A(2,3) =  cos(alpha2*deg2rad)
A(2,4) =  sin(beta2*deg2rad)

A(3,1) =  cos(2.d0*beta1*deg2rad)
A(3,2) = -sin(2.d0*beta1*deg2rad)*Vs1/Vp1
A(3,3) = -cos(2.d0*beta2*deg2rad)*rho_ratio*Vp2/Vp1
A(3,4) = -sin(2.d0*beta2*deg2rad)*rho_ratio*Vs2/Vp1

A(4,1) =  sin(2.d0*alpha1*deg2rad)*Vs1*Vs1/Vp1
A(4,2) =  cos(2.d0*beta1*deg2rad)*Vs1
A(4,3) =  sin(2.d0*alpha2*deg2rad)*rho_ratio*Vs2*Vs2/Vp2
A(4,4) = -cos(2.d0*beta2*deg2rad)*rho_ratio*Vs2

call inverse_martrix(4, A, E)

B(1,1) = -sin(alpha1*deg2rad)
B(2,1) =  cos(alpha1*deg2rad)
B(3,1) = -cos(2.d0*beta1*deg2rad)
B(4,1) =  sin(2.d0*alpha1*deg2rad)*Vs1*Vs1/Vp1

X = matmul(E, B)

Rpp = X(1,1)
Rps = X(2,1)
Tpp = X(3,1)
Tps = X(4,1)

end subroutine zoeppritz
!=================================================================!
subroutine inverse_martrix(n, a, b)

implicit none

integer(4), intent(in) :: n

real(8), intent(in) :: a(1:n,1:n)
real(8), intent(out) :: b(1:n,1:n)


integer(4) info
integer(4) ipiv(1:n)
 
real(8) work(1:n)
real(8) ac(1:n,1:n)


ac = a
call dgetrf(n, n, ac, n, ipiv, info)
call dgetri(n, ac, n, ipiv, work, n, info)
b = ac


end subroutine inverse_martrix
!=================================================================!
