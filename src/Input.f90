!input.for designed by liuyou-shan 2010-3-18
!=================================================================!!
subroutine get_file_num(nfile)

use FDFK_par, only: prefix, finpar

implicit none

integer, intent(out) :: nfile

integer iunit, ioer


finpar = trim(prefix)//'/'//trim(finpar)

call allocate_unit(iunit)
open(unit=iunit, file=finpar, status='old')

   nfile = 0

   do

      read(iunit, *, iostat=ioer)
      if (0 /= ioer) exit
      nfile = nfile + 1

   end do

close(iunit)


end subroutine get_file_num
!!=================================================================!!
subroutine read_filename(nfile)

use FDFK_par, only: finpar, fFDmodel, fsource, &
                    freceiver, fFKmodel_left, &
                    fFKmodel_right, fVpmodel, fVsmodel

implicit none

integer, intent(in) :: nfile


integer i, iunit

character(256), dimension(:), allocatable :: filename


allocate(filename(1:nfile))

call allocate_unit(iunit)
open(unit=iunit, file=finpar, status='old')

   do i = 1, nfile, 1
      read(iunit, *) filename(i)
   end do

close(iunit)


fFDmodel = trim(filename(1))
fsource = trim(filename(2))
freceiver = trim(filename(3))
fFKmodel_left = trim(filename(4))
fFKmodel_right = trim(filename(5))
fVPmodel = trim(filename(6))
fVSmodel = trim(filename(7))

deallocate(filename)


end subroutine read_filename
!!=================================================================!!
subroutine get_parameter(nfile, nx, nz, norder, npml, izt, is_PML_top, outsnap, src_type, nlayer_left, &
                    nlayer_right, nstep, nrcv, x0, z0, dx, dz, dt, ds, f0, tmax, tstp, rayp, evla, evlo)

use constants, only: PI
use FDFK_par, only: prefix, fFDmodel, fsource, freceiver, &
                    fFKmodel_left, fFKmodel_right, lat_org, &
                    lon_org, az_org

implicit none

integer, intent(in) :: nfile

integer, intent(out) :: izt, nstep
integer, intent(out) :: norder, npml
integer, intent(out) :: nx, nz, nrcv
integer, intent(out) :: nlayer_left, nlayer_right

logical, intent(out) :: is_PML_top, outsnap, src_type

real(8), intent(out) :: dt

real, intent(out) :: tstp
real, intent(out) :: tmax, ds, f0
real, intent(out) :: x0, z0, dx, dz

real(8), intent(out) :: rayp, evla, evlo


integer iunit

real xn, zn


character(256) path


! FDmodel.dat
path = trim(prefix)//'/'//trim(fFDmodel)


call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')

   read(iunit, *)
   read(iunit, *) norder                          ! the order of FD
   read(iunit, *)
   read(iunit, *) x0, xn, z0, zn, dx, dz, dt
   read(iunit, *)
   read(iunit, *) tmax, tstp, npml, is_PML_top    ! the end time
   read(iunit, *)
   read(iunit, *) lat_org, lon_org, az_org
   read(iunit, *)
   read(iunit, *) outsnap, nstep

close(iunit)
norder = int(norder/2)
if (2*norder > (npml+1)) norder = int((npml+1)/2)


! Source.dat
path = trim(prefix)//'/'//trim(fsource)
call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')

   read(iunit, *)
   read(iunit, *) f0, ds, src_type                ! the dominant frequency and the strength of the source
   read(iunit, *)
   read(iunit, *) rayp, evla, evlo

close(iunit)
! from s/deg to s/m
rayp = rayp * 180.d0/(PI*6371.d3)


! Receiver.dat
path = trim(prefix)//'/'//trim(freceiver)
call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')
   read(iunit,*)
   read(iunit,*) nrcv
close(iunit)


! FKmodel_left.dat
path = trim(prefix)//'/'//trim(fFKmodel_left)
call allocate_unit(iunit)
open(unit=iunit, file=path, status = 'old')
   read(iunit,*)
   read(iunit,*) nlayer_left
close(iunit)


! FKmodel_right.dat
path = trim(prefix)//'/'//trim(fFKmodel_right)
call allocate_unit(iunit)
open(unit=iunit, file=path, status = 'old')
   read(iunit,*)
   read(iunit,*) nlayer_right
close(iunit)



nx = nint((xn-x0)/dx) + 1
nz = nint((zn-z0)/dz) + 1
if (is_PML_top) then
   izt = 1-npml
else
   izt = 1
end if


end subroutine get_parameter
!!=================================================================!!
subroutine read_data(nfile, nx, nz, npml, izt, nlayer_left, nlayer_right, nrcv, &
                          x0, z0, dx, dz, FKmodel_left, FKmodel_right, rx, rz, c)

use FDFK_par, only: prefix, freceiver, &
                    fFKmodel_left, fFKmodel_right, &
                    fVPmodel, fVsmodel

implicit none

integer, intent(in) :: nx, nz, izt
integer, intent(in) :: nfile, nrcv, npml
integer, intent(in) :: nlayer_left, nlayer_right

real, intent(in) :: x0, z0, dx, dz

real, intent(out) :: rx(1:nrcv)
real, intent(out) :: rz(1:nrcv)

real, intent(out) :: FKmodel_left(1:4, 1:nlayer_left)
real, intent(out) :: FKmodel_right(1:4, 1:nlayer_right)

real, intent(out) :: c(1:3, izt:nz+npml, 1-npml:nx+npml)


integer ircv
integer ix, iz, l
integer iunit, ilayer

character(256) path, filename

real xmin, xmax
real zmin, zmax
real rho, vp, vs

integer(2) head(1:120)



! read receivers' coordinates
path = trim(prefix)//'/'//trim(freceiver)
call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')
   do l = 1, 3, 1
      read(iunit,*)
   end do
   do ircv = 1, nrcv, 1
      read(iunit,*) rx(ircv), rz(ircv)
   end do
close(iunit)

! check whether the receiver outside the model
xmin = x0
xmax = x0 + (nx-1)*dx
zmin = z0
zmax = z0 + (nz-1)*dz
do ircv = 1, nrcv, 1
   if ( (rx(ircv) < xmin) .or. (rx(ircv) > xmax) .or. (rz(ircv) < zmin) .or. (rz(ircv) > zmax) ) then
      write(*,"(A, I0, A, F0.4, A, F0.4, A)") 'Error: ', ircv, &
           ' -th receiver at (xr, zr) = ( ', rx(ircv), ', ', rz(ircv), ' ) locates outside the model'
      stop
   end if
end do


! read FKmodel_left
path = trim(prefix)//'/'//trim(fFKmodel_left)
call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   do ilayer = 1, nlayer_left, 1
      read(iunit,*) vp, vs, FKmodel_left(1,ilayer)
      rho = (vp + 980.0) / 2.760
      !vs = vp / 1.732d0
      FKmodel_left(2,ilayer) = vp
      FKmodel_left(3,ilayer) = vs
      FKmodel_left(4,ilayer) = rho
   end do
close(iunit)


! read FKmodel_right
path = trim(prefix)//'/'//trim(fFKmodel_right)
call allocate_unit(iunit)
open(unit=iunit, file=path, status='old')
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   do ilayer = 1, nlayer_right, 1
      read(iunit,*) vp, vs, FKmodel_right(1,ilayer)
      rho = (vp + 980.0) / 2.760
      !vs = vp / 1.732d0
      FKmodel_right(2,ilayer) = vp
      FKmodel_right(3,ilayer) = vs
      FKmodel_right(4,ilayer) = rho
   end do
close(iunit)

if (maxval(abs(FKmodel_left(2:4,nlayer_left) - FKmodel_right(2:4,nlayer_right))) > 1.e-6) then
   write(*,"(A)") 'The bottom half-space must be identical !'
   stop
end if



! read vp velocity model
filename = trim(fVpmodel)
path = trim(prefix)//'/'//trim(filename)
call allocate_unit(iunit)
open(unit=iunit, file=path, form='unformatted', access='stream', status = 'old')
   do ix = 1-npml, nx+npml, 1
      read(iunit) head, (c(2,iz,ix), iz = izt, nz+npml)
   end do
close(iunit)


head = 0
head(58) = nz
head(59) = 1e3
!l = index(filename,'.') - 1
l = len_trim(filename) - 3
path = trim(prefix)//'/'//trim(filename(1:l))//'_without_pml.su'
call allocate_unit(iunit)
open(unit=iunit, file=path, form='unformatted', access='stream', status='unknown')
   do ix = 1, nx, 1
      write(iunit) head, (c(2,iz,ix), iz=1, nz)
   end do
close(iunit)


! read vs velocity model
filename = trim(fVsmodel)
path = trim(prefix)//'/'//trim(filename)
call allocate_unit(iunit)
open(unit=iunit, file=path, form='unformatted', access='stream', status = 'old')
   do ix = 1-npml, nx+npml, 1
      read(iunit) head, (c(3,iz,ix), iz = izt, nz+npml)
   end do
close(iunit)


head = 0
head(58) = nz
head(59) = 1e3
!l = index(filename,'.') - 1
l = len_trim(filename) - 3
path = trim(prefix)//'/'//trim(filename(1:l))//'_without_pml.su'
call allocate_unit(iunit)
open(unit=iunit, file=path, form='unformatted', access='stream', status='unknown')
   do ix = 1, nx, 1
      write(iunit) head, (c(3,iz,ix), iz=1, nz)
   end do
close(iunit)


if (any(c < 0.0)) then
   write(*,"(A)") "Error: negative values exist in rho- vp- or vs-models, please check !"
end if


!write(*,*) minval(c(2,:,:)), maxval(c(2,:,:))
!write(*,*) minval(c(3,:,:)), maxval(c(3,:,:))
!do ix = 1, nx, 1
!   do iz = 1, nz, 1
!      if (c(2,iz,ix) < 0.0) then
!         write(*,*) iz, ix, c(2,iz,ix)
!      end if
!   end do
!end do
!write(*,*) c(2,1,20:21)
!pause


!write(*,*) nx, nz, (nx+2*npml), (nz+npml-izt+1)
!pause
!$omp parallel default(shared) private(ix, iz, rho, vp, vs)
!$omp do schedule(dynamic)
do ix = 1-npml, nx+npml, 1
   do iz = izt, nz+npml, 1
      vp = c(2,iz,ix)
      vs = c(3,iz,ix)
      rho = (vp + 980.0) / 2.760
      !vs = vp / 1.7320
      c(1,iz,ix) = rho
      c(2,iz,ix) = rho * (vp*vp - 2.0*vs*vs)
      c(3,iz,ix) = rho * vs*vs
   end do
end do
!$omp end do
!$omp end parallel


end subroutine read_data
!!=================================================================!!
!!=================================================================!!
