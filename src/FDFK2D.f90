!All copyright reserved (c)
!E-mail: ysliu@mail.iggcas.ac.cn
!QQ: 826897665
!Institute of Geology and Geophysics, Chinese Academy of Sciences
!Teleseismic wavefield modeling in 2-D istropic media using a hybrid method
!that coupling the finite difference (FD) method and frequency-wavenumber (FK) method
!
! Coordinates system
!
!     
!                x (R)
!                ^
!               /
!              /
!             /
!            /
!           /
!          /
!      O  +--------------------------> y (T)
!         |
!         |
!         |
!         |
!         |
!         |
!         |         
!         |
!         v
!         z
!     
!     
!!=================================================================!!
program FDFK2D

use FDFK_par

implicit none

integer nstep, izt
integer norder, npml
integer nargin, nrcv
integer nx, nz, nfile
integer nlayer_left, nlayer_right

logical is_PML_top, outsnap, src_type


real, allocatable :: rx(:)
real, allocatable :: rz(:)
real, allocatable :: c(:,:,:)
real, allocatable :: FKmodel_left(:,:)
real, allocatable :: FKmodel_right(:,:)

real(8) dt

real tstp
real tmax, f0, ds
real x0, z0, dx, dz

real(8) rayp, evla, evlo


integer(4) dtm_start(1:8)
integer(4) dtm_end(1:8)
integer(4) dtm_incr(1:8)

character(5)  myzone
character(8)  mydate
character(10) mytime


call date_and_time(mydate, mytime, myzone, dtm_start)



! parse command line parameters
nargin = command_argument_count()
if (nargin < 4) then
   write(*,"(A)") 'FDFK2D '//trim(ver)
   write(*,"(A)") './FDFK2D ./inpath inpar.dat ./outpath seism'
   write(*,"(A)")
   stop
end if

call get_command_argument(1, prefix)
call get_command_argument(2, finpar)
call get_command_argument(3, outpath)
call get_command_argument(4, fsis)


!! for test
!prefix = './input'
!finpar = 'inpar.dat'
!outpath = './seismograms'
!fsis = 'seis'


write(*,"(A)")
write(*,"(A)") 'This program simulates teleseismic wavefields with a hybrid method that '
write(*,"(A)") 'couples the Finite Difference Method (FD) with the Frequency-Wavenumber (FK) methods'
write(*,"(A)") 'Author: Youshan Liu [at] Institue of Geology and Geophysics, Chinese Academy of Sciences'
write(*,"(A)") 'All Copyright Reserved (C) ['//trim(ver)//']'
write(*,"(A)")

write(*,100) 'Start time: ', dtm_start(1), ' / ', dtm_start(2), ' / ', dtm_start(3), ' ; ', &
                  dtm_start(5), ' : ', dtm_start(6), ' : ', dtm_start(7), ' : ', dtm_start(8)


call get_file_num(nfile)

call read_filename(nfile)


call get_parameter(nfile, nx, nz, norder, npml, izt, is_PML_top, outsnap, src_type, nlayer_left, &
              nlayer_right, nstep, nrcv, x0, z0, dx, dz, dt, ds, f0, tmax, tstp, rayp, evla, evlo)


allocate(rx(1:nrcv))
allocate(rz(1:nrcv))

allocate(FKmodel_left(1:4, 1:nlayer_left))
allocate(FKmodel_right(1:4, 1:nlayer_right))

allocate(c(1:3, izt:nz+npml, 1-npml:nx+npml))

call read_data(nfile, nx, nz, npml, izt, nlayer_left, nlayer_right, nrcv, &
                    x0, z0, dx, dz, FKmodel_left, FKmodel_right, rx, rz, c)



call wavefieldsimulation(norder, nx, nz, npml, izt, is_PML_top, outsnap, src_type, nlayer_left, nlayer_right, nstep, &
               nrcv, x0, z0, dx, dz, dt, ds, f0, tmax, tstp, rayp, evla, evlo, FKmodel_left, FKmodel_right, rx, rz, c)


deallocate(FKmodel_left, FKmodel_right, rx, rz, c)


call date_and_time(mydate, mytime, myzone, dtm_end)

call elapsed_time(dtm_start, dtm_end, dtm_incr)


write(*,100) 'Start   time: ', dtm_start(1), ' / ', dtm_start(2), ' / ', dtm_start(3), ' ; ', &
                    dtm_start(5), ' : ', dtm_start(6), ' : ', dtm_start(7), ' : ', dtm_start(8)
write(*,200) 'End     time: ', dtm_end(1), ' / ', dtm_end(2), ' / ', dtm_end(3), ' ; ', &
                      dtm_end(5), ' : ', dtm_end(6), ' : ', dtm_end(7), ' : ', dtm_end(8)
write(*,300) 'Elapsed time: ', dtm_incr(1), ' / ', dtm_incr(2), ' / ', dtm_incr(3), ' ; ', &
                     dtm_incr(5), ' : ', dtm_incr(6), ' : ', dtm_incr(7), ' : ', dtm_incr(8)


100 format(A, I4, A3, I2, A3, I2, A3, I2, A3, I2, A3, I2, A3, I3)
200 format(A, I4, A3, I2, A3, I2, A3, I2, A3, I2, A3, I2, A3, I3)
300 format(A, I4, A3, I2, A3, I2, A3, I2, A3, I2, A3, I2, A3, I3)


end program FDFK2D
!!=================================================================!!
subroutine wavefieldsimulation(norder, nx, nz, npml, izt, is_PML_top, outsnap, src_type, nlayer_left, nlayer_right, nstep, &
                     nrcv, x0, z0, dx, dz, dt, ds, f0, tmax, tstp, rayp, evla, evlo, FKmodel_left, FKmodel_right, rx, rz, c)

use FDFK_par
use constants, only: PI, deg2rad, rad2deg

implicit none

integer, intent(in) :: nstep, izt
integer, intent(in) :: norder, npml
integer, intent(in) :: nx, nz, nrcv
integer, intent(in) :: nlayer_left, nlayer_right

real(8), intent(in) :: dt, evla, evlo

real, intent(in) :: tmax, f0, ds
real, intent(in) :: x0, z0, dx, dz

logical, intent(in) :: is_PML_top, outsnap, src_type

real, intent(in) :: rx(1:nrcv)
real, intent(in) :: rz(1:nrcv)

real, intent(in) :: FKmodel_left(1:4, 1:nlayer_left)
real, intent(in) :: FKmodel_right(1:4, 1:nlayer_right)

real, intent(in) :: c(1:3, izt:nz+npml, 1-npml:nx+npml)


real(8), intent(inout) :: rayp

real, intent(inout) :: tstp


integer ns, nsis, ircv
integer ix, iz, it, nt
integer isis, istp, nstp

integer(8) ntime1, ntime2, nrate

integer(2) head(1:120)

real vpml

character(8) answer
character(256) filename


! parameters for PMLs
integer npower, nbd, ier

real(8) epsil, R0, mem_size


real, allocatable, dimension(:) :: coordx

real, allocatable, dimension(:,:,:) :: ux, uz

real, allocatable, dimension(:,:) :: seisx, seisz

real, allocatable, dimension(:,:) :: bpxlx, bpzlx
real, allocatable, dimension(:,:) :: bpxrx, bpzrx

real, allocatable, dimension(:,:) :: bpxtx, bpxtz
real, allocatable, dimension(:,:) :: bpztx, bpztz
real, allocatable, dimension(:,:) :: bpxbx, bpxbz
real, allocatable, dimension(:,:) :: bpzbx, bpzbz

real, allocatable, dimension(:,:,:) :: bux1l, bux2l, bux3l
real, allocatable, dimension(:,:,:) :: buz1l, buz2l, buz3l
real, allocatable, dimension(:,:,:) :: bux1r, bux2r, bux3r
real, allocatable, dimension(:,:,:) :: buz1r, buz2r, buz3r

real, allocatable, dimension(:,:,:) :: bux1t, bux2t, bux3t
real, allocatable, dimension(:,:,:) :: buz1t, buz2t, buz3t
real, allocatable, dimension(:,:,:) :: bux1b, bux2b, bux3b
real, allocatable, dimension(:,:,:) :: buz1b, buz2b, buz3b


! arrays for background wavefields
real, allocatable, dimension(:,:,:) :: uxl_reg, uzl_reg
real, allocatable, dimension(:,:,:) :: uxr_reg, uzr_reg
real, allocatable, dimension(:,:,:) :: uxb_reg, uzb_reg
real, allocatable, dimension(:,:,:) :: uxl_pml, uzl_pml
real, allocatable, dimension(:,:,:) :: uxr_pml, uzr_pml
real, allocatable, dimension(:,:,:) :: uxb_pml, uzb_pml


real(8), allocatable, dimension(:) :: bazs

real(8), allocatable, dimension(:,:) :: coef0, coef1, coef2



allocate(coef0(1:norder, 1:norder))
allocate(coef1(1:norder, 1:norder))
allocate(coef2(1:norder, 1:norder))

allocate(ux(1:2, izt:nz+npml, 1-npml:nx+npml))
allocate(uz(1:2, izt:nz+npml, 1-npml:nx+npml))


allocate(bpxlx(1:nz, 1-npml:0))
allocate(bpzlx(1:nz, 1-npml:0))

allocate(bpxrx(1:nz, nx+1:nx+npml))
allocate(bpzrx(1:nz, nx+1:nx+npml))

allocate(bpxtx(1-npml:0, 1-npml:nx+npml))
allocate(bpxtz(1-npml:0, 1-npml:nx+npml))
allocate(bpztx(1-npml:0, 1-npml:nx+npml))
allocate(bpztz(1-npml:0, 1-npml:nx+npml))

allocate(bpxbx(nz+1:nz+npml, 1-npml:nx+npml))
allocate(bpxbz(nz+1:nz+npml, 1-npml:nx+npml))
allocate(bpzbx(nz+1:nz+npml, 1-npml:nx+npml))
allocate(bpzbz(nz+1:nz+npml, 1-npml:nx+npml))


allocate(bux1l(1:2, 1:nz, 1-npml:0))
allocate(bux2l(1:2, 1:nz, 1-npml:0))
allocate(bux3l(1:2, 1:nz, 1-npml:0))
allocate(buz1l(1:2, 1:nz, 1-npml:0))
allocate(buz2l(1:2, 1:nz, 1-npml:0))
allocate(buz3l(1:2, 1:nz, 1-npml:0))


allocate(bux1r(1:2, 1:nz, nx+1:nx+npml))
allocate(bux2r(1:2, 1:nz, nx+1:nx+npml))
allocate(bux3r(1:2, 1:nz, nx+1:nx+npml))
allocate(buz1r(1:2, 1:nz, nx+1:nx+npml))
allocate(buz2r(1:2, 1:nz, nx+1:nx+npml))
allocate(buz3r(1:2, 1:nz, nx+1:nx+npml))


allocate(bux1t(1:2, 1-npml:0, 1-npml:nx+npml))
allocate(bux2t(1:2, 1-npml:0, 1-npml:nx+npml))
allocate(bux3t(1:2, 1-npml:0, 1-npml:nx+npml))
allocate(buz1t(1:2, 1-npml:0, 1-npml:nx+npml))
allocate(buz2t(1:2, 1-npml:0, 1-npml:nx+npml))
allocate(buz3t(1:2, 1-npml:0, 1-npml:nx+npml))

allocate(bux1b(1:2, nz+1:nz+npml, 1-npml:nx+npml))
allocate(bux2b(1:2, nz+1:nz+npml, 1-npml:nx+npml))
allocate(bux3b(1:2, nz+1:nz+npml, 1-npml:nx+npml))
allocate(buz1b(1:2, nz+1:nz+npml, 1-npml:nx+npml))
allocate(buz2b(1:2, nz+1:nz+npml, 1-npml:nx+npml))
allocate(buz3b(1:2, nz+1:nz+npml, 1-npml:nx+npml))




nt = nint(tmax / dt)
if (tstp < dt) tstp = dt
nstp = nint(tstp / dt)
nsis = int(nt / nstp )
tstp = nstp * dt
isis = 0
istp = 0


nbd = 2*norder-1
mem_size = 4_8*int8(nt+1)*(4_8*int8((nz-nbd)*nbd) + 2_8*int8(nx*nbd) + &
                           4_8*int8(nz*nbd) + 2_8*int8((nx+2*nbd)*nbd)) / 1.0737418240d9

write(*,*)
write(*,"(A, F7.4, A)") 'INFO: it requires ', mem_size, ' GBs to store the background wavefields for Hybrid method !'
write(*,*)


! arrays for background wavefields
allocate(uxl_reg(1:nz-nbd, 1:nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxl_reg in line 253 of FDWFS.f90 !'
allocate(uzl_reg(1:nz-nbd, 1:nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uzl_reg in line 255 of FDWFS.f90 !'
allocate(uxr_reg(1:nz-nbd, nx-nbd+1:nx, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxr_reg in line 257 of FDWFS.f90 !'
allocate(uzr_reg(1:nz-nbd, nx-nbd+1:nx, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxr_reg in line 259 of FDWFS.f90 !'
allocate(uxb_reg(nz-nbd+1:nz, 1:nx, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxb_reg in line 261 of FDWFS.f90 !'
allocate(uzb_reg(nz-nbd+1:nz, 1:nx, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxb_reg in line 263 of FDWFS.f90 !'


allocate(uxl_pml(1:nz, 1-nbd:0, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxl_pml in line 267 of FDWFS.f90 !'
allocate(uzl_pml(1:nz, 1-nbd:0, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uzl_pml in line 269 of FDWFS.f90 !'
allocate(uxr_pml(1:nz, nx+1:nx+nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxr_pml in line 271 of FDWFS.f90 !'
allocate(uzr_pml(1:nz, nx+1:nx+nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxr_pml in line 273 of FDWFS.f90 !'
allocate(uxb_pml(nz+1:nz+nbd, 1-nbd:nx+nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uxb_pml in line 275 of FDWFS.f90 !'
allocate(uzb_pml(nz+1:nz+nbd, 1-nbd:nx+nbd, 0:nt),stat=ier)
if(0 /= ier) write(*,'(A)') 'ERROR: Fail to allocate array for uzb_pml in line 277 of FDWFS.f90 !'



! initalization
bpxlx = 0.0; bpzlx = 0.0
bpxrx = 0.0; bpzrx = 0.0
bpxtx = 0.0; bpxtz = 0.0
bpztx = 0.0; bpztz = 0.0
bpxbx = 0.0; bpxbz = 0.0
bpzbx = 0.0; bpzbz = 0.0

ux = 0.0; uz = 0.0
bux1l = 0.0; bux2l = 0.0; bux3l = 0.0
buz1l = 0.0; buz2l = 0.0; buz3l = 0.0
bux1r = 0.0; bux2r = 0.0; bux3r = 0.0
buz1r = 0.0; buz2r = 0.0; buz3r = 0.0

bux1t = 0.0; bux2t = 0.0; bux3t = 0.0
buz1t = 0.0; buz2t = 0.0; buz3t = 0.0
bux1b = 0.0; bux2b = 0.0; bux3b = 0.0
buz1b = 0.0; buz2b = 0.0; buz3b = 0.0



! compute backazimuths of the grids
allocate(coordx(1-npml:nx+npml))
do ix = 1-npml, nx+npml, 1
   coordx(ix) = x0 + (ix-1)*dx
end do
allocate(bazs(1-npml:nx+npml))
call compute_backazimuth_grids(nx+2*npml, coordx, evla, evlo, bazs)
bazs = bazs - az_org*deg2rad
if (abs(cos(sum(bazs)/(nx+2*npml+1))) < 0.5d0*sqrt(3.d0)) then
   write(*,"(A)") 'WRANING: It cannot accurately simulate the teleseismic wavefield using a 2D code'
   write(*,"(A)") 'when the angle between azimuth of the profile and backazimuth of station to event is large than 30 [deg].'
   write(*,"(A)", advance='no') 'Do you want to continue ? [Y/N]: '
   read(*,*) answer
   ix = index(answer,'N')
   if (ix > 0) then
      return
   end if
   ix = index(answer,'n')
   if (ix > 0) then
      return
   end if
end if
deallocate(coordx)



! compute FD cofficients
call fd_coef0(norder, coef0)
call fd_coef1(norder, coef1)
call fd_coef2(norder, coef2)
! check stability
call check_stability(nx, nz, npml, izt, norder, dx, dz, dt, coef2, c)


!================================================================================!
!================================================================================!
vpml = maxval( sqrt( (c(2,:,:) + 2.0*c(3,:,:)) / c(1,:,:) ) )
npower = 2
!epsil = 0.65d0
!R0 = dexp( -2.d0*(npml*dx)*(1.d0 - epsil) / (1.d0 + epsil) / ((npower+1)*dt*vpml) )
!R0 = 5.d-2
R0 = 1.d-3
write(*,"(A, ES12.4)") 'R0 = ', R0
!================================================================================!
!================================================================================!



write(*,*)
write(*,"(A)") 'computing background wavefield using FK method... '
write(*,*)
if (abs(rayp) < 1.d-9) rayp = 1.d-9*sign(1.d0, rayp)
call comput_background_wavefield(nx, nz, npml, izt, nbd, nt, nlayer_left, nlayer_right, src_type, &
           x0, z0, dx, dz, dt, ds, f0, rayp, bazs, FKmodel_left, FKmodel_right, uxl_reg, uzl_reg, &
           uxr_reg, uzr_reg, uxb_reg, uzb_reg, uxl_pml, uzl_pml, uxr_pml, uzr_pml, uxb_pml, uzb_pml)
deallocate(bazs)



write(*,*)
write(*,*)
write(*,"(A)") 'computing wavefield with FD-FK hybrid method ... '
write(*,*)
write(*,"(A, I0, 1x, A, 1x, I0)") 'nrcv = ', nx, ' ; nsamp = ', nsis+1
write(*,*)
write(*,"(A, F0.4, A, F0.4, A, F0.4, A, F9.6)") 'f0 = ', f0, ' ; dx = ', dx, ' ; dz = ', dz, ' ; dt =', dt
write(*,*)
write(*,"(A, I0, A, I0, A, I0)") 'nx = ', nx, ' ; nz = ', nz, ' ; nt = ', nt+1
write(*,*)
write(*,*)


allocate(seisx(0:nsis, 1:nrcv))
allocate(seisz(0:nsis, 1:nrcv))




head = 0
head(58) = nz
head(59) = dz
!head(58) = nz+npml
call system_clock(ntime1, nrate)
!open(41, file='seisux.dat')
!open(42, file='seisuz.dat')
do it = 0, nt, 1                       !iteration

   if (outsnap) then
      if( (it > 0) .and. (0 == mod(it,nstep)) )then
         write(filename,"(A,I5.5,A)") './snapshots/', it, 'ux.su'
         open(31, file=filename, form='unformatted', access='stream')
            do ix = 1, nx, 1
               write(31) head, (ux(2,iz,ix), iz = 1, nz)
            end do
         close(31)
         write(filename,"(A,I5.5,A)") './snapshots/', it ,'uz.su'
         open(32, file=filename, form='unformatted', access='stream')
            do ix = 1, nx, 1
               write(32) head, (uz(2,iz,ix), iz = 1, nz)
            end do
         close(32)
         !write(filename,"(A,A,I5.5,A)") './snapshots/', '/', it, 'ux.su'
         !open(31, file=filename, form='unformatted', access='stream')
         !   do ix = 1-npml, nx+npml, 1
         !      write(31) head, (ux(2,iz,ix), iz = 1, nz+npml)
         !   end do
         !close(31)
         !write(filename,"(A,A,I5.5,A)") './snapshots/', '/', it ,'uz.su'
         !open(32, file=filename, form='unformatted', access='stream')
         !   do ix = 1-npml, nx+npml, 1
         !      write(32) head, (uz(2,iz,ix), iz = 1, nz+npml)
         !   end do
         !close(32)
      end if
   end if


   if ( (isis == it) .and. (istp <= nsis) ) then
      !seisx(istp,1:nx) = ux(1,izrcv,1:nx)
      !seisz(istp,1:nx) = uz(1,izrcv,1:nx)
      !$omp parallel default(shared) private(ircv)
      !$omp do schedule(dynamic)
      do ircv = 1, nrcv, 1
         call loadseis(nx, nz, npml, izt, x0, z0, dx, dz, rx(ircv), rz(ircv), ux, uz, seisx(istp,ircv), seisz(istp,ircv))
      end do
      !$omp end do
      !$omp end parallel
      istp = istp + 1
      isis = isis + nstp
   end if

   !do ix = 1, nx, 1
   !   write(41,*) (ix-1)*dx, it*dt, ux(2,izrcv,ix)
   !   write(42,*) (ix-1)*dx, it*dt, uz(2,izrcv,ix)
   !end do

   call updatedispl(norder, nx, nz, npml, izt, is_PML_top, npower, nbd, it, nt, &
             dx, dz, dt, vpml, R0, uxl_reg, uzl_reg, uxr_reg, uzr_reg, uxb_reg, &
             uzb_reg, uxl_pml, uzl_pml, uxr_pml, uzr_pml, uxb_pml, uzb_pml, c, &
             coef0, coef1, ux, uz, bpxlx, bpzlx, bpxrx, bpzrx, bpxtx, bpxtz, &
             bpztx, bpztz, bpxbx, bpxbz, bpzbx, bpzbz, bux1l, bux2l, bux3l, buz1l, &
             buz2l, buz3l, bux1r, bux2r, bux3r, buz1r, buz2r, buz3r, bux1b, bux2b, &
             bux3b, buz1b, buz2b, buz3b, bux1t, bux2t, bux3t, buz1t, buz2t, buz3t)

   if(0 == mod(it, 50))then
      call system_clock(ntime2)
      write(*, '(A, F10.5, A, ES12.4, ES12.4, 1x, A3, F8.3, A)') 't =  ', &
              it*dt, ' s', maxval(abs(ux(1,:,:))), maxval(abs(uz(1,:,:))), ' ; ', &
              (nt-it)*real(ntime2-ntime1) / real(nrate) / 3000.d0, ' min. remained'
      ntime1 = ntime2
   end if

end do
!close(41)
!close(42)


allocate(bazs(1:nrcv))
call compute_backazimuth_grids(nrcv, rx, evla, evlo, bazs)
bazs = bazs - az_org*deg2rad
do ircv = 1, nrcv, 1
   seisx(:,ircv) = seisx(:,ircv)*cos(bazs(ircv))
end do
deallocate(bazs)



head = 0
head(58) = nsis+1
head(59) = tstp*1e6
open(31, file=trim(outpath)//'/'//trim(fsis)//'x.su', form='unformatted', access='stream')
close(31, status='delete')
open(32, file=trim(outpath)//'/'//trim(fsis)//'z.su', form='unformatted', access='stream')
close(32, status='delete')
open(31, file=trim(outpath)//'/'//trim(fsis)//'x.su', form='unformatted', access='stream')
open(32, file=trim(outpath)//'/'//trim(fsis)//'z.su', form='unformatted', access='stream')
   do ircv = 1, nrcv, 1
      write(31) head, (seisx(it,ircv), it = 0, nsis)
      write(32) head, (seisz(it,ircv), it = 0, nsis)
   end do
close(31)
close(32)
pause 3


deallocate(ux, uz)
deallocate(seisx, seisz)

deallocate(bpxlx, bpzlx)
deallocate(bpxrx, bpzrx)
deallocate(bpxtx, bpxtz)
deallocate(bpztx, bpztz)
deallocate(bpxbx, bpxbz)
deallocate(bpzbx, bpzbz)

deallocate(coef0, coef1, coef2)

deallocate(bux1l, bux2l, bux3l)
deallocate(buz1l, buz2l, buz3l)
deallocate(bux1r, bux2r, bux3r)
deallocate(buz1r, buz2r, buz3r)

deallocate(bux1t, bux2t, bux3t)
deallocate(buz1t, buz2t, buz3t)
deallocate(bux1b, bux2b, bux3b)
deallocate(buz1b, buz2b, buz3b)


end subroutine wavefieldsimulation
!!=================================================================!!
