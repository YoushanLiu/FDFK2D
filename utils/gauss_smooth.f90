program gauss_smooth

!use dflib

implicit none

integer imethod, nlen
integer i1, j1, k1, k2
integer narg, nz0, nz1
integer i, j, k, ix, iz
integer nx0, nx1, iunit
integer nx, nz, nwx, nwz

integer(2) head(1:120)

real, allocatable :: w(:,:)
real, allocatable :: vel(:,:)
real, allocatable :: vel2(:,:)

character(256) infile
character(256) str_dx
character(256) str_dz
character(256) outfile
character(256) str_nz0
character(256) str_nz1
character(256) str_nx0
character(256) str_nx1
character(256) str_method
character(256) str_lambdax
character(256) str_lambdaz

real x, z, sigma2
real dx, dz, mysum
real dx2, dz2, mywgt
real lambdax, lambdaz

real(8), parameter :: PI = 4.d0*datan(1.d0)


logical existed


narg = iargc()
!narg = get_command_count()


if(narg < 1)then
   !write(*,"(A)") 'At least the filename is required '
   !write(*,'(A)') ' filename,  lambdax,  lambdaz,  dx,  dz,  iz_first,  iz_end,  ix_first,  ix_end,  method '
   print*, 'Usage: '
   print*, 'Author: Youshan Liu'
   write(*,'(A)') ' gauss_smooth(filename, lambdax, lambdaz, dx, dz, iz_first, iz_end, ix_first, ix_end, method)'
   print*, 'filename: input filename (mandatory)'
   print*, 'lambdax: horizontal correlation length (optional)'
   print*, 'lambdaz: vertical correlation length (optional)'
   print*, 'dx: horizontal interval (optional)'
   print*, 'dz: vertical interval (optional)'
   print*, 'iz_first: lower bound of the smooth window in vertical (default 1, optional)'
   print*, '    when iz_first <= 0 or iz_first > nz, iz_first will be set 1'
   print*, 'iz_end: upper bound of the smooth window in vertical (default nz, optional)'
   print*, '    when iz_end   <= 0 or iz_end   > nz, iz_first will be set nz'
   print*, 'ix_first: lower bound of the smooth window in horizontal (default 1, optional)'
   print*, '    when ix_first <= 0 or ix_first > nx, iz_first will be set 1'
   print*, 'ix_end: lower bound of the smooth window in horizontal (default nx, optional)'
   print*, '    when ix_end   <= 0 or ix_end   > nx, iz_first will be set nx'
   print*, 'method: gauss or triangular (default 0 for gaussian smooth otherwise 1 for triangular smooth, optional)'
   write(*,'(A)')
   stop
end if
!call getarg(1, infile)
call get_command_argument(1, infile)


call get_parameter(nx, nz, infile)
allocate(vel(1:nx,1:nz))
call read_data(nx, nz, dz, vel, infile)

if(narg <= 3)then
   dx = dz
end if

if (1 ==narg) then
   lambdax = dz
   lambdaz = lambdax
   nz0 = 1
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (2 == narg) then
   !call getarg(2, str_lambdax)
   call get_command_argument(1, str_lambdax)
   read(str_lambdax,*) lambdax
   lambdaz = lambdax
   nz0 = 1
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (3 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   nz0 = 1
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (4 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx,*) dx
   nz0 = 1
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (5 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx,*) dx
   read(str_dz,*) dz
   nz0 = 1
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (6 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   !call getarg(6, str_nz0)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   call get_command_argument(6, str_nz0)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx, *) dx
   read(str_dz, *) dz
   read(str_nz0,*) nz0
   nz1 = nz
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (7 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   !call getarg(6, str_nz0)
   !call getarg(7, str_nz1)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   call get_command_argument(6, str_nz0)
   call get_command_argument(7, str_nz1)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx, *) dx
   read(str_dz, *) dz
   read(str_nz0,*) nz0
   read(str_nz1,*) nz1
   nx0 = 1
   nx1 = nx
   imethod = 0
else if (8 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   !call getarg(6, str_nz0)
   !call getarg(7, str_nz1)
   !call getarg(8, str_nx0)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   call get_command_argument(6, str_nz0)
   call get_command_argument(7, str_nz1)
   call get_command_argument(8, str_nx0)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx, *) dx
   read(str_dz, *) dz
   read(str_nz0,*) nz0
   read(str_nz1,*) nz1
   read(str_nx0,*) nx0
   nx1 = nx
   imethod = 0
else if (9 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   !call getarg(6, str_nz0)
   !call getarg(7, str_nz1)
   !call getarg(8, str_nx0)
   !call getarg(9, str_nx1)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   call get_command_argument(6, str_nz0)
   call get_command_argument(7, str_nz1)
   call get_command_argument(8, str_nx0)
   call get_command_argument(9, str_nx1)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx, *) dx
   read(str_dz, *) dz
   read(str_nz0,*) nz0
   read(str_nz1,*) nz1
   read(str_nx0,*) nx0
   read(str_nx1,*) nx1
   imethod = 0
else if (10 == narg) then
   !call getarg(2, str_lambdax)
   !call getarg(3, str_lambdaz)
   !call getarg(4, str_dx)
   !call getarg(5, str_dz)
   !call getarg(6, str_nz0)
   !call getarg(7, str_nz1)
   !call getarg(8, str_nx0)
   !call getarg(9, str_nx1)
   !call getarg(10, str_method)
   call get_command_argument(2, str_lambdax)
   call get_command_argument(3, str_lambdaz)
   call get_command_argument(4, str_dx)
   call get_command_argument(5, str_dz)
   call get_command_argument(6, str_nz0)
   call get_command_argument(7, str_nz1)
   call get_command_argument(8, str_nx0)
   call get_command_argument(9, str_nx1)
   call get_command_argument(10, str_method)
   read(str_lambdax,*) lambdax
   read(str_lambdaz,*) lambdaz
   read(str_dx, *) dx
   read(str_dz, *) dz
   read(str_nz0,*) nz0
   read(str_nz1,*) nz1
   read(str_nx0,*) nx0
   read(str_nx1,*) nx1
   read(str_method,*) imethod
end if


if (nx0 <= 0) nx0 = 1
if (nx0 > nx) nx0 = 1
if (nx1 <= 0) nx1 = nx
if (nx1 > nx) nx1 = nx
if (nz0 <= 0) nz0 = 1
if (nz0 > nz) nz0 = 1
if (nz1 <= 0) nz1 = nz
if (nz1 > nz) nz1 = nz

if (nx0 > nx1) then
   ix = nx0
   nx0 = nx1
   nx1 = ix
   nx0 = max(1 , nx0)
   nx1 = min(nx, nx1)
end if
if (nz0 > nz1) then
   iz = nz0
   nz0 = nz1
   nz1 = iz
   nz0 = max(1 , nz0)
   nz1 = min(nz, nz1)
end if

if (imethod > 1) imethod = 0



if (0 == imethod) then
   write(*,100) trim(adjustl(infile)), lambdax, lambdaz, dx, dz, nz0, nz1, nx0, nx1, ' gauss'
else
   write(*,100) trim(adjustl(infile)), lambdax, lambdaz, dx, dz, nz0, nz1, nx0, nx1, ' triangular'
end if

write(*,'(A)')



dx2 = dx*dx
dz2 = dz*dz
!nwx = 4*ceiling(lambdax/dx)
!nwz = 4*ceiling(lambdaz/dz)

!allocate(w(-nwx:nwx,-nwz:nwz))


if (0 == imethod) then

   nwx = ceiling(3.6*lambdax/dx)
   nwz = ceiling(3.6*lambdaz/dz)
   allocate(w(-nwx:nwx,-nwz:nwz))

   w = 0.d0
   !$omp parallel default(shared) private(j, k, x, z)
   !$omp do schedule(dynamic)
   do j = -nwx, nwx, 1
      x = j*dx
      do k = -nwz, nwz, 1
         z = k*dz
         w(j,k) = dexp( -0.5d0*( (x/lambdax)*(x/lambdax) + (z/lambdaz)*(z/lambdaz) ) )
         !w(j,k) = exp(-dble(j*j*dx2+k*k*dz2)*sigma)
      end do
   end do
   !$omp end do
   !$omp end parallel
   w = w/sum(w)

else if (1 == imethod) then

   nwx = ceiling(2.1*lambdax/dx)
   nwz = ceiling(2.1*lambdaz/dz)
   allocate(w(-nwx:nwx,-nwz:nwz))

   w = 0.d0
   sigma2 = dble(1.8d0*ceiling(lambdaz/dz))
   !$omp parallel default(shared) private(j, k, x, z)
   !$omp do schedule(dynamic)
   do j = -nwx, nwx, 1
      do k = -nwz, nwz, 1
         w(j,k) = dble(nwz - abs(k)) / sigma2
      end do
      !do k = -nwz, 0, 1
      !   w(j,k) = dble(nwz+k)/sigma2
      !end do
      !do k = 1, nwz, 1
      !   w(j,k) = dble(nwz-k)/sigma2
      !end do
   end do
   !$omp end do
   sigma2 = dble(1.8d0*ceiling(lambdax/dx))
   !$omp end parallel
   do k = -nwz, nwz, 1
      do j = -nwx, nwx, 1
         w(j,k) = dble(nwx - abs(j)) / sigma2 * w(j,k)
      end do
      !do j = -nwx, 0, 1
      !   w(j,k) = dble(j+nwx)/dble(nwx)*w(j,k)
      !end do
      !do j = 1, nwx, 1
      !   w(j,k) = dble(nwx-j)/dble(nwx)*w(j,k)
      !end do
   end do
   w = w/sum(w)

end if



vel2 = vel
!allocate(vel2(1:nx,1:nz))

!$omp parallel default(shared) private(i, j, i1, j1, k1, k2, mysum, mywgt)
!$omp do schedule(dynamic)
do i = nx0, nx1, 1
   do j = nz0, nz1, 1
      mysum = 0.d0
      mywgt = 0.d0
      do i1 = -nwx, nwx, 1
         do j1 = -nwz, nwz, 1
            k1 = i+i1
            k2 = j+j1
            if(k1 < 1)  k1 = 1-k1
            if(k1 > nx) k1 = 2*nx-k1+1
            if(k2 < 1)  k2 = 1-k2
            if(k2 > nz) k2 = 2*nz-k2+1
            mywgt=mywgt + w(i1,j1)
            mysum=mysum + w(i1,j1)*vel(k1,k2)
         end do
      end do
      vel2(i,j) = mysum / mywgt
   end do
end do
!$omp end do
!$omp end parallel


nlen = index(infile,'.') - 1
write(outfile,'(A, A)') infile(1:nlen), '_smooth.su'

inquire(file=outfile, exist=existed)
if(existed)then
   open(unit=11, file=outfile)
   close(11, status='delete')
end if

head = 0
head(58) = nz
head(59) = dz*1e3


iunit = 25
open(unit=iunit, file=outfile, form='unformatted', access='stream', status='unknown')
   do j = 1, nx, 1
      write(iunit) head, (sngl(vel2(j,k)), k = 1, nz)
   end do
close(iunit)

write(*,200) maxval(abs(vel-vel2)), sqrt(sum(((vel-vel2)/vel)**2/size(vel)))*100.0
write(*,'(A)')

deallocate(w,vel,vel2)

!100 format(A, F12.2, F12.2, I8, F10.3, F10.3)
100 format('file = ', A, ' ; lambdax = ', F0.2, ' ; lambdaz = ', F0.2, ' ; dx = ', F0.2, ' ; dz = ', F0.2, / &
'iz_first = ', I0, ' ; iz_end = ', I0, ' ; ix_first = ', I0, ' ; ix_end = ', I0, ' ; method = ', A)
200 format('After gauss_smooth: maximum absolute error = ', f0.2,' ; maxmum relative error = ', f0.2,' %')

end program gauss_smooth
!=============================================================
!=============================================================
subroutine get_parameter(ntrace, ntime, filename)

!use ifport

implicit none

integer ntrace, ntime
integer ier, iskip, iunit

integer fseek, ftell

integer(8) nbytes

integer(2) head(1:120)

character*(*) filename


head = 0
iunit = 21
open(unit=iunit, file=filename, status='unknown', form='unformatted', access='stream')
   read(iunit) head
   ntime = head(58)
   !ier = fseeki8(iunit, 0, 2)
   !nbytes = ftelli8(iunit)
   ier = fseek(iunit, 0, 2)
   nbytes = ftell(iunit)
   ntrace = int8(nbytes/(240 + 4*int8(ntime)))
close(iunit)

write(*,'(A, I0)') 'The number of trace is: ', ntrace
write(*,'(A, I0)') 'The number of sampling points is: ', ntime
write(*,*)

end subroutine get_parameter
!=============================================================
subroutine read_data(ntrace, ntime, dt, seis, filename)

implicit none

integer itr, it, iunit
integer ntrace, ntime

integer(2) head(1:120)

character*(*) filename

real dt
real seis(1:ntrace, 1:ntime)


iunit = 21
open(unit=iunit, file=filename, status='unknown', form='unformatted', access='stream')

   do itr = 1, ntrace, 1
      read(iunit) head, (seis(itr,it), it = 1, ntime)
   end do

close(iunit)

dt = head(59) ! *1.e-3

end subroutine read_data
!=============================================================
!=============================================================