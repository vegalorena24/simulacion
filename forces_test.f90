program thest

implicit none

integer                                 :: N, i
real, dimension(:,:), allocatable       :: pos, forces_mat
real                                    :: epot, L, start, end

open(unit=123, file='forces_positions.xyz', status='old', action='read')
read(unit=123, fmt=*) N, L

allocate(pos(N,3),forces_mat(N,3))
do i = 1, N
    read(unit=123, fmt=*)   pos(i,:)
enddo
close(unit=123)

call cpu_time(start)
call    forces(pos,L,forces_mat,L,epot)
call cpu_time(end)

print*, "cpu_time = ", end-start

open(unit=124,file='forces_serie.dat',status='replace',action='write')

do i = 1, N
    write(unit=124,fmt='(3f20.10)') forces_mat(i,:)
enddo
close(unit=124)
contains

subroutine forces(pos,boxlength,accel,rc,epot)
    real, dimension(:,:), intent(in)  :: pos
    real, dimension(:,:), intent(out) :: accel
    real, intent(in)                  :: boxlength, rc
    real, intent(out)                 :: epot
    integer                             :: is, js, natoms
    real                             :: pot

    natoms = size(pos,1)

    accel = 0.0d0
    epot = 0.d0

    !      atom-atom interactions

    do is = 1,natoms-1
       do js = is+1,natoms
          call lj(is,js,pos,boxlength,accel,rc,pot)
          epot = epot + pot
       end do
    end do


end subroutine

subroutine lj(is,js,pos,boxlength,accel,rc,pot)
real, dimension(:,:), intent(in)      :: pos
real, dimension(:,:), intent(inout)   :: accel
real, intent(in)                      :: boxlength, rc
integer, intent(in)                     :: is, js
real, intent(out)                     :: pot
real, dimension(size(pos,2))          :: rij
real                                  :: rr2, rijl, rr, forcedist
integer                                 :: l, dim

dim = size(pos,2)

rr2 = 0.d0
pot = 0.d0
do l = 1,dim
   rijl = pos(js,l) - pos(is,l)
   rij(l) = rijl - boxlength*nint(rijl/boxlength)
   rr2 = rr2 + rij(l)*rij(l)
end do

rr = sqrt(rr2)

if (rr.lt.rc) then
    forcedist = 24.d0*(2.d0/rr**14-1.0d0/rr**8)
    pot = 4.d0*(1.0d0/rr**12-1.0d0/rr**6)
    do l = 1,dim
        accel(is,l) = accel(is,l) - forcedist*rij(l)
        accel(js,l) = accel(js,l) + forcedist*rij(l)
    end do
end if

end subroutine
end program
