program hola

implicit none
integer,parameter       :: snaps=3,dimen=3,part=2
real(4)                 :: input(snaps*part,dimen)
real(4)                 :: desplacement(snaps,part)
integer                 :: i,j

! Read coordinates of the trajectory (restart xyz file)
open(unit=123,file='restart.rst',status='unknown',action='read')
  do i=1,snaps*part
    read(123,*), input(i,1), input(i,2), input(i,3)
  end do
close(123)


! Calculate the rmsd for each particle
call postvisual(input,snaps,dimen,part,desplacement)


! Write output file
open(unit=124,file='desplazamientos.dat',status='unknown',action='write')
  do i=1,snaps
    write(124,*), desplacement(i,:)
  end do
close(124)



contains

subroutine postvisual(input,snaps,dimen,part,desplacement)
!!!
!!! Computes the RMSD for each particle.
!!!
  integer,intent(in)      :: snaps,dimen,part
  real(4),intent(in)      :: input(snaps*part,dimen)
  real(4),intent(out)     :: desplacement(snaps,part)
  integer                 :: frame
  integer                 :: i,j

  desplacement=0.

  do j=1,part
    do i=1,snaps
      frame=j+((i-1)*part)
      desplacement(i,j)=dist(input(j,:),input(frame,:))
    end do
  end do


end subroutine postvisual


function dist(x,y)
!!!
!!! Computes the distance between two points.
!!!
  real(4),dimension(dimen)     :: x,y    ! Input vectors
  real(4),dimension(dimen)     :: vector ! Output vectors
  real(4)                      :: dist

  vector=x-y

  dist=vector(1)**2+vector(2)**2+vector(3)**2

end function dist


end program hola
