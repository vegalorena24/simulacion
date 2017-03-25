program hola

implicit none
integer,parameter       :: snaps=3,dimen=3,part=2
real(4)                 :: input(snaps*part,dimen)
real(4)                 :: acumulado(part)
integer                 :: i,j


open(unit=123,file='restart.rst',status='unknown',action='read')
  do i=1,snaps*part
    read(123,*), input(i,1), input(i,2), input(i,3)
  end do
close(123)

!do i=0,snaps-1
!  print *, " "
!  print *, "Snpashot ", i
!  do j=1+(i*part),part+(i*part)
!    print *, input(j,:)
!  end do
!end do

call postvisual(input,snaps,dimen,part,acumulado)

open(unit=124,file='desplazamientos.dat',status='unknown',action='write')
  do i=1,part
    write(124,*), acumulado(i)
  end do
close(124)


contains


subroutine postvisual(input,snaps,dimen,part,acumulado)
  integer,intent(in)      :: snaps,dimen,part
  real(4),intent(in)      :: input(snaps*part,dimen)
  real(4),intent(out)     :: acumulado(part)
  real(4)                 :: desplacement(part)
  integer                 :: frame
  integer                 :: i,j

  acumulado=0.
  desplacement=0.

  do i=0,snaps-1
    do j=1,part
      frame=j+(i*part)
      desplacement(j)=dist(input(j,:),input(frame,:))
      acumulado(j)=acumulado(j)+desplacement(j)
!      print *, "Desplazamiento de la particula ",j," en el frame ", i," : ", desplacement(j)
!      print *, "Desp. acum. de la particula ",j," en el frame ", i," : ", acumulado(j)
!      print *, " "
    end do
  end do

end subroutine postvisual


function dist(x,y)

  real(4),dimension(dimen)     :: x,y    ! Input vectors
  real(4),dimension(dimen)     :: vector ! Output vectors
  real(4)                      :: dist

  vector=x-y
!  print *, "Vector: ", vector

  dist=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!  print *, "Distancia: ", dist

end function dist


end program hola
