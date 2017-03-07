!********************************************************************************
!************************ENVIRONMENT FOR SUBROUTINE PBC**************************
!********************************************************************************
!program pbc
!implicit none
!integer::t,i,j,N,DIM
!real:: BoxSize
!real,dimension(:,:),allocatable:: x
!real,dimension(:), allocatable:: distancia
!BoxSize=2.0
!DIM=3                      !dimension
!N=4                        !number of part.
!allocate(x(DIM,N))
!allocate(distancia(DIM))
! x=reshape((/(t,t=1,12)/),(/3,4/))

!do i=1,N-1                                 !loop over all pairs
 !do j=i+1,N
 !print*, x(:,i), "y" ,x(:,j)
 !call Refold_Positions(x(:,i),x(:,j),BoxSize,DIM,distancia)
 !end do
!end do 

!end program

subroutine Refold_Positions(posi,posj,BoxSize,DIM,distancia)
implicit none
integer::DIM
real:: BoxSize
real,dimension(DIM):: posi, posj,distancia

distancia=posi-posj !vector distancia entre particula i y j

distancia=distancia-BoxSize*nint(distancia/BoxSize) !condiciones periodicas

end subroutine Refold_Positions
