module albert
use paralelizar
contains
! ==== POSTVISUALIZATION ================================================================
subroutine postvisual(snaps,dimen,part)


! Routine variables
integer,intent(in)      :: snaps,dimen             ! Number of snapshots; dimensions 
integer,intent(in)      :: part                    ! Number of particles            
real(8)                 :: input(snaps*part,dimen) ! Input data    
real(8)                 :: desplacement(snaps,part)   ! Dummy variable
integer                 :: i,j,ip                  ! Iterators
integer                 :: a,b                     ! Index variable
! MPI variables
integer                 :: partxproc               ! Number of particles per processor


 MASTER = 0

! Read coordinates of the trajectory (restart xyz file)
  open(unit=126,file='restart.rst',status='unknown',action='read')
    print *, "Leyendo"
    do i=1,snaps*part
      read(126,*), input(i,1), input(i,2), input(i,3)
    end do
  close(126)

! Distribution of the work per each processor
  partxproc=nint(real(part)/real(numproc))
  do i=0,numproc-2
    ini(i)=i*partxproc+1
    fin(i)=(i+1)*partxproc
  end do
  ini(numproc-1)=(numproc-1)*partxproc+1
  fin(numproc-1)=part


! Calculate the rmsd for each particle
print *, "Calculamos", rank
call desplace(input,snaps,dimen,part,desplacement,ini(rank),fin(rank),part)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! Send the calculated distances to the MASTER
if (rank /= MASTER) then
  print *, "Rank --> master"
  a=ini(rank)
  b=fin(rank)
  call MPI_ISEND(desplacement(:,a:b),snaps*(b-a+1),MPI_DOUBLE_PRECISION,MASTER,2,&
MPI_COMM_WORLD,request,ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD,ierror)

if ( rank == MASTER ) then
  do ip=1,numproc-1
    print *, "Master <-- rank", rank
    a=ini(ip)
    b=fin(ip)
    call MPI_RECV(desplacement(:,a:b),(snaps*(b-a+1)),MPI_DOUBLE_PRECISION,ip,2,&
MPI_COMM_WORLD,stat,ierror)
  end do

! Write output file
  open(unit=125,file='desplazamientos.dat',status='unknown',action='write')
    print *, "Imprimiendo"
    do i=1,snaps
      write(125,*), desplacement(i,:)
    end do
  close(125)
end if
call MPI_BARRIER(MPI_COMM_WORLD,ierror)



end subroutine postvisual


subroutine desplace(input,snaps,dimen,part,desplacement,from_i,to_j,n)
!!!
!!! Computes the RMSD for each particle.
!!!
  integer,intent(in)      :: n                          ! Total particles
  integer,intent(in)      :: from_i, to_j               ! Initial particle; final particle
  integer,intent(in)      :: snaps,dimen,part           ! Number of snapshots; dimensions; number of particles
  real(8),intent(in)      :: input(snaps*n,dimen)       ! Coordinates array
  real(8),intent(out)     :: desplacement(snaps,part)   ! Dummy variable
  integer                 :: frame                      ! Frame computed
  integer                 :: i,j                        ! Iterators

  desplacement=0.

  do j=from_i,to_j
    do i=1,snaps
      frame=j+((i-1)*part)
      desplacement(i,j)=dist(input(j,:),input(frame,:),dimen)
    end do
  end do

end subroutine desplace


function dist(x,y,dimen)
!!!
!!! Computes the distance between two points.
!!!
  integer                      :: dimen  ! dimension 
  real(8),dimension(dimen)     :: x,y    ! Input vectors
  real(8),dimension(dimen)     :: vector ! Output vectors
  real(8)                      :: dist

  vector=x-y

  dist=vector(1)**2+vector(2)**2+vector(3)**2

end function dist
end module albert

! ======================================================================================
