program postvis

  include 'mpif.h'

! Routine variables
integer,parameter       :: snaps=3,dimen=3         ! Number of snapshots; dimensions | To be
integer,parameter       :: part=2                  ! Number of particles             | read.
real(8)                 :: input(snaps*part,dimen) ! Input data    
real(8)                 :: desplacement(snaps,part)   ! Dummy variable
integer                 :: i,j,ip                  ! Iterators
integer                 :: a,b                     ! Index variable
! MPI variables
integer,allocatable     :: ini(:),fin(:)           ! Initial index; final index
integer                 :: MASTER
integer                 :: partxproc               ! Number of particles per processor
integer                 :: ierror,iproc,rank
integer                 :: request,numproc
integer                 :: stat(MPI_STATUS_SIZE)



! MPI initialization
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


MASTER=0

! Read coordinates of the trajectory (restart xyz file)
  open(unit=123,file='restart.rst',status='unknown',action='read')
    do i=1,snaps*part
      read(123,*), input(i,1), input(i,2), input(i,3)
    end do
  close(123)

! Distribution of the work per each processor
  allocate(ini(0:numproc),fin(0:numproc))
  partxproc=nint(real(part)/real(numproc))
  do i=0,numproc-2
    ini(i)=i*partxproc+1
    fin(i)=(i+1)*partxproc
  end do
  ini(numproc-1)=(numproc-1)*partxproc+1
  fin(numproc-1)=part


! Calculate the rmsd for each particle
call postvisual(input,snaps,dimen,part,desplacement,ini(rank),fin(rank),part)
!desplacement=1.

! Send the calculated distances to the MASTER
if (rank/=MASTER) then
print *, "Mando a Master", rank
  a=ini(rank)
  b=fin(rank)
print *, "A: ", a
print *, "B: ", b
  call MPI_ISEND(desplacement(:,a:b),snaps*(b-a+1),MPI_DOUBLE_PRECISION,MASTER,2,&
MPI_COMM_WORLD,request,ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD,ierror) !
call MPI_BARRIER(MPI_COMM_WORLD,ierror) ! ESPERATE COÑIO YA
call MPI_BARRIER(MPI_COMM_WORLD,ierror) ! 

if ( rank == MASTER ) then
print *, "Master recibe"
  do ip=1,numproc-1
    a=ini(ip)
    b=fin(ip)
    call MPI_RECV(desplacement(:,a:b),(snaps*(b-a+1)),MPI_DOUBLE_PRECISION,ip,2,&
MPI_COMM_WORLD,stat,ierror)
  end do

! Write output file
  open(unit=124,file='desplazamientos.dat',status='unknown',action='write')
    do i=1,snaps
      write(124,*), desplacement(i,:)
    end do
  close(124)
end if
call MPI_BARRIER(MPI_COMM_WORLD,ierror)



! MPI finalization
call MPI_FINALIZE(ierror)





contains

subroutine postvisual(input,snaps,dimen,part,desplacement,from_i,to_j,n)
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
      desplacement(i,j)=dist(input(j,:),input(frame,:))
    end do
  end do

end subroutine postvisual


function dist(x,y)
!!!
!!! Computes the distance between two points.
!!!
  real(8),dimension(dimen)     :: x,y    ! Input vectors
  real(8),dimension(dimen)     :: vector ! Output vectors
  real(8)                      :: dist

  vector=x-y

  dist=vector(1)**2+vector(2)**2+vector(3)**2

end function dist


end program postvis
