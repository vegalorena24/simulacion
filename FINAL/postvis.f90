module postvisualization

use paralelizar

contains

subroutine postvisual(snaps,dimen,part,ini,fin)
! Routine variables
  integer,intent(in)      :: snaps,dimen             ! Number of snapshots; dimensions 
  integer,intent(in)      :: part                    ! Number of particles            
  integer,intent(in)      :: ini(0:),fin(0:)         ! Rank index arrays            
  real(8)                 :: input(snaps*part,dimen) ! Input data    
  real(8)                 :: desplacement(snaps,part)   ! Dummy variable
  integer                 :: i,j,ip                  ! Iterators
  integer                 :: a,b                     ! Index variable

! Read coordinates of the trajectory (restart xyz file)
  open(unit=126,file='restart.rst',status='unknown',action='read')
    do i=1,snaps*part
      read(126,*), input(i,1), input(i,2), input(i,3)
    end do
  close(126)

! Calculate the rmsd for each particle
  call desplace(input,snaps,dimen,part,desplacement,ini(rank),fin(rank))
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! Send the calculated distances to the MASTER
  if (rank /= MASTER) then
    a=ini(rank)
    b=fin(rank)
    call MPI_ISEND(desplacement(:,a:b),snaps*(b-a+1),MPI_DOUBLE_PRECISION,MASTER,2,&
MPI_COMM_WORLD,request,ierror)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  if ( rank == MASTER ) then
    do ip=1,numproc-1
      a=ini(ip)
      b=fin(ip)
      call MPI_RECV(desplacement(:,a:b),(snaps*(b-a+1)),MPI_DOUBLE_PRECISION,ip,2,&
MPI_COMM_WORLD,stat,ierror)
    end do

! Write output file
    open(unit=125,file='desplazamientos.dat',status='unknown',action='write')
      do i=1,snaps
        write(125,*), desplacement(i,:)
      end do
    close(125)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

end subroutine postvisual


subroutine desplace(input,snaps,dimen,part,desplacement,from_i,to_j)
!!!
!!! Computes the RMSD for each particle.
!!!
  integer,intent(in)      :: from_i, to_j               ! Initial particle; final particle
  integer,intent(in)      :: snaps,dimen,part           ! Number of snapshots; dimensions; number of particles
  real(8),intent(in)      :: input(snaps*part,dimen)    ! Coordinates array
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


end module postvisualization

! ======================================================================================
