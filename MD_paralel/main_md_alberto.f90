module paralelizar

include 'mpif.h'
integer, public:: ierror,rank,numproc, request, stat, MASTER
integer, dimension(:),allocatable, public:: ini, fin
end module paralelizar


program m
use paralelizar
use forces_mod
implicit none

real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc
integer:: N,dimnsion,Nsteps,i,j,step
integer:: Nrestart,frame
real*8, dimension(:,:), allocatable:: positions,accel,velocities
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
!datos de entrada de prueba
deltat=0.001
Nsteps=100
N=4 
dimnsion=3
BoxSize=6.1984d0
mass=1.0d0
rc=2.5d0

 !!! Post-Vis
Nrestart=50
frame=Nsteps/Nrestart
 !!! Post-Vis

open (unit=10, File='coordenadas.dat')
!read(10,*) N, BoxSize
allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))

do i=1,N
 read(10,*) positions(i,:)
end do
close (10)

if ( rank == MASTER ) then
    print *, 'Initialised'
endif

velocities = 0.0d0

!MAIN

!call init()

allocate (ini(0:numproc-1),fin(0:numproc-1))
   !distribucion de particulas
    partxproc=nint(real(N)/real(numproc))
    do j=0,numproc-2
       ini(j)=j*partxproc+1
       fin(j)=(j+1)*partxproc
    end do
    ini(numproc-1)=(numproc-1)*partxproc+1
    fin(numproc-1)=N

open(unit=123,file='energy.dat',status='replace',action='write')
 !!! Post-Vis
open(unit=124,file='restart.rst',status='unknown',action='write')
 !!! Post-Vis

!accel = 10.0

do step=1,Nsteps

 call forces(positions,BoxSize,ini,fin,accel)

 call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)

 if ( rank == MASTER ) then
     print *, '******************************************'
     print *, '  step : ', step
     print *, 'r = ', positions(1,:)
     print *, 'v = ', velocities(1,:)
     print *, 'F = ', accel(1,:)
  endif

 !!! Post-Vis
  if (mod(step,Nrestart)==0) then
    do i=1,N
      write(124,*), positions(i,:)
    end do
  end if
 !!! Post-Vis
enddo
 !!! Post-Vis
close(124)
call postvisual(frame,dimnsion,N)
 !!! Post-Vis

call MPI_FINALIZE(ierror)

contains


! ======= EULER POSITIONS ==============================================================================
    subroutine EulerPositions(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    !use paralelizar
    !include 'mpif.h'
    integer::dimnsion,N,i,request,MASTER=0,iproc,partxproc
    real*8:: BoxSize, deltat, mass
    real*8,dimension(N,dimnsion):: pos, vel, forces
    double precision:: start_time, lapso_time,end_time

    ! Update of the positions using Euler
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
     pos(i,:) = pos(i,:) + ( vel(i,:) + deltat*forces(i,:)/mass ) * deltat
    end do
    lapso_time=MPI_Wtime()

    ! Send the positions from workers to Master
    if (rank /= MASTER) then
    call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Wait for all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receives and merge all the calculations
    if ( rank == MASTER ) then
     do iproc=1,numproc-1
      call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
     end do

   ! UPDATE
    do iproc=1,numproc-1
     call MPI_ISEND(pos(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
    end do
    end if

    ! Synchronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! All workers receive the positions once updated and merged
    if (rank /= MASTER ) then
      call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if
    end_time=MPI_Wtime()

    ! Apply periodic boundary conditions
    call Refold_Positions(pos,N,dimnsion,BoxSize)

    end subroutine EulerPositions

! ==== EULER VelocITIES =================================================================
    subroutine EulerVelocities(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    !use paralelizar
    !include 'mpif.h'
    integer::dimnsion,N,i,request,MASTER=0,iproc,partxproc !N=Number of part.
    real*8:: BoxSize, deltat, mass
    real*8,dimension(N,dimnsion):: pos, vel, forces !positions
    double precision:: start_time, lapso_time,end_time

    ! Update of the velocities using Euler
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
     vel(i,:) = vel(i,:) + deltat*forces(i,:)/mass
    end do
    lapso_time=MPI_Wtime()

    ! Send velocities from the workers to the Main
    if (rank /= MASTER) then
    call MPI_ISEND(vel(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Waiting all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receive and merge
    if ( rank == MASTER ) then
     do iproc=1,numproc-1
      call MPI_RECV(vel(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
     end do

    ! UPDATE
    do iproc=1,numproc-1
     call MPI_ISEND(vel(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
    end do
    end if

    ! Synchronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! All workers receive the velocities once updated
    if (rank /= MASTER ) then
      call MPI_RECV(vel(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if
    end_time=MPI_Wtime()

    end subroutine EulerVelocities

! ======================================================================================
    subroutine Refold_Positions(pos,N,dimnsion,BoxSize)

    !use paralelizar
    !include 'mpif.h'
    integer::dimnsion,N,i,request,MASTER=0,iproc,partxproc !N=Number of part.
    real*8:: BoxSize
    real*8,dimension(N,dimnsion):: pos !positions
    double precision:: start_time, lapso_time,end_time

    !start calculation
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
     pos(i,:) = pos(i,:) - BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
    end do
    lapso_time=MPI_Wtime()

    !sending work of worker(i) to master
    if (rank /= MASTER) then
    call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    !waiting all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    !master receive and merge
    if ( rank == MASTER ) then
     do iproc=1,numproc-1
      call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
     end do

   !UPDATE
    do iproc=1,numproc-1
     call MPI_ISEND(pos(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
    end do
    end if

    !Sincronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !all workers receive the coord
    if (rank /= MASTER ) then
      call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if
    end_time=MPI_Wtime()

    end subroutine Refold_Positions
! ======================================================================================


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
call desplace(input,snaps,dimen,part,desplacement,ini(rank),fin(rank),part)
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


! ======================================================================================

end program m
