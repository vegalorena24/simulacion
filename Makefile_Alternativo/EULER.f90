module integradores
use paralelizar
contains
! ======= EULER POSITIONS ==============================================================================
    subroutine EulerPositions(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    
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
      call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
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
      call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    ! Apply periodic boundary conditions
    call Refold_Positions(pos,N,dimnsion,BoxSize)

    end subroutine EulerPositions

! ==== EULER VelocITIES =================================================================
    subroutine EulerVelocities(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    
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
      call MPI_RECV(vel(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
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
      call MPI_RECV(vel(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    end subroutine EulerVelocities

! ======================================================================================
 subroutine IntegrationVerletPositions(pos,vel,acc,deltat,npart,mass,dim,lenght)

!include "mpif.h"
!use paralelizar

real*8,dimension(npart,3),intent(inout)::pos
real*8,dimension(npart,3),intent(in)::acc,vel
real*8,intent(in)::deltat,mass,lenght
integer::npart,dim,MASTER=0,iproc,partxproc,request
integer::i,ip
double precision:: start_time, lapso_time,end_time

! START CALCULATION
start_time=MPI_Wtime()
! Each worker computes a pice of translation of coordinates
do i=ini(rank),fin(rank)
   pos(i,:) = pos(i,:) + deltat*vel(i,:) + 0.5*(deltat**2)*acc(i,:)/mass      ! r(t+dt)
end do
lapso_time=MPI_Wtime()
! END CALCULATION

! sending work of WORKER(i) to MASTER
if (rank /= MASTER) then
        call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),&
MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)
end if
! Waiting all WORKERS
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! MASTER receive and merge
if (rank == MASTER) then
  do iproc=1,numproc-1
    call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),&
MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,stat,ierror)
  end do

!!! UPDATE the coordinations of MASTER to WORKERS
  do iproc=1,numproc-1
    call MPI_ISEND(pos(:,:),3*npart,MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
  end do
end if

! Sincronize all processors
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! all WORKERS receive coordinates
if (rank /= MASTER) then
    call MPI_RECV(pos(:,:),3*npart,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
end if
end_time=MPI_Wtime()

!Condiciones periódicas de contorno
call Refold_Positions(pos,npart,dim,lenght)

end subroutine IntegrationVerletPositions

! ======================================================================================
subroutine IntegrationVerletVelocities(vel,acc,deltat,npart,mass)

!include "mpif.h"
!use paralelizar

real*8,dimension(npart,3),intent(inout)::vel
real*8,dimension(npart,3),intent(in)::acc
real*8,intent(in)::deltat,mass
integer::npart,MASTER=0,iproc,partxproc,request
double precision:: start_time, lapso_time,end_time

! START CALCULATION
start_time=MPI_Wtime()
do i=ini(rank),fin(rank)
   vel(i,:) = vel(i,:) + 0.5*deltat*acc(i,:)/mass  ! v(t+dt/2)
end do
lapso_time=MPI_Wtime()
! END CALCULATION

! sending work of WORKER(i) to MASTER
if (rank /= MASTER) then
        call MPI_ISEND(vel(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),&
MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)
end if
! Waiting all WORKERS
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! MASTER receive and merge
if (rank == MASTER) then
  do iproc=1,numproc-1
    call MPI_RECV(vel(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),&
MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,stat,ierror)
  end do

!!! UPDATE the coordinations of MASTER to WORKERS
  do iproc=1,numproc-1
    call MPI_ISEND(vel(:,:),3*npart,MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
  end do
end if

! Sincronize all processors
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! all WORKERS receive coordinates
if (rank /= MASTER) then
    call MPI_RECV(vel(:,:),3*npart,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
end if
end_time=MPI_Wtime()

end subroutine IntegrationVerletVelocities
! ======================================================================================
subroutine Refold_Positions(pos,N,dimnsion,BoxSize)


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
      call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
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
      call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    end subroutine Refold_Positions


end module