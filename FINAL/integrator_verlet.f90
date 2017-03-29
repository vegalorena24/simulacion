Module integrator_verlet
use paralelizar
contains

subroutine IntegrationVerletPositions(pos,vel,acc,deltat,npart,mass,dim,lenght,ini,fin)
integer::npart,dim
integer, dimension(0:), intent(in)::ini, fin
real*8,dimension(npart,3),intent(inout)::pos
real*8,dimension(npart,3),intent(in)::acc,vel
real*8,intent(in)::deltat,mass,lenght
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

if (rank == MASTER) then
   print *,"Verlet positions Time : ",end_time-start_time," seconds"
endif

end subroutine IntegrationVerletPositions


! ======================================================================================
subroutine IntegrationVerletVelocities(vel,acc,deltat,npart,mass,ini,fin)

integer::npart
integer, dimension(0:), intent(in)::ini, fin
real*8,dimension(npart,3),intent(inout)::vel
real*8,dimension(npart,3),intent(in)::acc
real*8,intent(in)::deltat,mass
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

if (rank == MASTER) then
   print *,"Verlet Velocities Time : ",end_time-start_time," seconds"
endif

end subroutine IntegrationVerletVelocities

end module
