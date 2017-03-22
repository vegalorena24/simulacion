!Integración de las posiciones y las velocidades con Velocity Verlet

!Xavier Marugan


!Se tiene que llamar al principio o al final de cada paso de la MD
!Incluye refold positions

subroutine IntegrationVerletPositions(pos,vel,forces,deltat,npart,mass,dim,lenght)

include "mpif.h"
use paralelizar

real,dimension(npart,3),intent(inout)::pos
real,dimension(npart,3),intent(in)::forces,vel
real,intent(in)::deltat,mass,lenght
integer,intent(in)::npart,dim
integer::i,ip
real::partxproc


! START CALCULATION
! Each worker computes a pice of translation of coordinates
do i=ini(rank),fin(rank)
   pos(i,:) = pos(i,:) + deltat*vel(i,:) + 0.5*(deltat**2)*forces(i,:)/mass      ! r(t+dt)
end do
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
    call MPI_RECV(pos(:,:),3*npart,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
end if

!Condiciones periódicas de contorno
call Refold_Positions(pos,npart,dim,lenght)

end subroutine IntegrationVerletPositions



!Se tiene que llamar 2 veces, antes y después de calcular las fuerzas
!porque tiene en cuenta la aceleración del tiempo anterior y la actual

subroutine IntegrationVerletVelocities(vel,forces,deltat,npart,mass)

include "mpif.h"
use paralelizar

real,dimension(npart,3),intent(inout)::vel
real,dimension(npart,3),intent(in)::forces
real,intent(in)::deltat,mass
integer,intent(in)::npart
integer::i

! START CALCULATION
do i=ini(rank),fin(rank)
   vel(i,:) = vel(i,:) + 0.5*deltat*forces(i,:)/mass  ! v(t+dt/2)
end do
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
    call MPI_RECV(vel(:,:),3*npart,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
end if


end subroutine IntegrationVerletVelocities
                                                         
