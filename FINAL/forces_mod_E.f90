module forces_mod
use paralelizar

public:: forces
private:: lj 

contains

subroutine forces(pos,L,ini,fin,forces_mat_sum)
real*8, dimension(:,:),intent(in):: pos
real*8, dimension(:,:),intent(out)::forces_mat_sum
integer, dimension(0:):: ini, fin
real*8,intent(in):: L
real*8,dimension(:,:,:),allocatable:: forces_mat
real*8 :: start_time,lapso_time,end_time,start_sum_time,end_sum_time,epot
integer:: particles,iproc,i,j,k,MASTER



particles = size(pos,1)
allocate(forces_mat(particles,particles,3))

MASTER=0

start_time=MPI_Wtime()


!start forces calculation
do i=ini(rank),fin(rank)
   forces_mat(i,:,:)=0.d0
   do j=1,particles
      if(i == j) then
         forces_mat(i,j,:) = 0.
      else
         call lj(i,j,pos,L,forces_mat,L,epot)
      end if
   end do
end do
!end forces calculation

lapso_time=MPI_Wtime()

!send work from all workers to MASTER

if(rank /= MASTER) then
    call MPI_SEND(forces_mat(ini(rank):fin(rank),:,:),3*(fin(rank)-ini(rank)+1)*particles, &
    MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)
end if


! MASTER receive and merge
if (rank == MASTER) then
   do iproc=1,numproc-1
      call MPI_RECV(forces_mat(ini(iproc):fin(iproc),:,:),3*(fin(iproc)-ini(iproc)+1)*particles,&
        MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,stat,ierror)
   end do


   start_sum_time = MPI_Wtime()
   forces_mat_sum = 0.
   !sum forces all particles do to one
   do i=1,particles
      do k=1,3
         do j=1,particles
            forces_mat_sum(i,k) = forces_mat_sum(i,k) + forces_mat(i,j,k)
         end do
      end do
   end do
   end_sum_time = MPI_Wtime()



   !!! UPDATE the coordinations of MASTER to WORKERS
   do iproc=1,numproc-1
      call MPI_ISEND(forces_mat_sum(:,:),3*particles,MPI_DOUBLE_PRECISION,iproc,1, &
            MPI_COMM_WORLD,request,ierror)
   end do
end if

! Sincronize all processors
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! all WORKERS receive coordinates
if (rank /= MASTER) then
   call MPI_RECV(forces_mat_sum(:,:),3*particles,MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,stat,ierror)
end if

end_time=MPI_Wtime()


!print *, "forces calculation time", rank, lapso_time-start_time
!print *, 'time to sum total force', end_sum_time-start_sum_time
!print *, "comunication time", rank, end_time-lapso_time



end subroutine




subroutine lj(is,js,pos,boxlength,accel,rc,pot)
real*8, dimension(:,:), intent(in):: pos
real*8, dimension(:,:,:), intent(inout):: accel
real*8, intent(in):: boxlength, rc
integer, intent(in):: is, js
real*8, intent(out):: pot
real*8, dimension(size(pos,2)):: rij
real*8:: rr2, rijl, rr, forcedist
integer:: l, dimen

!accel = 0.0d0

dimen = size(pos,2)

rr2 = 0.d0
pot = 0.d0
do l = 1,dimen
   rijl = pos(js,l) - pos(is,l)
   rij(l) = rijl - boxlength*nint(rijl/boxlength)
   rr2 = rr2 + rij(l)*rij(l)
end do

rr = sqrt(rr2)

if (rr.lt.rc) then
   forcedist = 24.d0*(2.d0/rr**14-1.0d0/rr**8)
   pot = 4.d0*(1.0d0/rr**12-1.0d0/rr**6)
   do l = 1,dimen
      accel(is,js,l) = accel(is,js,l) - forcedist*rij(l)
   end do
end if

end subroutine


end module
