module forces_mod

use paralelizar

public :: forces
private :: lj

contains

subroutine forces(pos,L,ini,fin,forces_mat)
real*8, dimension(:,:), intent(in)		:: pos
real*8, intent(in)						:: L
real*8, dimension(:,:), intent(out)		:: forces_mat
integer, dimension(0:), intent(in)		:: ini, fin
integer									:: i, N, j, iproc, MASTER
real*8, dimension(:,:,:), allocatable	:: interaction_mat

MASTER =0
N = size(pos,1)
allocate(interaction_mat(N,N,3))

!CALCULATION OF THE INTERACTIONS BETWEEN PARTICLES
do i = ini(rank), fin(rank)

	interaction_mat(i,:,:) = 0.0

    !Particles lower diagonal
    do j = mod(i,2) + 1 , i - 1 , 2
        call lj(pos(i,:),pos(j,:),L,interaction_mat(i,j,:))
		!print*, 'lower',i, j, interaction_mat(i,j,1)
    enddo

    !Particles upper diagonal
    do j = i + 2, N, 2
        call lj(pos(i,:),pos(j,:),L,interaction_mat(i,j,:))
		!print*, 'upper', i, j, interaction_mat(i,j,1)
    enddo

enddo

!SENDING DATA FROM WORKERS TO MASTER
if (rank /= MASTER) then
	call MPI_ISEND(interaction_mat(ini(rank):fin(rank),:,:),N*(fin(rank)-ini(rank)+1)*3,&
		MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)
endif

!WAITING
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!MASTER RECIEVE
if (rank == MASTER) then
	do iproc = 1, numproc - 1
		call MPI_RECV(interaction_mat(ini(iproc):fin(iproc),:,:),N*3*(fin(iproc)-ini(iproc)+1),&
    		MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,stat,ierror)
	enddo

	!print*, "interaction", interaction_mat(1,1,:), rank
!	do i = 1, N
!			print*, 'mat',interaction_mat(i,:,1)
!	enddo

	!SUMMING ALL THE FORCES
	do i = 1, N
		forces_mat(i,:) = sum(interaction_mat(i,:,:),1)
		!print*, i, interaction_mat(i,:,:)
	enddo

	do i = 1, N
		forces_mat(i,:) = forces_mat(i,:) - sum(interaction_mat(:,i,:),1)
		!print*, i, forces_mat(i,:)
	enddo
endif

!SENDING FORCES_MAC TO ALL WORKERS
if (rank == MASTER) then
	do iproc = 1, numproc - 1
		call MPI_ISEND(forces_mat(:,:),N*3,&
		MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,request,ierror)
	enddo
endif

!WAITING
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

!RECIEVE
if (rank /= MASTER) then
	call MPI_RECV(forces_mat(:,:),N*3,&
    		MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,stat,ierror)
endif
deallocate(interaction_mat)
end subroutine



subroutine lj(vec1,vec2,L,force_vec)
real*8, dimension(:), intent(in)    :: vec1, vec2
real*8, intent(in)                  :: L
real*8, dimension(:), intent(out)   :: force_vec
real*8, dimension(size(vec1))       :: dist
real*8                              :: rr

dist = vec1 - vec2
dist = dist - L*nint(dist/L)
rr = sqrt(sum(dist**2))


force_vec = 24.d0*(2.d0/rr**14-1.d0/rr**8)*dist
end subroutine






end module
