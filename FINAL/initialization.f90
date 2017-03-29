module initialization
use paralelizar
contains
subroutine initialize_system(particles,density,positions,velocities,randomize_initial_positions,random_initial_velocities)
implicit none


integer                                         :: stat(MPI_STATUS_SIZE)
integer                                         :: partxproc, iproc
double precision                                :: start_time, lapso_time1, lapso_time2, end_time

double precision                                :: density
integer                                         :: particles
double precision, dimension (:,:), allocatable  :: v_cm
double precision, dimension (particles,3)       :: positions, velocities
double precision, dimension (3)                 :: displacement, Vcm=0.0
double precision                                :: L, dr
logical                                         :: randomize_initial_positions
logical                                         :: random_initial_velocities
integer                                         :: lin_dim, N, i, last


            !RANDOM NUMBER VARIABLES
            INTEGER :: ii, nn, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            !RANDOM NUMBER SEED GENERATION
            CALL RANDOM_SEED(size = nn)
            ALLOCATE(seed(nn))
    
            CALL SYSTEM_CLOCK(COUNT=clock)
    
            seed = clock + 37 * (/ (ii - 1, ii = 1, nn) /)
            CALL RANDOM_SEED(PUT = seed)
    
            DEALLOCATE(seed)






MASTER=0
partxproc=nint(real(particles)/real(numproc))

L=(particles/density)**(1.0/3.0)
lin_dim=ceiling(particles**(1.0/3.0))
dr=L/(lin_dim+1)


if (rank==MASTER) then
  print*
  print*, "particles per processor:",partxproc
  print*, "box length (L):",L
  print*, "max particles in each dimension:",lin_dim
  print*, "ininital inter-particle distance (dr):",dr
  print*, "number of processors:",numproc
  print*
end if

! Defining dimensions
allocate (v_cm(0:numproc-1,3))
v_cm=0.0


start_time=MPI_Wtime()

if (rank /= numproc-1) then 
  do i=1, partxproc
    
    N=partxproc*rank+i

    positions(N,1)=mod(N,lin_dim)
    if (positions(N,1)==0) positions(N,1)=lin_dim

    positions(N,2)=mod((N-1)/lin_dim+1,lin_dim)
    if (positions(N,2)==0) positions(N,2)=lin_dim

    positions(N,3)=((N-1)/(lin_dim**2)+1)


    if (randomize_initial_positions) then
      call random_number(displacement)
      positions(N,:)=positions(N,:)*dr-L/2+(displacement-0.5)*dr*(3.0/4.0)
    else
      positions(N,:)=positions(N,:)*dr-L/2
    end if

    if (random_initial_velocities) then
      call random_number(velocities(N,:))
      velocities(N,:)=(velocities(N,:)-0.5)*dr/6
    else
      velocities=0.0
    end if
    v_cm(rank,:)=v_cm(rank,:)+velocities(N,:)

  end do
  lapso_time1=MPI_Wtime()
  if (rank /= MASTER) then
    call MPI_SEND(positions(partxproc*rank+1:partxproc*(rank+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)

    call MPI_SEND(velocities(partxproc*rank+1:partxproc*(rank+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,MASTER,2,MPI_COMM_WORLD,request,ierror)

    call MPI_SEND(v_cm(rank,:),3,MPI_DOUBLE_PRECISION,MASTER,3,MPI_COMM_WORLD,request,ierror)
  end if
end if



if (rank==numproc-1) then
last=particles - partxproc*(numproc-1)
  do i=1, last

    N=partxproc*rank+i

    positions(N,1)=mod(N,lin_dim)
    if (positions(N,1)==0) positions(N,1)=lin_dim

    positions(N,2)=mod((N-1)/lin_dim+1,lin_dim)
    if (positions(N,2)==0) positions(N,2)=lin_dim

    positions(N,3)=((N-1)/(lin_dim**2)+1)

    if (randomize_initial_positions) then
      call random_number(displacement)
      positions(N,:)=positions(N,:)*dr-L/2+(displacement-0.5)*dr*(3.0/4.0)
    else
      positions(N,:)=positions(N,:)*dr-L/2
    end if

    if (random_initial_velocities) then
      call random_number(velocities(N,:))
      velocities(N,:)=(velocities(N,:)-0.5)*dr/6
    else
      velocities=0.0
    end if
    v_cm(rank,:)=v_cm(rank,:)+velocities(N,:)

  end do
  lapso_time2=MPI_Wtime()

  call MPI_SEND(positions(partxproc*rank+1:partxproc*rank+last,:),3*last,&
  MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)

  call MPI_SEND(velocities(partxproc*rank+1:partxproc*rank+last,:),3*last,&
  MPI_DOUBLE_PRECISION,MASTER,2,MPI_COMM_WORLD,request,ierror)
  
  call MPI_SEND(v_cm(rank,:),3,MPI_DOUBLE_PRECISION,MASTER,3,MPI_COMM_WORLD,request,ierror)
end if



!! Waiting all WORKERS
!call MPI_BARRIER(MPI_COMM_WORLD, ierror)



! MASTER receive and merge
if (rank == MASTER) then
  last=particles - partxproc*(numproc-1)

  do iproc=1,numproc-2
    call MPI_RECV(positions(partxproc*iproc+1:partxproc*(iproc+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,stat,ierror)
    call MPI_RECV(velocities(partxproc*iproc+1:partxproc*(iproc+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,iproc,2,MPI_COMM_WORLD,stat,ierror)
    call MPI_RECV(v_cm(iproc,:),3,MPI_DOUBLE_PRECISION,iproc,3,MPI_COMM_WORLD,stat,ierror)
  end do

  call MPI_RECV(positions(partxproc*(numproc-1)+1:partxproc*(numproc-1)+last,:),3*last,&
  MPI_DOUBLE_PRECISION,(numproc-1),1,MPI_COMM_WORLD,stat,ierror)
  call MPI_RECV(velocities(partxproc*(numproc-1)+1:partxproc*(numproc-1)+last,:),3*last,&
  MPI_DOUBLE_PRECISION,(numproc-1),2,MPI_COMM_WORLD,stat,ierror)
  call MPI_RECV(v_cm(numproc-1,:),3,MPI_DOUBLE_PRECISION,(numproc-1),3,MPI_COMM_WORLD,stat,ierror)

  do i=1,3
    Vcm(i)=sum(v_cm(:,i))
  end do

  do i=1,particles
    velocities(i,:)=velocities(i,:) - Vcm/particles
  end do

  do iproc=1,numproc-1
    call MPI_ISEND(positions(:,:),3*particles,MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
    call MPI_ISEND(velocities(:,:),3*particles,MPI_REAL8,iproc,1,MPI_COMM_WORLD,request,ierror)
  end do
end if  

! Sincronize all processors
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! all WORKERS receive coordinates
if (rank /= MASTER) then
    call MPI_RECV(positions(:,:),3*particles,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    call MPI_RECV(velocities(:,:),3*particles,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
end if


end_time=MPI_Wtime()
!if (rank/=numproc-1) then
!print *,"rank: ",rank," time calculation: ",lapso_time1-start_time,"seconds"
!end if
!if (rank==numproc-1) then
!print *,"rank: ",rank," time calculation: ",lapso_time2-start_time,"seconds"
!end if
!if (rank== MASTER) then
!  print *,"MASTER. Updating information. Time: ",end_time-lapso_time1," seconds"
!endif



if(rank==MASTER) then

! Writing final coordinates
open(unit=1, file='initial_positions.xyz')
open(unit=2, file='initial_velocities.xyz')
do i=1,particles
  write(1,*) positions(i,:) 
  write(2,*) velocities(i,:)
end do
close(1)
close(2)
end if


! Waiting all WORKERS
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

!do i=1,particles
!  print*, positions(i,:)
!end do

!print*
!do i=1,particles
!  print*, velocities(i,:)
!end do


end subroutine initialize_system

!======================================================================================================

 
end module initialization
!======================================================================================================
