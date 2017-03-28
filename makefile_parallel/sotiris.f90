module s
use paralelizar
contains
subroutine initialize_system(particles,density,positions,velocities)
implicit none


integer                                   :: ierror, request, rank, numproc, MASTER, iproc
integer                                         :: stat(MPI_STATUS_SIZE)
integer                                         :: partxproc
double precision                                :: start_time, lapso_time1, lapso_time2, end_time

double precision                                :: density
integer                                         :: particles
double precision, dimension (:,:), allocatable  :: v_cm!positions, velocities,
double precision, dimension (:,:)               :: positions, velocities !v_cm
double precision, dimension (3)                 :: displacement, Vcm=0.0
double precision                                :: L, dr
integer                                         :: lin_dim, N, i, last
!character(len=50)                               :: parameter_name


            !RANDOM NUMBER VARIABLES
            INTEGER :: ii, nn, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed




!call MPI_INIT(ierror)
!call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierror)
!call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)


!open(100,file="input_data")

 !   read(100,*) parameter_name, particles        
 !   read(100,*) parameter_name, density         

!close(100)



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


!if (rank==MASTER) print*, "partxproc:",partxproc, "L:",L, "lin_dim:",lin_dim, "dr:",dr, "numproc:",numproc


! Defining dimensions
!---------!allocate (positions(particles,3),velocities(particles,3),v_cm(0:numproc-1,3))
allocate (v_cm(0:numproc-1,3))
v_cm=0.0


start_time=MPI_Wtime()

if (rank /= numproc-1) then 
  do i=1, partxproc
    
    N=partxproc*rank+i
    !print*, N, "from rank:", rank

    positions(N,1)=mod(N,lin_dim)
    if (positions(N,1)==0) positions(N,1)=lin_dim

    positions(N,2)=mod((N-1)/lin_dim+1,lin_dim)
    if (positions(N,2)==0) positions(N,2)=lin_dim

    positions(N,3)=((N-1)/(lin_dim**2)+1)


    call random_number(displacement)
    positions(N,:)=positions(N,:)*dr-L/2+(displacement-0.5)*dr*(3.0/4.0)
    call random_number(velocities(N,:))
    velocities(N,:)=velocities(N,:)-0.5
    v_cm(rank,:)=v_cm(rank,:)+velocities(N,:)

  end do
  lapso_time1=MPI_Wtime()
  if (rank /= MASTER) then
    call MPI_ISEND(positions(partxproc*rank+1:partxproc*(rank+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)

    call MPI_ISEND(velocities(partxproc*rank+1:partxproc*(rank+1),:),3*partxproc,&
    MPI_DOUBLE_PRECISION,MASTER,2,MPI_COMM_WORLD,request,ierror)

    call MPI_ISEND(v_cm(rank,:),3,MPI_DOUBLE_PRECISION,MASTER,3,MPI_COMM_WORLD,request,ierror)
  end if
end if



if (rank==numproc-1) then
last=particles - partxproc*(numproc-1)
  do i=1, last

    N=partxproc*rank+i
    !print*, N, "from rank:", rank

    positions(N,1)=mod(N,lin_dim)
    if (positions(N,1)==0) positions(N,1)=lin_dim

    positions(N,2)=mod((N-1)/lin_dim+1,lin_dim)
    if (positions(N,2)==0) positions(N,2)=lin_dim

    positions(N,3)=((N-1)/(lin_dim**2)+1)


    call random_number(displacement)
    positions(N,:)=positions(N,:)*dr-L/2+(displacement-0.5)*dr*(3.0/4.0)
    call random_number(velocities(N,:))
    velocities(N,:)=velocities(N,:)-0.5
    v_cm(rank,:)=v_cm(rank,:)+velocities(N,:)

  end do
  lapso_time2=MPI_Wtime()

  call MPI_ISEND(positions(partxproc*rank+1:partxproc*rank+last,:),3*last,&
  MPI_DOUBLE_PRECISION,MASTER,1,MPI_COMM_WORLD,request,ierror)

  call MPI_ISEND(velocities(partxproc*rank+1:partxproc*rank+last,:),3*last,&
  MPI_DOUBLE_PRECISION,MASTER,2,MPI_COMM_WORLD,request,ierror)
  
  call MPI_ISEND(v_cm(rank,:),3,MPI_DOUBLE_PRECISION,MASTER,3,MPI_COMM_WORLD,request,ierror)
end if



! Waiting all WORKERS
call MPI_BARRIER(MPI_COMM_WORLD, ierror)



! MASTER receive and merge
if (rank == MASTER) then
  last=particles - partxproc*(numproc-1)

  do iproc=1,numproc-2
  !print*, "from processor:", iproc, "-------------------------"
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
end if  



end_time=MPI_Wtime()
if (rank/=numproc-1) then
print *,"rank: ",rank," time calculation: ",lapso_time1-start_time,"seconds"
end if
if (rank==numproc-1) then
print *,"rank: ",rank," time calculation: ",lapso_time2-start_time,"seconds"
end if
if (rank== MASTER) then
  print *,"MASTER. Updating information. Time: ",end_time-lapso_time1," seconds"
endif



if(rank==MASTER) then

do i=1,3
  Vcm(i)=sum(v_cm(:,i))
end do

!print*, Vcm/particles

! Writing final coordinates
open(unit=1, file='initial_positions.xyz')
open(unit=2, file='initial_velocities.xyz')
!--------!write(1,*) particles
!--------!write(1,*) 'iniital positions'
!--------!write(2,*) particles
!--------!write(2,*) 'initial velocities'
do i=1,particles
  write(1,*) 'H', positions(i,:) 
  write(2,*) 'H', velocities(i,:) - Vcm/particles
end do
close(1)
close(2)
end if


!! Waiting all WORKERS
!call MPI_BARRIER(MPI_COMM_WORLD, ierror)


end subroutine initialize_system

!======================================================================================================

 
end module
!======================================================================================================
