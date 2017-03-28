module paralelizar

include 'mpif.h'
integer, public:: ierror,rank,numproc, request, stat, MASTER
integer, dimension(:),allocatable, public:: ini, fin
end module paralelizar


program m
use paralelizar
use forces_mod
implicit none




real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc, density
integer:: N,dimnsion,Nsteps,i,j,step
integer:: Nrestart,frame
real*8, dimension(:,:), allocatable:: positions,accel,velocities
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
!datos de entrada de prueba
deltat=0.001
Nsteps=500
N=18
dimnsion=3
BoxSize=6.1984d0
mass=1.0d0
rc=2.5d0
density=N/(Boxsize**3)

allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))

call initialize_system(N,density,positions,velocities)

 !!! Post-Vis
Nrestart=50
frame=Nsteps/Nrestart
 !!! Post-Vis

!-------------------------------!open (unit=10, File='coordenadas.dat')
!-------------------------------!open (unit=10, File='initial_positions.xyz')
!-------------------------------!read(10,*) N, BoxSize

!-------------------------------!do i=1,N
!-------------------------------! read(10,*) positions(i,:)
!-------------------------------!end do
!-------------------------------!close (10)

if ( rank == MASTER ) then
    print *, 'Initialised'
endif

!-------------------------------!velocities = 0.0d0

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
call forces(positions,BoxSize,ini,fin,accel)

do step=1,Nsteps

 call IntegrationVerletPositions(positions,velocities,accel,deltat,N,mass,dimnsion,BoxSize)
 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)

 call forces(positions,BoxSize,ini,fin,accel)

 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)
 
 !call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 !call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)

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


print*,""
print*,""
print*,"                     FUCK YEAH!!!        "
print*,""
print*,""
print*,""
print*,""
print*,"      Sexy?Sex"
print*," ____?Sexy?Sexy"
print*," ___y?Sexy?Sexy?"
print*," ___?Sexy?Sexy?S"
print*," ___?Sexy?Sexy?S"
print*," __?Sexy?Sexy?Se"
print*," _?Sexy?Sexy?Se"
print*," _?Sexy?Sexy?Se"
print*," _?Sexy?Sexy?Sexy?"
print*," ?Sexy?Sexy?Sexy?Sexy"
print*," ?Sexy?Sexy?Sexy?Sexy?Se"
print*," ?Sexy?Sexy?Sexy?Sexy?Sex"
print*," _?Sexy?__?Sexy?Sexy?Sex"
print*," ___?Sex____?Sexy?Sexy?"
print*," ___?Sex_____?Sexy?Sexy"
print*," ___?Sex_____?Sexy?Sexy"
print*," ____?Sex____?Sexy?Sexy"
print*," _____?Se____?Sexy?Sex"
print*," ______?Se__?Sexy?Sexy"
print*," _______?Sexy?Sexy?Sex"
print*," ________?Sexy?Sexy?sex"
print*," _______?Sexy?Sexy?Sexy?Se"
print*," _______?Sexy?Sexy?Sexy?Sexy?"
print*," _______?Sexy?Sexy?Sexy?Sexy?Sexy"
print*," _______?Sexy?Sexy?Sexy?Sexy?Sexy?S"
print*," ________?Sexy?Sexy____?Sexy?Sexy?se"
print*," _________?Sexy?Se_______?Sexy?Sexy?"
print*," _________?Sexy?Se_____?Sexy?Sexy?"
print*," _________?Sexy?S____?Sexy?Sexy"
print*," _________?Sexy?S_?Sexy?Sexy"
print*," ________?Sexy?Sexy?Sexy"
print*," ________?Sexy?Sexy?S"
print*," ________?Sexy?Sexy"
print*," _______?Sexy?Se"
print*," _______?Sexy?"
print*," ______?Sexy?"
print*," ______?Sexy?"
print*," ______?Sexy?"
print*," ______?Sexy"
print*," ______?Sexy"
print*," _______?Sex"
print*," _______?Sex"
print*," _______?Sex"
print*," ______?Sexy"
print*," ______?Sexy"
print*," _______Sexy"
print*," _______ Sexy?"
print*," ________SexY"







contains

!==============================SUBRUTINA DE INICIALIZACION=======================================

subroutine initialize_system(particles,density,positions,velocities)
implicit none


!----!integer                                         :: ierror, request, rank, numproc, MASTER, iproc
integer                                         :: stat(MPI_STATUS_SIZE)
integer                                         :: partxproc,iproc
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
MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,request,ierror)
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
MPI_DOUBLE_PRECISION,iproc,1,MPI_COMM_WORLD,request,ierror)
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
end_time=MPI_Wtime()

end subroutine IntegrationVerletVelocities
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
    print *, "Leyendo"
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
print *, "Calculamos", rank
call desplace(input,snaps,dimen,part,desplacement,ini(rank),fin(rank),part)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! Send the calculated distances to the MASTER
if (rank /= MASTER) then
  print *, "Rank --> master"
  a=ini(rank)
  b=fin(rank)
  call MPI_ISEND(desplacement(:,a:b),snaps*(b-a+1),MPI_DOUBLE_PRECISION,MASTER,2,&
MPI_COMM_WORLD,request,ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD,ierror)

if ( rank == MASTER ) then
  do ip=1,numproc-1
    print *, "Master <-- rank", rank
    a=ini(ip)
    b=fin(ip)
    call MPI_RECV(desplacement(:,a:b),(snaps*(b-a+1)),MPI_DOUBLE_PRECISION,ip,2,&
MPI_COMM_WORLD,stat,ierror)
  end do

! Write output file
  open(unit=125,file='desplazamientos.dat',status='unknown',action='write')
    print *, "Imprimiendo"
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
