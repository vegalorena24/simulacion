
program m

!**********************MODULOS************************
use paralelizar
use initialization
use forces_mod
use integrators
use integrator_verlet
use postvisualization  
use analisis
implicit none
!*****************************************************



!*************DEFINICION VARIABLES********************
!Variables de MPI definidas en parametros.f90
real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc, density
integer:: N,dimnsion,Nsteps,i,j,step
integer:: Nrestart,frame,integrator
real(8) :: temperatura,presion
real*8, dimension(:,:), allocatable  :: positions,accel,velocities
logical                              :: randomize_initial_positions
logical                              :: random_initial_velocities
character(len=50)                    :: parameter_name
!*****************************************************



!****************INICIALIZACION MPI*******************
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
!*****************************************************




!******************INPUT******************************
open(100,file="input_data")
 
    read(100,*) parameter_name, deltat
    read(100,*) parameter_name, Nsteps
    read(100,*) parameter_name, N
    read(100,*) parameter_name, dimnsion
    read(100,*) parameter_name, BoxSize
    read(100,*) parameter_name, mass
    read(100,*) parameter_name, rc
    read(100,*) parameter_name, Nrestart
    read(100,*) parameter_name, integrator 
    read(100,*) parameter_name, randomize_initial_positions
    read(100,*) parameter_name, random_initial_velocities
    
close(100)
 
 
density=dble(N)/(Boxsize**3)
frame=Nsteps/Nrestart


allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))
!******************************************************




!******************INICIALIZACION************************
call initialize_system(N,density,positions,velocities,randomize_initial_positions,random_initial_velocities)
!********************************************************




!*********DISTRIBUCION DE PARTICULAS EN WORKERS**********
!ini y fin definidos en parametros.f90
allocate(ini(0:numproc-1),fin(0:numproc-1))
partxproc=nint(real(N)/real(numproc))
do i=0,numproc-2
        ini(i)=i*partxproc+1
        fin(i)=(i+1)*partxproc
end do
ini(numproc-1)=(numproc-1)*partxproc+1
fin(numproc-1)=N
!********************************************************



if (rank==MASTER) then
  open(unit=999,file='restart.rst',status='unknown',action='write')
end if
!**************MAIN LOOP*********************************
if (rank==MASTER) print*, "COMPUTING MAIN LOOP OF MOLECULAR DYNAMICS"

!Verlet integrator
  if (integrator==1) then
    call forces(positions,BoxSize,ini,fin,accel)
    do i = 1, Nsteps  
      call IntegrationVerletPositions(positions,velocities,accel,deltat,N,mass,dimnsion,BoxSize,ini,fin)
      call Refold_Positions(positions,N,dimnsion,BoxSize,ini,fin)
      call IntegrationVerletVelocities(velocities,accel,deltat,N,mass,ini,fin)
      call forces(positions,BoxSize,ini,fin,accel)
      call IntegrationVerletVelocities(velocities,accel,deltat,N,mass,ini,fin)
     temperatura = T_compute_paralel(N,velocities)
     presion = P_compute_paralel(N,BoxSize,positions,accel,temperatura)
      if (rank==MASTER) then
        if (mod(i,Nrestart)==0) then
          do j=1,N
            write(999,*), positions(i,:)
          end do
        end if
      end if
    end do
!Euler integrator
  else if (integrator==0) then    
    do i = 1, Nsteps
      call forces(positions,BoxSize,ini,fin,accel)
      call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat,ini,fin)
      call Refold_Positions(positions,N,dimnsion,BoxSize,ini,fin)
      call EulerVelocities(velocities,accel,N,dimnsion,BoxSize,mass,deltat,ini, fin)
     temperatura = T_compute_paralel(N,velocities)
     presion = P_compute_paralel(N,BoxSize,positions,accel,temperatura)
      if (rank==MASTER) then
        if (mod(i,Nrestart)==0) then
          do j=1,N
            write(999,*), positions(i,:)
          end do
        end if
      end if
    enddo
  end if
if (rank==MASTER) print*, "MAIN LOOP DONE"
!********************************************************
if (rank==MASTER) then
  close(999)
end if


!******************POSTVISUALIZACIÓN*********************
call postvisual(frame,dimnsion,N,ini,fin)
!********************************************************




!************FINALIZACION MPI****************************
call MPI_FINALIZE(ierror)
!********************************************************
end program m
