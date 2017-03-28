
program m

!**********************MODULOS************************
use paralelizar
use s
use forces_mod
use integradores
use albert
use analisis
implicit none
!*****************************************************



!*************DEFINICION VARIABLES********************
!Variables de MPI definidas en parametros.f90
real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc, density, presion,temperatura
integer:: N,dimnsion,Nsteps,i,j,step
integer:: Nrestart,frame
real*8, dimension(:,:), allocatable:: positions,accel,velocities
!*****************************************************



!****************INICIALIZACION MPI*******************
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
!*****************************************************




!******************INPUT******************************
!Este debería de estar en un fichero
!datos de entrada de prueba
deltat=0.001
Nsteps=100
N=18
dimnsion=3
BoxSize=6.1984d0
mass=1.0d0
rc=2.5d0
density=dble(N)/(Boxsize**3)

allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))
!******************************************************




!******************INICIALIZACION************************
call initialize_system(N,density,positions,velocities)
!********************************************************
!!! Post-Vis
Nrestart=50
frame=Nsteps/Nrestart
 !!! Post-Vis



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
if ( rank == MASTER ) then
!!! Post-Vis
open(unit=124,file='restart.rst',status='unknown',action='write')
 !!! Post-Vis
end if



!**************MAIN LOOP*********************************
do i = 1, Nsteps
    call forces(positions,BoxSize,ini,fin,accel) 
  call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
  call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
   temperatura = T_compute_paralel(N,velocities)
 presion = P_compute_paralel(N,BoxSize,positions,accel,temperatura)
  print *, '******************************************'
     print *, '  step : ', i
     print *, 'r = ', positions(1,:)
     print *, 'v = ', velocities(1,:)
     print *, 'F = ', accel(1,:)
     print *, 'T = ', temperatura
     print *, 'P = ', presion

   !!! Post-Vis


  if (mod(i,Nrestart)==0) then
    do j=1,N
      write(124,*), positions(j,:)
    end do
  end if

enddo


!********************************************************
 !!! Post-Vis
close(124)
call postvisual(frame,dimnsion,N)
 !!! Post-Vis


!************FINALIZACION MPI****************************
call MPI_FINALIZE(ierror)
!********************************************************
end program m
