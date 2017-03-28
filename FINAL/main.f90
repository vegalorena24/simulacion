
program m

!**********************MODULOS************************
use paralelizar
use s
use forces_mod
use integrators
implicit none
!*****************************************************



!*************DEFINICION VARIABLES********************
!Variables de MPI definidas en parametros.f90
real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc, density
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
!deltat=0.001
!Nsteps=500
!N=18
!dimnsion=3
!BoxSize=6.1984d0
!mass=1.0d0
!rc=2.5d0
!density=dble(N)/(Boxsize**3)

!datos de entrada de prueba
open(100,file="input_data")
 
    read(100,*) parameter_name, deltat
    read(100,*) parameter_name, Nsteps
    read(100,*) parameter_name, N
    read(100,*) parameter_name, dimnsion
    read(100,*) parameter_name, BoxSize
    read(100,*) parameter_name, mass
    read(100,*) parameter_name, rc
 
close(100)
 
 
density=dble(N)/(Boxsize**3)


allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))
!******************************************************




!******************INICIALIZACION************************
call initialize_system(N,density,positions,velocities)
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




!**************MAIN LOOP*********************************
if (rank==MASTER) print*, "COMPUTING MAIN LOOP OF MOLECULAR DYNAMICS"
do i = 1, Nsteps
    call forces(positions,BoxSize,ini,fin,accel)
    call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat,ini,fin)
    call Refold_Positions(positions,N,dimnsion,BoxSize,ini,fin)
    call EulerVelocities(velocities,accel,N,dimnsion,BoxSize,mass,deltat,ini, fin)
enddo
if (rank==MASTER) print*, "MAIN LOOP DONE"
!********************************************************



!************FINALIZACION MPI****************************
call MPI_FINALIZE(ierror)
!********************************************************
end program m
