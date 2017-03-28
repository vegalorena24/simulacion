
program m
use paralelizar
use forces_mod
use albert
use s
use integradores
use analisis
implicit none




real*8 :: deltat, BoxSize, mass,rc,epot, ekin,partxproc, density
integer:: N,dimnsion,Nsteps,i,j,step
integer:: Nrestart,frame
real*8, dimension(:,:), allocatable:: positions,accel,velocities
real(8) :: temperatura,presion
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
density=dble(N)/(Boxsize**3)

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

do step=1,Nsteps

 call IntegrationVerletPositions(positions,velocities,accel,deltat,N,mass,dimnsion,BoxSize)
 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)

 call forces(positions,BoxSize,ini,fin,accel)

 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)
 temperatura = T_compute_paralel(N,velocities)
 presion = P_compute_paralel(N,BoxSize,positions,accel,temperatura)
 
 !call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 !call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)

 if ( rank == MASTER ) then
     print *, '******************************************'
     print *, '  step : ', step
     print *, 'r = ', positions(1,:)
     print *, 'v = ', velocities(1,:)
     print *, 'F = ', accel(1,:)
     print *, 'T = ', temperatura
     print *, 'P = ', presion
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

end program m
