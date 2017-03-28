program m
use xabi
use postvisualization
use analisis2
implicit none

real*8:: deltat, BoxSize, mass,rc,epot, ekin
real(8) :: temperatura,presion,density
integer:: N,dimnsion,Nsteps,i,j,step,Nrestart,frame,integrar
real*8, dimension(:,:), allocatable:: positions,accel,velocities

print*, "Para integrador Verlet pulsa 1, cualquier otro numero tomara como integrador Euler."
read(5,*) integrar

!datos de entrada de prueba
deltat=0.001
Nsteps=1000
N=200
dimnsion=3
BoxSize=6.1984
mass=1.0
rc=BoxSize
density=dble(N)/(BoxSize**3)

!!! Post-Vis
Nrestart=50
frame=Nsteps/Nrestart
 !!! Post-Vis

allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))
!open (unit=10, File='coordenadas.dat')

!!! Post-Vis
open(unit=124,file='restart.rst',status='unknown',action='write')
 !!! Post-Vis


!do i=1,N
! read(10,*) positions(i,:)
!end do
!close (10)

!velocities = 0.0

!MAIN

!call init() 


open(unit=123,file='energy.dat',status='replace',action='write')

call initialize_system(N,density,positions,velocities)

call forces(positions,BoxSize,accel,rc,epot)
do step=1,Nsteps
 if ( integrar .eq. 1 ) then
 call IntegrationVerletPositions(positions,velocities,accel,deltat,N,mass,dimnsion,BoxSize)
 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)
 end if
 call forces(positions,BoxSize,accel,rc,epot)
 if ( integrar .eq. 1 ) then
 call IntegrationVerletVelocities(velocities,accel,deltat,N,mass)
 else
 call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 temperatura = T_compute(N,velocities)
 presion = P_compute(N,BoxSize,positions,accel,temperatura)
 end if

!     print *, '******************************************'
!     print *, '  step : ', step
!     print *, 'r = ', positions(1,:)
!     print *, 'v = ', velocities(1,:)
!     print *, 'F = ', accel(1,:)
!     print *, 'T = ', temperatura
!     print *, 'P = ', presion


 !!! Post-Vis
  if (mod(step,Nrestart)==0) then
    do i=1,N
      write(124,*), positions(i,:)
    end do
  end if
 !!! Post-Vis
 ekin = sum(mass*sqrt(sum(velocities**2,2))/2.0)

 write(unit=123,fmt='(i10,3f20.10)') step, ekin+epot, ekin, epot

enddo
close(123)
 !!! Post-Vis
close(124)
call postvisual(frame,dimnsion,N)
 !!! Post-Vis


end program m
