
program m
use paralelizar
use s
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
density=dble(N)/(Boxsize**3)

allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))










call initialize_system(N,density,positions,velocities)


if ( rank == MASTER ) then
    print *, 'Initialised'
endif



call MPI_FINALIZE(ierror)

end program m
