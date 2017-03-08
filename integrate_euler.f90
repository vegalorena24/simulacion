! ********************************************************************************
!           Integration of the equations of motion : Euler algorithm
! ********************************************************************************

!                                                  Sergi Roca Bonet, Barcelona 2017

include "mpif.h"

subroutine IntegrationEuler(positions, velocities, forces, deltat, particles, &
	                        mass, dim, length)
	real, dimension(particles,3)  ::  positions, velocities, forces
	real, intent(in)     ::  deltat, mass, length
	integer, intent(in)  ::  particles, dim
	integer              ::  i, j, counter, ierror

! Define the communicator
	!comm = MPI_COMM_WORLD
! Find out the number of processes
	!call MPI_COMM_SIZE(comm, numproc, ierror)
! Obtain the taskid, the id of an individual process
	!call MPI_COMM_RANK(comm, taskid, ierror)

	counter = particles / numproc

! Integrate positions
	do k = 1, particles
		positions(k,1) = positions(k,1) + ( velocities(k,1) + forces(k,1)*deltat ) * deltat
		positions(k,2) = positions(k,2) + ( velocities(k,2) + forces(k,2)*deltat ) * deltat
		positions(k,3) = positions(k,3) + ( velocities(k,3) + forces(k,3)*deltat ) * deltat
	enddo

! Call PBC subroutine (Lorena)
	call Refold_Positions(positions, particles, dim, length)

! Integrate velocities
	do k = 1, particles
		velocities(k,1) = velocities(k,1) + ( forces(k,1) / mass ) * deltat
		velocities(k,2) = velocities(k,2) + ( forces(k,2) / mass ) * deltat
		velocities(k,3) = velocities(k,3) + ( forces(k,3) / mass ) * deltat
	enddo

end subroutine IntegrationEuler
