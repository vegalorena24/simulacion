! ********************************************************************************
!           Integration of the equations of motion : Euler algorithm
! ********************************************************************************

!                                                  Sergi Roca Bonet, Barcelona 2017

subroutine IntegrationEuler(positions, velocities, forces, deltat, particles, &
	                        mass, dim, length)
	real, intent(in)     ::  deltat, mass, length
	integer, intent(in)  ::  particles, dim
	integer              ::  k
	real, dimension(particles,3)  ::  positions, velocities, forces

	! Integrate positions and velocities
	do k = 1, particles
		positions(k,:) = positions(k,:) + ( velocities(k,:) + forces(k,:) * deltat ) * deltat
		velocities(k,:) = velocities(k,:) + ( forces(k,:) / mass ) * deltat
	enddo

	! Apply the PBC : Reput particles inside the box
    call Refold_Positions(positions,particles,dim,length)

end subroutine IntegrationEuler
