!Integración de las posiciones y las velocidades con Velocity Verlet

!Xavier Marugan



!Se tiene que llamar al principio o al final de cada paso de la MD
!Incluye refold positions

subroutine IntegrationVerletPositions(pos,vel,forces,deltat,npart,mass,dim,lenght)

real,dimension(npart,3),intent(inout)::pos
real,dimension(npart,3),intent(in)::forces,vel
real,intent(in)::deltat,mass,lenght
integer,intent(in)::npart,dim
integer::i

do i=1,npart
   pos = pos + deltat*vel + 0.5*(deltat**2)*forces/mass      ! r(t+dt)
end do

!Condiciones periódicas de contorno
call refold_positions(pos,npart,dim,lenght)
   
end subroutine IntegrationVerletPositions



!Se tiene que llamar 2 veces, antes y después de calcular las fuerzas
!porque tiene en cuenta la aceleración del tiempo anterior y la actual

subroutine IntegrationVerletVelocities(vel,forces,deltat,npart,mass)

real,dimension(npart,3),intent(inout)::vel
real,dimension(npart,3),intent(in)::forces
real,intent(in)::deltat,mass
integer,intent(in)::npart
integer::i

do i=1,npart
   vel = vel + 0.5*deltat*forces/mass  ! v(t+dt/2)
end do

end subroutine IntegrationVerletVelocities
