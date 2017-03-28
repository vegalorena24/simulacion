!Integración de las posiciones y las velocidades con Velocity Verlet

!Xavier Marugan Ferrer


!Se tiene que llamar al principio o al final de cada paso de la MD
!Incluye refold positions

subroutine IntegrationVerletPositions(pos,vel,forces,deltat,npart,mass,dim,lenght)

integer,intent(in)::npart,dim
real*8,dimension(npart,3),intent(inout)::pos
real*8,dimension(npart,3),intent(in)::forces,vel
real*8,intent(in)::deltat,mass,lenght
integer::i

   pos = pos + deltat*vel + 0.5*(deltat**2)*forces/mass      ! r(t+dt)

!Condiciones periódicas de contorno
call Refold_Positions(pos,npart,dim,lenght)
   
end subroutine IntegrationVerletPositions



!Se tiene que llamar 2 veces, antes y después de calcular las fuerzas
!porque tiene en cuenta la aceleración del tiempo anterior y la actual

subroutine IntegrationVerletVelocities(vel,forces,deltat,npart,mass)

integer,intent(in)::npart
real*8,dimension(npart,3),intent(inout)::vel
real*8,dimension(npart,3),intent(in)::forces
real*8,intent(in)::deltat,mass
integer::i

   vel = vel + 0.5*deltat*forces/mass  ! v(t+dt/2)

end subroutine IntegrationVerletVelocities