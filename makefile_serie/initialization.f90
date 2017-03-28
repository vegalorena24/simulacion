
! ======= INICIALIZACION  ===========================================================    ===================

subroutine initialize_system(particles,density,positions,velocities)
implicit none


double precision                                :: start_time, end_time
double precision                                :: density
integer                                         :: particles
double precision, dimension (particles,3)       :: positions, velocities
double precision, dimension (:,:), allocatable  :: displacement
double precision, dimension (3)                 :: Vcm=0.0
double precision                                :: L, dr
integer                                         :: lin_dim, N, i, x, y, z, cnt 
!------!character(len=50)                               :: parameter_name



            INTEGER :: ii, nn, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = nn) 
            ALLOCATE(seed(nn))
    
            CALL SYSTEM_CLOCK(COUNT=clock)
    
            seed = clock + 37 * (/ (ii - 1, ii = 1, nn) /)
            CALL RANDOM_SEED(PUT = seed)
    
            DEALLOCATE(seed)



!open(100,file="input_data")
!
!    read(100,*) parameter_name, particles
!    read(100,*) parameter_name, density      
!
!close(100)


L=(particles/density)**(1.0/3.0)
lin_dim=ceiling(particles**(1.0/3.0))
dr=L/(lin_dim+1)
cnt=1

! Defining dimensions
!------!allocate (positions(particles,3),velocities(particles,3),displacement(particles,3))
allocate (displacement(particles,3))


!------!call cpu_time(start_time)


do z=1,lin_dim
  do y=1,lin_dim
    do x=1,lin_dim
      if (cnt .gt. particles) exit !Avoid allocating unexisting particles
        positions(cnt,1) = dr*x-L/2
        positions(cnt,2) = dr*y-L/2
        positions(cnt,3) = dr*z-L/2
        cnt=cnt+1
    enddo
  enddo
enddo

call random_number(displacement)
!positions=positions-(displacement-0.5)*dr*(3.0/4.0)
call random_number(velocities)
velocities=0.0!velocities-0.5

do i=1,3
  Vcm(i)=sum(velocities(:,i))
  velocities(:,i)=velocities(:,i)-Vcm(i)/particles
end do


!------!call cpu_time(end_time)
!------!print*, "calculation time:", end_time-start_time


! Writing final coordinates
open(unit=1, file='initial_positions_serie.xyz')
open(unit=2, file='initial_velocities_serie.xyz')
write(1,*) particles
write(1,*) 'iniital positions'
write(2,*) particles
write(2,*) 'initial velocities'
do i=1,particles
  write(1,*) 'H', positions(i,:)
  write(2,*) 'H', velocities(i,:)
end do
close(1)
close(2)


end subroutine initialize_system
! =====================================================================================

