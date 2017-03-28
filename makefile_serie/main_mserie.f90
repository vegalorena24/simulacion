
program m

implicit none

real*8:: deltat, BoxSize, mass,rc,epot, ekin, density
integer:: N,dimnsion,Nsteps,i,j,step
real*8, dimension(:,:), allocatable:: positions,accel,velocities

!datos de entrada de prueba
deltat=0.001
Nsteps=100
N=18
dimnsion=3
BoxSize=6.1984
mass=1.0
rc=2.5
density=dble(N)/(BoxSize**3)

allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))
!--------------------------!open (unit=10, File='coordenadas.dat')
!--------------------------!do i=1,N
!--------------------------! read(10,*) positions(i,:)
!--------------------------!end do
!--------------------------!close (10)

!--------------------------!velocities = 0.0
print*, "hoola"
!MAIN

!call init() 

call initialize_system(N,density,positions,velocities)

open(unit=123,file='energy.dat',status='replace',action='write')

do step=1,Nsteps

 call forces(positions,BoxSize,accel,rc,epot)

 call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 call EulerVelocities(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)

 !call sample

! ekin = sum(mass*norm2(velocities,2)/2.0)

 write(unit=123,fmt='(i10,3f20.10)') step, ekin+epot, ekin, epot

enddo



contains


! ======= INICIALIZACION  ============================================================================

subroutine initialize_system(particles,density,positions,velocities)
implicit none


double precision                                :: start_time, end_time
double precision                                :: density
integer                                         :: particles
double precision, dimension (:,:)               :: positions, velocities
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
positions=positions-(displacement-0.5)*dr*(3.0/4.0)
call random_number(velocities)
velocities=velocities-0.5

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





subroutine forces(positions,boxlength,accel,rc,epot)
real*8, dimension(:,:), intent(in)  :: positions
real*8, dimension(:,:), intent(out) :: accel
real*8, intent(in)                  :: boxlength, rc
real*8, intent(out)                 :: epot
integer                             :: is, js, natoms
real*8                             :: pot

       natoms = size(positions,1)

        accel = 0.0d0
        epot = 0.d0

        !      atom-atom interactions

        do is = 1,natoms-1
           do js = is+1,natoms
              call lj(is,js,positions,boxlength,accel,rc,pot)
              epot = epot + pot
           end do
        end do


    end subroutine

    subroutine lj(is,js,positions,boxlength,accel,rc,pot)
    real*8, dimension(:,:), intent(in)      :: positions
    real*8,dimension(:,:), intent(inout)   :: accel
    real*8, intent(in)                      :: boxlength, rc
    integer, intent(in)                     :: is, js
    real*8, intent(out)                     :: pot
    real*8, dimension(size(positions,2))         :: rij
    real*8                                  :: rr2, rijl, rr, forcedist
    integer                                 :: l, dim

    dim = size(positions,2)

    rr2 = 0.d0
    pot = 0.d0
    do l = 1,dim
       rijl = positions(js,l) - positions(is,l)
       rij(l) = rijl - boxlength*nint(rijl/boxlength)
       rr2 = rr2 + rij(l)*rij(l)
    end do

    rr = sqrt(rr2)

    if (rr.lt.rc) then
       forcedist = 24*(2/rr**14-1.0/rr**8)
        pot = 4.0*(1.0/rr**12-1.0/rr**6)
        do l = 1,dim
            accel(is,l) = accel(is,l) - forcedist*rij(l)
            accel(js,l) = accel(js,l) + forcedist*rij(l)
        end do
    end if

    end subroutine

! ======= EULER POSITIONS ==============================================================================
    subroutine EulerPositions(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    
    integer::dimnsion,N,i !N=Number of part.
    real*8:: BoxSize, deltat, mass
    real*8,dimension(N,dimnsion):: pos, vel, forces !positions
   

    !start calculation
  
    do i=1,N
     pos(i,:) = pos(i,:) + ( vel(i,:) + deltat*forces(i,:)/mass ) * deltat
    end do

    call Refold_Positions(pos,N,dimnsion,BoxSize)

    end subroutine EulerPositions

! ==== EULER VELOCITIES =================================================================
    subroutine EulerVelocities(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat)

    integer::dimnsion,N,i !N=Number of part.
    real*8:: BoxSize, deltat, mass
    real*8,dimension(N,dimnsion):: pos, vel, forces !positions
  
    do i=1,N
     vel(i,:) = vel(i,:) + deltat*forces(i,:)/mass
    end do
  


    end subroutine EulerVelocities

! ======================================================================================
    subroutine Refold_Positions(pos,N,dimnsion,BoxSize)

    integer::dimnsion,N,i !N=Number of part.
    real*8:: BoxSize
    real*8,dimension(N,dimnsion):: pos !positions
 

    do i=1,N
     pos(i,:) = pos(i,:) - BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
    end do
   open (unit=11, File='coord_serie.dat')
    do i=1,N
     write(11,*) pos(i,:)
    end do
   close(11)
    end subroutine Refold_Positions
! ======================================================================================


end program m
