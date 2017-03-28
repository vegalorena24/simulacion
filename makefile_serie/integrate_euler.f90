! ********************************************************************************
!           Integration of the equations of motion : Euler algorithm
! ********************************************************************************

!                                                  Sergi Roca Bonet, Barcelona 2017

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
