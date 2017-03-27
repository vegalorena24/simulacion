module postvisualization


contains
!!!
!!! 
!!!    call postvisual(snaps,dimen,part)
!!!
!!!

subroutine postvisual(snaps,dimen,part)


! Routine variables
integer,intent(in)      :: snaps,dimen             ! Number of snapshots; dimensions 
integer,intent(in)      :: part                    ! Number of particles            
real(8)                 :: input(snaps*part,dimen) ! Input data    
real(8)                 :: desplacement(snaps,part)   ! Dummy variable
integer                 :: i,j,ip                  ! Iterators
integer                 :: a,b                     ! Index variable



! Read coordinates of the trajectory (restart xyz file)
  open(unit=123,file='restart.rst',status='unknown',action='read')
    do i=1,snaps*part
      read(123,*), input(i,1), input(i,2), input(i,3)
    end do
  close(123)



! Calculate the rmsd for each particle
call desplace(input,snaps,dimen,part,desplacement,ini(rank),fin(rank),part)


! Write output file
  open(unit=124,file='desplazamientos.dat',status='unknown',action='write')
    do i=1,snaps
      write(124,*), desplacement(i,:)
    end do
  close(124)
end if


end subroutine postvisual


subroutine desplace(input,snaps,dimen,part,desplacement)
!!!
!!! Computes the RMSD for each particle.
!!!
  integer,intent(in)      :: snaps,dimen,part           ! Number of snapshots; dimensions; number of particles
  real(8),intent(in)      :: input(snaps*part,dimen)    ! Coordinates array
  real(8),intent(out)     :: desplacement(snaps,part)   ! Dummy variable
  integer                 :: frame                      ! Frame computed
  integer                 :: i,j                        ! Iterators

  desplacement=0.

  do j=1,part
    do i=1,snaps
      frame=j+((i-1)*part)
      desplacement(i,j)=dist(input(j,:),input(frame,:),dimen)
    end do
  end do

end subroutine desplace


function dist(x,y,dimen)
!!!
!!! Computes the distance between two points.
!!!
  integer                      :: dimen  ! dimension 
  real(8),dimension(dimen)     :: x,y    ! Input vectors
  real(8),dimension(dimen)     :: vector ! Output vectors
  real(8)                      :: dist

  vector=x-y

  dist=vector(1)**2+vector(2)**2+vector(3)**2

end function dist


end module postvisualization
