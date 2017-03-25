program initilize

implicit none

real*8                                 :: L
real*8, dimension(:,:), allocatable    :: pos
integer                                :: i, N

N = 1000
L = 40
allocate(pos(N,3))

call initialize(L,pos)

open(unit=123, file="forces_positions.xyz", status="replace", action="write")

write(unit=123,fmt='(i10,f20.10)') N, L
do i = 1, N
        write(unit=123, fmt='(3f20.10)') pos(i,:)
end do

close(unit=123)
contains
subroutine initialize(L,mat)
real*8, intent(in)					:: L
real*8, dimension(:,:), intent(out)	:: mat
integer								:: n, i
real*8								:: min_distance

n = size(mat,1)
mat = 0.0
do i = 1, n
    min_distance = -1.0 ! Any number smaller than T
    do while (min_distance < 1.0)
		call random_number(mat(i,:))
		mat(i,:) = L*mat(i,:) - L/2.0d0
        min_distance = distancemin(mat,i)
    enddo
enddo
end subroutine

function distancemin(mat,i) result(em)

! This function calculates the distances between all atoms and
! returns the value of the smaller distance.

real*8, dimension(:,:), intent(out)  :: mat
integer, intent(in)  				 :: i
integer                				 :: j
real*8, dimension(3)            	 :: vec
real*8                 				 :: em, distance

em = 1000.0
do j = 1, i

    if (i /= j) then

        vec = mat(i,:) - mat(j,:)
        distance = sqrt(sum(vec**2))
        if (distance < em) then
            em = distance
        endif

    endif

enddo

end function
end program
