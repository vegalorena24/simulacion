module paralelizar

include 'mpif.h'
integer, public:: ierror,rank,numproc, request, stat, MASTER=0
integer, dimension(:),allocatable, public:: ini, fin
end module paralelizar