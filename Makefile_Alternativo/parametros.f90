module paralelizar
include 'mpif.h'
integer, public:: ierror,rank,numproc, request, MASTER=0
integer, dimension(MPI_STATUS_SIZE), public :: stat
integer, dimension(:),allocatable, public:: ini, fin
end module paralelizar