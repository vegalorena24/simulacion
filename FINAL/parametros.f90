module paralelizar

include 'mpif.h'
integer, public:: ierror,rank,numproc, request, MASTER=0, iproc
integer, dimension(:),allocatable, public:: ini, fin
integer, dimension(MPI_STATUS_SIZE), public :: stat
end module paralelizar
