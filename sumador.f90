program hw

include 'mpif.h'
integer :: comm,rank,numproc,ierror,MASTER

integer(8),dimension(101) :: v
integer(8) :: suma_part
integer(8) :: suma_total
integer(8) :: interval_index
MASTER = 0

call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_Size(MPI_COMM_WORLD,numproc,ierror)
v = (/ (i , i=1,101) /)

interval_index = size(v)/numproc
if (rank/=numproc-1)then
    suma_part = sum(v((rank)*interval_index+1:(rank+1)*interval_index))
else
    suma_part = sum(v((rank)*interval_index+1:))
endif
call MPI_REDUCE(suma_part,suma_total,1,MPI_INTEGER,MPI_SUM,MASTER,MPI_COMM_WORLD,ierror)
print*,'Hola mon de processador',rank,'de',numproc,'suma',suma_part,'interval',(rank)*interval_index+1,(rank+1)*interval_index

if(rank==MASTER)then
print*,'MASTER diu, resultat:',suma_total
endif

call MPI_FINALIZE(ierror)
end program
