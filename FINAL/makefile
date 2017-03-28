MAIN : parametros.o main.o sotiris.o
	mpif90 -o MAIN parametros.o main.o sotiris.o
sotiris.o : sotiris.f90 parametros.o
	mpif90 -c sotiris.f90
parametros.o : parametros.f90
	mpif90 -c parametros.f90
main.o : main.f90 parametros.o sotiris.o
	mpif90 -c main.f90
