MAIN : parametros.o main.o sotiris.o forces_mod.o integrators.o postvis.o
	mpif90 -o MAIN parametros.o main.o sotiris.o forces_mod.o integrators.o postvis.o
forces_mod.o : forces_mod.f90 parametros.o
	mpif90 -c forces_mod.f90
integrators.o : integrators.f90 parametros.o
	mpif90 -c integrators.f90
sotiris.o : sotiris.f90 parametros.o
	mpif90 -c sotiris.f90
postvis.o : postvis.f90 parametros.o
	mpif90 -c postvis.f90
parametros.o : parametros.f90
	mpif90 -c parametros.f90
main.o : main.f90 parametros.o sotiris.o forces_mod.o integrators.o
	mpif90 -c main.f90
