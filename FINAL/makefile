MAIN : parametros.o main.o initialization.o forces_mod.o integrators.o integrator_verlet.o postvis.o analisis_paralel.o
	mpif90 -o MAIN parametros.o main.o initialization.o forces_mod.o integrators.o integrator_verlet.o postvis.o analisis_paralel.o
	
forces_mod.o : forces_mod.f90 parametros.o
	mpif90 -c forces_mod.f90
	
integrators.o : integrators.f90 parametros.o
	mpif90 -c integrators.f90
	
integrator_verlet.o : integrator_verlet.f90 parametros.o
	mpif90 -c integrator_verlet.f90
	
initialization.o : initialization.f90 parametros.o
	mpif90 -c initialization.f90
	
postvis.o : postvis.f90 parametros.o
	mpif90 -c postvis.f90
	
parametros.o : parametros.f90
	mpif90 -c parametros.f90
	

analisis_paralel.o : analisis_paralel.f90 parametros.o
	mpif90 -c analisis_paralel.f90
	
main.o : main.f90 parametros.o initialization.o forces_mod.o integrators.o integrator_verlet.o postvis.o analisis_paralel.o
	mpif90 -c main.f90

.PHOPNY: clean
clean:
	@rm *.o
	@rm *.mod

