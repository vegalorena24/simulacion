MAIN : parametros.o main_md_sotiris.o postv.o forces_mod.o sotiris.o EULER.o analisis_paralel.o
	mpif90 -o MAIN parametros.o main_md_sotiris.o postv.o forces_mod.o sotiris.o  EULER.o analisis_paralel.o
	
EULER.o : EULER.f90 parametros.o
	mpif90 -c EULER.f90
	
sotiris.o : sotiris.f90 parametros.o
	mpif90 -c sotiris.f90
	
forces_mod.o : forces_mod.f90 parametros.o
	mpif90 -c forces_mod.f90
	
postv.o : postv.f90 parametros.o
	mpif90 -c postv.f90
	
parametros.o : parametros.f90
	mpif90 -c parametros.f90
	

analisis_paralel.o : analisis_paralel.f90 parametros.o
	mpif90 -c analisis_paralel.f90
	

main_md_sotiris.o : main_md_sotiris.f90 forces_mod.o postv.o parametros.o sotiris.o EULER.o analisis_paralel.o
	mpif90 -c main_md_sotiris.f90
	
clean:
	rm *.mod
	rm *.o
