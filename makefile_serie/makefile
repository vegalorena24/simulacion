mainMD : intento_main.o pbc.o integrate_euler.o integrate_verlet.o f_xabi_2.o mod_postvis_serie.o initialization.o analisis.o
	gfortran -o mainMD intento_main.o pbc.o integrate_euler.o integrate_verlet.o f_xabi_2.o mod_postvis_serie.o initialization.o analisis.o
initialization.o : initialization.f90
	gfortran -c initialization.f90
analisis.o : analisis.f90
	gfortran -c analisis.f90
pbc.o : pbc.f90
	gfortran -c pbc.f90
integrate_euler.o : integrate_euler.f90
	gfortran -c integrate_euler.f90
integrate_verlet.o : integrate_verlet.f90
	gfortran -c integrate_verlet.f90
f_xabi_2.o : f_xabi_2.f90
	gfortran -c f_xabi_2.f90
mod_postvis_serie.o : mod_postvis_serie.f90
	gfortran -c mod_postvis_serie.f90
intento_main.o : intento_main.f90 f_xabi_2.o mod_postvis_serie.o analisis.o
	gfortran -c intento_main.f90
