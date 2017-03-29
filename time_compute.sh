#!/bin/sh
##################################################################
#                  AUTOR: Oscar Fajardo Fontiveros               #
#                                                                #
#    Programa que calcula el tiempo de ejecucion total de un     #
#    programa en serie y en paralelo en funcion del numero       #
#         de procesadores y del numero de particulas             #
#                                                                #
#                                                                #
#                                                                #
#                                                                #
#                                                                #
#                                                                #
#----------------------------COMO USARLO-------------------------#
#                                                                #
# Al vuestro programa de testeo le teneis que poner antes de     #
# compilar, en el contains el comando                            #
#                      include 'get_int.inc'                     #
#                                                                #
#  Luego le debeis decir a vuestro programa que el número de     #
#  partículas lo inicialice de la siguiente manera:              #
#                                                                #
#    N =  get_int_arg(1,int(256,8))                              #
#                                                                #
#                                                                #
#  Si trabajais con integer(4) teneis que poner:                 #
#                                                                #
#   N =  int(get_int_arg(1,int(256,8)),4)                        #
#                                                                #
#                                                                #
# Luego teneis que aseguraros que el archivo get_int.inc esté    #
# junto vuestros programas (tanto en serie como en paralelo)     #
# y que vuestros ejecutables se llamen main                      #
#                                                                #
# Luego teneis que ir al directorio donde está el script que     #
# estais leyendo y ejecutarlo así:                               #
#                                                                #
#        sh time_compute.sh maxProcesadores maxParticles         #
#                                                                #
#  AVISO: Debido a la naturaleza del comando time, no sé si      #
#         este script funciona en el cluster                     #
#                                                                #
#################################################################
#la comanda es....: module load openmpi/1.4.2_intel-11.1.072
thisDir=$(pwd)
maxProcesadors=$1
maxParticles=$2
if [ $# -ne 2 ]; then
    echo $0: usage: $0 maxProcesadores maxParticles
    exit 1
fi
exec 4> plot_times.gp
echo "set xlabel 'Num particules'">&4
echo "set ylabel 'Temps(s)'">&4
echo "set title 'Temps de computació(# particules)'">&4
for NProcesadors in $(seq 2 $maxProcesadors)
do
	cd $thisDir
	
	if [ "$NProcesadors" -eq 2 ]
	then
		echo "plot '"time_$NProcesadors.dat"' u 1:2 t 'paral·lel amb "$NProcesadors" procs.',\\">&4
	else
		echo "'"time_$NProcesadors.dat"' u 1:2 t 'paral·lel amb "$NProcesadors" procs.',\\">&4
	fi
	exec 3> time_$NProcesadors.dat
	for NParticles in $(seq 1 $maxParticles)
	do

		cd $thisDir/MD_paralel

		time_paralel=$( { /usr/bin/time -f "%e" mpirun -np $NProcesadors ./main $NParticles; } 2>&1)

		echo -e $NParticles' \t '$time_paralel >&3 #'\t'$time_series >&3
	done 
	exec 3<&- 
done

echo "'time_series.dat' u 1:2 t 'serie'\\">&4
exec 3> time_series.dat
for NParticles in $(seq 1 $maxParticles)
do

	cd $thisDir/MD_serie

	time_series=$( { /usr/bin/time -f "%e"  ./main $NParticles; } 2>&1)

	echo -e $NParticles' \t '$time_series #'\t'$time_series >&3
done 

	cd $thisDir
echo 'pause(-1)'>&4
