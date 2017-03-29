#!/bin/csh
##########################################
# SGE options and parameters
##########################################
# (1) Name of the job
#$ -N paralel_en_Menorqui
# (2) Requested resources
# Parallel Environment and number of cores
#$ -pe omp* 8
# Queue
#$ -q cerqt2.q
# Shell
#$ -S /bin/csh
# (3) Output files
#$ -cwd
#$ -o gaussian_-METHOD-_-ANGLE-.out
#$ -e gaussian_-METHOD-_-ANGLE-.err
# (4) Remove the first '#' of the following 2 lines if you want to receive an email when the job ends.
##$ -m e
##$ -M  oscarcapote1618@gmail.com

##########################################
# User environment.
##########################################
# Load the modules needed

source /etc/profile.d/modules.csh
#module load gaussian/g09b01
module load openmpi/1.4.2_intel-11.1.072
##########################################
# Copying files needed
##########################################
# We copy the inputs to the directory where the jobs will run

setenv old `pwd`


cp -r MD_paralel $TMPDIR
cp -r MD_serie $TMPDIR
cp time_compute.sh $TMPDIR

cd $TMPDIR

# Set some variables for gaussian

#setenv GAUSS_SCRDIR $TMPDIR

##########################################
# Run the job
##########################################
# We run gaussian g09
sh time_compute.sh 8 1500


##########################################
# Copy the results to our home directory
##########################################


cp *.dat $old
cp *.gp $old

