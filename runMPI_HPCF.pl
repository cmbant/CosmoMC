#!/usr/local/bin/perl

#Example from Cambridge HPCF machine
#Sets 2 processes per node 
#(so quad-core each process can use 2 openmp threads if set in .ini or environment)
#perl runMPI_HPCF.pl params 3
#should run 6 chains, 2 per the three nodes.

#also set number of threads (=2 in example below) in params.ini
#or maybe just export OMP_NUM_THREADS=2 in ~/.bashrc. 

#Number of cores per node (property of cluster)
$ppn = 4;

#Number of chains to run per mode (each $ppm/$chainspn openmp theads)
$chainspn = 2;

use Cwd;

#Use current directory as root
$cosmomc = cwd;

$params =  $ARGV[0];
$numnodes = $ARGV[1];

$nchains = $numnodes*$chainspn;

$ini = $params;
if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

$path = $cosmomc;

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/csh -f
#PBS -N cosmomc
#PBS -l nodes=$numnodes:ppn=$ppn,walltime=10:00:00
#PBS -r n
##PBS -W x=NACCESSPOLICY:SINGLEJOB
cd $cosmomc

echo Running on host \`hostname\`
echo Time is \`date\`
echo Directory is \`pwd\`
echo PBS job ID is \$PBS_JOBID
echo This jobs runs on the following machines:
echo \`cat \$PBS_NODEFILE | uniq\`

#! Create a machine file for MPI
cat \$PBS_NODEFILE | uniq > scripts/machine.file.\$PBS_JOBID

time mpirun -np $nchains -machinefile scripts/machine.file.\$PBS_JOBID ./cosmomc $ini > ./scripts/$params.log

EMP
close(Fout);

chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
chdir("../");
