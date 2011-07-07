#!/usr/bin/perl -w

# Example from Cambridge HPCS machine.
#  Sets 2 MPI processes (chains) per node (see value of $chainspn)
#  each generating 2 OpenMP threads (see value of $omp).
#  E.g. perl runMPI_HPCS.pl params 3
#  should run 6 chains, 2 per the three nodes (12 threads in total).

use Cwd;

# Use current directory as root
$cosmomc = cwd;

$params =  $ARGV[0];
$numnodes = $ARGV[1];

$chainspn = 2;  # number of chains per node
$omp = 2;       # value of OMP_NUM_THREADS
$ppn = ( $chainspn * $omp ) ; # NB this must be <= 4
$nchains = $numnodes * $chainspn ;
$mem = ( 7972 * $numnodes ) ;  # MB 

$ini = $params;
if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/bash
#PBS -q woodcrest
#PBS -N cosmomc
#PBS -l nodes=$numnodes:ppn=$ppn,mem=${mem}mb,walltime=10:00:00
#PBS -m abe
#PBS -r y
##PBS -W x=NACCESSPOLICY:SINGLEJOB
cd $cosmomc

echo Running on host \`hostname\`
echo Time is \`date\`
echo Directory is \`pwd\`
echo PBS job ID is \$PBS_JOBID
echo This jobs runs on the following machines:
echo \`cat \$PBS_NODEFILE | uniq\`

#! Create a machine file for MPI
#cat \$PBS_NODEFILE | uniq > scripts/machine.file.\$PBS_JOBID

export OMP_NUM_THREADS=$omp
export IPATH_NO_CPUAFFINITY=1
time mpiexec -npernode $chainspn ./cosmomc $ini > ./scripts/$params.log

#time mpirun -np $nchains -machinefile scripts/machine.file.\$PBS_JOBID ./cosmomc $ini > ./scripts/$params.log

EMP
close(Fout);

chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
chdir("../");


