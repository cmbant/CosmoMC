#!/usr/bin/perl -w

#  Attempting to change to SLURM...

# SG Now sets 4 MPI processes per node, 4 OMP threads per process

# Example from Cambridge HPCS machine.
#  Sets 4 MPI processes (chains) per node (see value of $chainspn)
#  each generating 4 OpenMP threads (see value of $omp).
#  E.g. perl runMPI_HPCS.pl params 2
#  should run 8 chains, 4 per the two nodes (32 threads in total).
# for action =2 runs in one mpi process (e.g. for minimization) use
# E.g. perl runMPI_HPCS.pl params 0

# Note: the runbatch.py script now uses the job_script template instead

use Cwd;

# Use current directory as root
$cosmomc = cwd;

$params =  $ARGV[0];
$numnodes = $ARGV[1];

#if $numnodes=0 just run one MPI process (e.g. for action=2) using all cores
if ($numnodes == 0){
$chainspn=1;
$omp=16;
$numnodes=1;
$walltime='08:00:00';
} else
{
$chainspn = 4;  # number of chains per node
$omp = 4;       # value of OMP_NUM_THREADS
$walltime='24:00:00';
}

$ppn = ( $chainspn * $omp ) ; # NB this must be <= 16
$nchains = $numnodes * $chainspn ;
$mem = ( 63900 * $numnodes ) ;  # MB 

$ini = $params;
if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

$params =~ s/\//_/g;

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/bash
#SBATCH -p sandybridge
#SBATCH -A PLANCK
#SBATCH --nodes=$numnodes
#SBATCH --ntasks=$nchains
#SBATCH --time=$walltime 
#SBATCH --mail-type=FAIL
#SBATCH --mem=63900
##PBS -N $params
##PBS -l nodes=$numnodes:ppn=$ppn,mem=${mem}mb,walltime=$walltime
##PBS -m n
##PBS -r n

cd $cosmomc

. /etc/profile.d/modules.sh
#module purge
module load default-impi
module load cfitsio

echo Running on host \`hostname\`
echo Time is \`date\`
echo Directory is \`pwd\`
echo Slurm job ID is \$SLURM_JOBID
echo This jobs runs on the following machines:
export NODEFILE=`generate_pbs_nodefile`
echo \`cat \$NODEFILE | uniq\`

#! Create a machine file for MPI
cat \$NODEFILE | uniq > scripts/machine.file.\$SLURM_JOBID

export OMP_NUM_THREADS=$omp
export I_MPI_PIN_DOMAIN=omp:compact
export I_MPI_PIN_ORDER=scatter
export I_MPI_CPUINFO=proc

time mpirun -ppn $chainspn -np $nchains ./cosmomc $ini > ./scripts/$params.log 2>&1

EMP
close(Fout);

chdir("./scripts");
@args=("sbatch","./script_MPI");
system(@args);
chdir("../");


