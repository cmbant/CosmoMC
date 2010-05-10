#!/usr/local/bin/perl

#Example from Cambridge Leda machine
#Sets 2 processes per node 
#(so quad-core each process can use 2 openmp threads if set in .ini or environment)
#perl runMPI_leda.pl params 3
#should run 6 chains, 2 per the three nodes.

use Cwd;

#Use current directory as root
$cosmomc = cwd;
$cosmomc =~ s/export\///;

$params =  $ARGV[0];
$num = $ARGV[1];

$ini = $params;
if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

$path = $cosmomc;

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/csh -f
#PBS -N cosmomc
#PBS -l nodes=$num:ppn=2,walltime=10:00:00
#PBS -r n
#PBS -W x=NACCESSPOLICY:SINGLEJOB
cd $cosmomc

echo Running on host \`hostname\`
echo Time is \`date\`
echo Directory is \`pwd\`
echo PBS job ID is \$PBS_JOBID
echo This jobs runs on the following machines:
echo \`cat \$PBS_NODEFILE | uniq\`

#! Create a machine file for MPI
cat \$PBS_NODEFILE | uniq > scripts/machine.file.\$PBS_JOBID

time /usr/local/openmpi/intel10/64/1.2.2/bin/mpirun ./cosmomc $ini > ./scripts/$params.log
EMP
close(Fout);

chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
chdir("../");
