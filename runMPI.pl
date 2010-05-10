#!/usr/local/bin/perl
use Cwd;

#Use current directory as root
$cosmomc = cwd;

$params =  $ARGV[0];
$num = $ARGV[1];

$ini = $params;
if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

$path = $cosmomc;

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/csh -f
#PBS -N cosmomc
#PBS -l nodes=$num:ppn=2
#PBS -q workq
#PBS -r n
lamboot
cd $cosmomc
time mpirun N -O ./cosmomc $ini > ./scripts/$params.log
lamhalt
EMP
close(Fout);

chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
chdir("../");
