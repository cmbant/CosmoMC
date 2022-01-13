git gfortran -v

cd /cosmomc


make

mpirun -np 1 --allow-run-as-root ./cosmomc test.ini
rc=$?



exit $rc

