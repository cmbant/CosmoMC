
cd /cosmomc

if ["$GCCTAG" != "devel"]
then

wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Code-v2.0.R2.00.tar.bz2
tar xvfj COM_Likelihood_Code-v2.*.tar.bz2
cd plc-2.0
./waf configure --install_all_deps
./waf install
source ./bin/clik_profile.sh
cd ..
wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Data-baseline_R2.00.tar.gz
tar xvfz COM_Likelihood_Data-baseline_R2.00.tar.gz
ln -s $(pwd)/plc_2.0 ./data/clik
rm -f COM_Likelihood_Data-baseline_R2.00.tar.gz

fi

make

mpirun -np 1 --allow-run-as-root ./cosmomc test.ini

if ["$GCCTAG" != "devel"]
then

mpirun -np 1 --allow-run-as-root ./cosmomc test_planck.ini

fi

rc=$?


exit $rc

