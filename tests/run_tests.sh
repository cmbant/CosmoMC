
cd /cosmomc
#wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Code-v2.0.R2.00.tar.bz2
wget https://cdn.cosmologist.info/cosmobox/plc-2.1_py3.tar.bz2
tar xvfj *.tar.bz2
mv plc-2.1_py3 plc-2.0
cd plc-2.0
./waf configure --install_all_deps
./waf install
source ./bin/clik_profile.sh
cd ..
wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Data-baseline_R2.00.tar.gz
tar xvfz COM_Likelihood_Data-baseline_R2.00.tar.gz
ln -s $(pwd)/plc_2.0 ./data/clik
rm -f COM_Likelihood_Data-baseline_R2.00.tar.gz


make

mpirun -np 1 --allow-run-as-root ./cosmomc test.ini
rc=$?

mpirun -np 1 --allow-run-as-root ./cosmomc test_planck.ini
rc=$?


exit $rc

