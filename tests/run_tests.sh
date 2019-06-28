cd /cosmomc/forutils
make ReleaseMPI

cd /cosmomc
#wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Code-v2.0.R2.00.tar.bz2
wget https://cdn.cosmologist.info/cosmobox/plc-3.0.tar.bz2
tar xvfj *.tar.bz2
cd plc-3.0
./waf configure --install_all_deps
./waf install
source ./bin/clik_profile.sh
cd ..
#wget http://irsa.ipac.caltech.edu/data/Planck/release_2/software/COM_Likelihood_Data-baseline_R2.00.tar.gz
wget https://cdn.cosmologist.info/cosmobox/test2019_kaml/baseline.tar.gz
#tar xvfz COM_Likelihood_Data-baseline_R2.00.tar.gz
tar xvfz baseline.tar.gz
ln -s $(pwd)/baseline/plc_3.0 ./data/clik
#rm -f COM_Likelihood_Data-baseline_R2.00.tar.gz
rm -f baseline.tar.gz

make

mpirun -np 1 --allow-run-as-root ./cosmomc test.ini

mpirun -np 1 --allow-run-as-root ./cosmomc test_planck.ini
rc=$?


exit $rc

