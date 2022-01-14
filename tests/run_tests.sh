gfortran -v

# Getting Choleksy errors: bug in openBLAS? https://github.com/xianyi/OpenBLAS/pull/2286
export OPENBLAS_CORETYPE=Haswell
cd /cosmomc
wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Code-v3.0_R3.01.tar.gz
tar xvf *.tar.gz
cd code/plc_3.0/plc-3.01
./waf configure --install_all_deps
./waf install
source ./bin/clik_profile.sh
cd ../../..
wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz
tar xvfz *baseline*.tar.gz
ln -s $(pwd)/baseline/plc_3.0 ./data/clik_14.0
rm -f *baseline*.tar.gz


make

mpirun -np 1 --allow-run-as-root ./cosmomc test.ini
rc=$?

mpirun -np 1 --allow-run-as-root ./cosmomc test_planck.ini
rc=$?


exit $rc
