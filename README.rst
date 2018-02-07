===================
CosmoMC
===================
:CosmoMC:  Fortran 2008 parallelized MCMC sampler (general and cosmology)
:Homepage: http://cosmologist.info/cosmomc/

Description and installation
=============================

For full details see the `ReadMe <http://cosmologist.info/cosmomc/readme.html>`_.

Algorithm details
==================

See the latest `paper <http://arxiv.org/abs/1304.4473>`_.

GetDist
===================

CosmoMC includes the GetDist python sample analysis and plotting package, which is
also `available separately <http://getdist.readthedocs.org/en/latest/>`_.

Branches
=============================

The master branch contains latest changes to the main release version.

.. image:: https://secure.travis-ci.org/cmbant/CosmoMC.png?branch=master
  :target: https://secure.travis-ci.org/cmbant/CosmoMC/builds

The devel branch is a development version, using CAMB devel branch which integrates 
CAMB and CAMB sources (though CAMB sources functions are not available via CosmoMC yet).
Includes run-time changing of dark energy model between fluid and PPF modes (easily extended).
Shared general function now taken from the `forutils <https://github.com/cmbant/forutils>`_ library.

.. image:: https://secure.travis-ci.org/cmbant/CosmoMC.png?branch=devel
  :target: https://secure.travis-ci.org/cmbant/CosmoMC/builds

Both branches now have a travis unit test to check they work with the Planck 2015 data. The test
does a test install of forutils, CosmoMC and the Planck likelihood code and checks the likelihood is as expected.
See tests/run_tests.sh for the setup and test code. There are small changes in the absolute likelihood value between branches
due to small changes in the CAMB version (e.g. implied optical depth changes very slightly due to changes in time sampling).

=============

.. raw:: html

    <a href="http://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex.png" height="170px"></a>
    <a href="http://erc.europa.eu/"><img src="https://erc.europa.eu/sites/default/files/content/erc_banner-vertical.jpg" height="200px"></a>
