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

Related code
==================

The new Python `Cobaya <https://github.com/CobayaSampler/cobaya>`_ sampling package incorporates a 
version of CosmoMC's sampler and most other CosmoMC features, but has more general speed optimization and
general support of multiple inter-dependent theory and likelihood codes.


Branches
=============================

The master branch contains latest changes to the main release version, using latest CAMB 1.x.

.. image:: https://travis-ci.org/cmbant/CosmoMC.svg?branch=master
  :target: https://travis-ci.org/cmbant/CosmoMC/builds

The planck2018 branch contains the configuration used for the final Planck 2018 analysis, with 
corresponding CAMB version.

The devel branch is a development branch.

=============

.. raw:: html

    <a href="http://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex.png" height="170px"></a>
    <a href="http://erc.europa.eu/"><img src="https://erc.europa.eu/sites/default/files/content/erc_banner-vertical.jpg" height="200px"></a>
