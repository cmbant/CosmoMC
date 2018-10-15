import copy
import re
from paramgrid import batchjob

# plik foregrounds have to be calculated a posteriori as before (before zre6p5 filter).

ini_dir = 'batch3/'
cov_dir = 'planck_covmats/'

defaults = ['common.ini']
getdist_options = {'ignore_rows': 0.3, 'marker[nrun]': 0,
                   'marker[w]': -1}

importanceDefaults = ['importance_sampling.ini']

# ranges for parameters when they are varied
params = dict()
params['w'] = '-0.99 -3. 1 0.02 0.02'
params['wa'] = '0 -3 2 0.05 0.05'
params['mnu'] = '0.02 0 5 0.1 0.03'
params['omegak'] = '-0.0008 -0.3 0.3 0.001 0.001'  # starting exactly on flat seems to confuse minimizer
params['nnu'] = '3.046 0.05 10 0.05 0.05'
params['nrun'] = '0 -1 1 0.005 0.001'
params['nrunrun'] = '0 -1 1 0.002 0.001'
params['r'] = '0 0 3 0.03 0.03'
params['Alens'] = '1 0 10 0.05 0.05'
params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
params['alpha1'] = '0 -1 1 0.0003 0.0003'
params['meffsterile'] = '0.1 0 3 0.1 0.03'
params['Aphiphi'] = '1 0 10 0.02 0.02'
params['Alensf'] = '1 0 10 0.03 0.03'

# extra parameters that are set only when specific parameters are varied.
param_extra_opts = {
    'mnu': {'num_massive_neutrinos': 3, 'neutrino_hierarchy': 'degenerate'},
    'meffsterile': {'param[mnu]': '0.06', 'param[nnu]': '3.1 3.046 10 0.05 0.05', 'num_massive_neutrinos': 1,
                    'accuracy_level': 1.2},
    'yhe': {'bbn_consistency': False},
    'r': {'compute_tensors': True},
    'nt': {'inflation_consistency': False, 'lmax_tensor': 1000}
}

# dataset names
lowl = 'lowl'
lensing = 'lensing'
lensonly = 'lensonly'
lowE = 'lowE'

WMAP = 'WMAP'
BAO = 'BAO'
HST = 'Riess18'
JLA = 'JLA'
Pantheon = 'Pantheon'
BAORSD = 'BAORSD'

lowEdata = 'simall_EE.ini'

BAOdata = 'BAO.ini'
RSDdata = 'BAODR12_RSD.ini'
HSTdata = 'HST_Riess2018.ini'

theta_prior = {'prior[theta]': '1.0409 0.0006'}

# dataset names
tauprior = {'prior[tau]': '0.055 0.009'}
tauname = 'tau055'

variant_tag = ['TT', 'TTTEEE']
variant_pol_tag = ['TE', 'EE']
variants = variant_tag

planck_highL_sets = []
planck_pol_sets = []
planck_vars = ['plikHM', 'CamSpecHM']

planck_ini = ['plik_rd12_HM_v22_%s.ini', 'nonclik_v10_7_%s.ini']
clean_ini = ['nonclik_v10_7_TT_clean.ini']
# planck_ini = ['plik_rd12_HM_v22_%s.ini', 'CAMspec_%s_clik14.ini']
planck_base = [[], []]

for planck, ini, base in zip(planck_vars, planck_ini, planck_base):
    for name, var in zip(variant_tag, variants):
        planck_highL_sets.append(batchjob.dataSet([planck, name], base + [ini % var]))
    for var in variant_pol_tag:
        planck_pol_sets.append(batchjob.dataSet([planck, var], base + [ini % var]))

baseTT = planck_highL_sets[0]
baseTTTEEE = planck_highL_sets[1]

WMAP9 = [[WMAP], ['WMAP.ini']]

likechecks = []

newCovmats = False


# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return [planck for planck in planck_vars if planck in jobItem.data_set.names] and (
                not 'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))


class importanceFilterSN:
    def wantImportance(self, jobItem):
        return 'JLA' not in jobItem.data_set.names and 'Pantheon' not in jobItem.data_set.names


class reion_importance(batchjob.importanceSetting):
    def wantImportance(self, jobItem):
        return [planck for planck in planck_vars if
                planck in jobItem.data_set.names] and not 'reion' in jobItem.data_set.names


class zre_importance(batchjob.importanceFilter):
    def wantImportance(self, jobItem):
        return [planck for planck in planck_vars if
                planck in jobItem.data_set.names] and not 'reion' in jobItem.data_set.names

    def filter(self, batch, jobItem):
        samples = jobItem.parent.getMCSamples(settings=getdist_options)
        pars = samples.getParams()
        samples.filter(pars.zrei >= 6.5)
        samples.ranges.setRange('zrei', [6.5, None])
        samples.saveChainsAsText(jobItem.chainRoot, properties={'burn_removed': True})


class importanceFilterNotOmegak:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set)


class importanceFilterBAO:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set) and jobItem.data_set.hasName('BAO')


class importanceFilterNnu:
    def wantImportance(self, jobItem):
        return 'nnu' in jobItem.param_set


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_lensingBAO = [[BAO, lensing], [BAOdata, 'lensing.ini'], importanceFilterNotOmegak()]
post_lensingPantheon = [[lensing, Pantheon], ['lensing.ini', 'Pantheon.ini'], importanceFilterSN()]
post_lensingJLA = [[lensing, JLA], ['lensing.ini', 'JLA_marge.ini'], importanceFilterNotOmegak()]

post_BAO = [[BAO], [BAOdata], importanceFilterNotOmegak()]
post_HST = [[HST], [HSTdata], importanceFilterNotOmegak()]
post_BAOJLA = [[BAO, JLA], [BAOdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_BAOPantheon = [[BAO, Pantheon], [BAOdata, 'Pantheon.ini'], importanceFilterNotOmegak()]
post_BAOHST = [[BAO, HST], [BAOdata, HSTdata], importanceFilterNotOmegak()]

post_BAOHSTJLA = [[BAO, JLA, HST], [BAOdata, 'JLA_marge.ini', HSTdata], importanceFilterNotOmegak()]
post_BAOHSTPantheon = [[BAO, Pantheon, HST], [BAOdata, 'Pantheon.ini', HSTdata], importanceFilterNotOmegak()]
post_BAOlensingPantheon = [[BAO, lensing, Pantheon], [BAOdata, 'lensing.ini', 'Pantheon.ini'],
                           importanceFilterNotOmegak()]

post_Pantheon = [[Pantheon], ['Pantheon.ini'], importanceFilterNotOmegak()]

post_CookeBBN = ['Cooke17']
post_Aver15 = [['Aver15'], ['Aver15BBN.ini'], importanceFilterNnu()]
post_BBN = [['Cooke17', 'Aver15'], ['Aver15BBN.ini', 'Cooke17BBN.ini'], importanceFilterNnu()]

# set up groups of parameters and data sets

groups = []

g = batchjob.jobGroup('main')
g.datasets = copy.deepcopy(planck_highL_sets)
for d in g.datasets:
    d.add(lowl)
    d.add(lowE, lowEdata)

g.params = [[], ['omegak'], ['mnu'], ['r'], ['nrun', 'r'], ['nnu'], ['nrun'], ['Alens'], ['yhe'], ['w'], ['alpha1']]
g.importanceRuns = [post_BAO, post_lensing, post_lensingBAO, post_HST, post_BBN]
groups.append(g)

gpol = batchjob.jobGroup('mainpol')
gpol.datasets = copy.deepcopy(planck_pol_sets)
for d in gpol.datasets:
    d.add(lowE, lowEdata)
gpol.params = [[], ['mnu'], ['nnu'], ['nrun'], ['Alens'], ['yhe'], ['r']]
gpol.importanceRuns = [post_BAO]
groups.append(gpol)

gpol = batchjob.jobGroup('polbao')
gpol.datasets = copy.deepcopy(planck_pol_sets)
for d in gpol.datasets:
    d.add(lowE, lowEdata)
    d.add(BAO, BAOdata)
gpol.params = [[], ['mnu'], ['nnu']]
gpol.importanceRuns = [post_lensing]
groups.append(gpol)

gpol = batchjob.jobGroup('pollensing')
gpol.datasets = copy.deepcopy(planck_pol_sets)
for d in gpol.datasets:
    d.add(lowE, lowEdata)
    d.add(lensing)
for d in list(gpol.datasets):
    d = d.copy().add(BAO, BAOdata).add('CookeDH', 'baryon_density.ini')
    gpol.datasets.append(d)
for d in copy.deepcopy(planck_pol_sets):
    d.add(lowE, lowEdata)
    d.add(lensing)
    d.add('CookeDH', 'baryon_density.ini')
    gpol.datasets.append(d)
gpol.params = [[]]
gpol.importanceRuns = []
groups.append(gpol)

gnotau = batchjob.jobGroup('nopoltau')
gnotau.params = [[]]
gnotau.datasets = copy.deepcopy(planck_highL_sets)
for d in gnotau.datasets:
    d.add(lowl)
for d in copy.deepcopy(planck_highL_sets):
    d.add(lowl)
    d.add(lensing)
    gnotau.datasets.append(d)
for d in copy.deepcopy(planck_highL_sets):
    d.add(lowl)
    d.add('reion', 'reion_tau.ini')
    gnotau.datasets.append(d)
gnotau.importanceRuns = [post_BAO]
groups.append(gnotau)

gnotau = batchjob.jobGroup('nopoltaumnu')
gnotau.params = [['mnu'], ['Alens']]
gnotau.datasets = []
for d in copy.deepcopy(planck_highL_sets):
    d.add(lowl)
    d.add(lensing)
    gnotau.datasets.append(d)
gnotau.importanceRuns = [post_BAO]
groups.append(gnotau)

g2 = batchjob.jobGroup('ext')
g2.datasets = copy.deepcopy(g.datasets)
g2.params = [['nnu', 'meffsterile'], ['nnu', 'mnu'], ['nnu', 'yhe']]
g2.importanceRuns = [post_BAO, post_HST, post_lensing, post_lensingBAO]
groups.append(g2)

g3 = batchjob.jobGroup('geom')
g3.params = [['omegak'], ['w'], ['w', 'wa']]
g3.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    g3.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    d.add(HST, HSTdata)
    d.add(JLA)
    g3.datasets.append(d)

g3.importanceRuns = [post_lensing, post_lensingPantheon]
groups.append(g3)

g3 = batchjob.jobGroup('de')
g3.params = [['w'], ['w', 'wa']]
g3.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    d.add(Pantheon)
    g3.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    d.add(HST, HSTdata)
    d.add(Pantheon)
    g3.datasets.append(d)
g3.importanceRuns = [post_lensing]
groups.append(g3)

g6 = batchjob.jobGroup('lensing')
g6.datasets = copy.deepcopy(g.datasets)
for d in g6.datasets:
    d.add(lensing)

g6.params = [['omegak'], ['mnu'], ['Alens']]
g6.importanceRuns = [post_BAO]

groups.append(g6)

inflat = batchjob.jobGroup('inflat')
inflat.datasets = copy.deepcopy(g.datasets)
for d in inflat.datasets:
    d.add(lensing)
inflat.params = [['r'], ['nrun', 'r']]
inflat.importanceRuns = [post_BAO]
groups.append(inflat)

gbest = batchjob.jobGroup('basebest')
gbest.datasets = copy.deepcopy(g.datasets)
for d in gbest.datasets:
    d.add(lensing)

gbest.params = [[]]
gbest.importanceRuns = [post_BAO, post_HST, post_BAOHST, post_Pantheon, post_BAOHSTJLA, post_BAOPantheon,
                        post_BAOHSTPantheon]
groups.append(gbest)

g7 = batchjob.jobGroup('mnu')
g7.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    g7.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    d.add(BAO, BAOdata)
    g7.datasets.append(d)

g7.params = [['mnu'], ['nnu', 'meffsterile'], ['nnu', 'mnu']]
g7.importanceRuns = [post_Pantheon, post_Aver15, post_BBN]
groups.append(g7)

gnonbbn = batchjob.jobGroup('nonbbn')
gnonbbn.datasets = copy.deepcopy(g.datasets)
for d in gnonbbn.datasets:
    d.add('Aver15', 'Aver15BBN')
gnonbbn.params = [['yhe'], ['nnu', 'yhe']]
gnonbbn.importanceRuns = [post_BAO, post_lensing, post_lensingBAO]
groups.append(gnonbbn)

gnnu = batchjob.jobGroup('nnu')
gnnu.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    gnnu.datasets.append(d)
gnnu.params = [['nnu']]
gnnu.importanceRuns = [post_lensingJLA, post_lensingPantheon, post_lensing, post_Aver15, post_BBN]
groups.append(gnnu)

gHST = batchjob.jobGroup('HST')
gHST.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(HST, HSTdata)
    gHST.datasets.append(d)
gHST.params = [['nnu']]
gHST.importanceRuns = [post_BAO, post_BAOPantheon, post_lensing, post_lensingBAO, post_BAOlensingPantheon]
groups.append(gHST)

gclean = batchjob.jobGroup('cleaned')
d = batchjob.dataSet(['CleanedCamSpecHM', 'TT'], [clean_ini])
d.add(lowl)
d.add(lowE, lowEdata)
gclean.datasets = [d]
gclean.params = [[], ['mnu'], ['nnu'], ['yhe'], ['Alens'], ['omegak'], ['nrun'], ['r'], ['w']]
groups.append(gclean)

# Things mainly for the lensing paper
glens = batchjob.jobGroup('lensonly')
# get this well converged so importance sampling might work
lensdata = batchjob.dataSet(['lensing', 'lenspriors'], [lensonly, 'lensonly_priors', {'MPI_Converge_Stop': 0.002}])
glens.datasets = [lensdata]
glens.datasets += [lensdata.copy().add('theta', theta_prior)]
lensdata = lensdata.copy().add(BAO, BAOdata)
glens.datasets.append(lensdata)
lensdata = lensdata.copy().add('theta', theta_prior)
glens.datasets.append(lensdata)
glens.params = [[], ['mnu']]
glens.importanceRuns = [post_Pantheon]

varnames = ['agr2', 'conslmin40', 'agrlmax425', 'ptt', 'pttagr2', 'bfcl', 'agr2bfcl', 'linear', 'acc', 'agr2acc',
            'takahashi', 'agr2takahashi', 'Apr6']
if True:
    # Consistency checks mainly for the lensing paper
    base = 'smicadx12_Dec5_ftl_mv2_n'
    vars = [('dclpp_p_teb_agr2_CMBmarged', {}),
            ('dclpp_p_teb_consext8_CMBmarged', {'cmb_dataset[lensing,use_min]': 2}),
            ('dclpp_p_teb_agr2_CMBmarged', {'cmb_dataset[lensing,use_max]': 9}),
            ('dclpttptt_p_teb_consext8_CMBmarged', {}),
            ('dclpttptt_p_teb_agr2_CMBmarged', {}),
            ('dclpp_p_teb_consext8_lensonly', {}),
            ('dclpp_p_teb_agr2_lensonly', {}),
            #           ('dclpp_p_teb_consext8_lensonly', {
            #               'cmb_dataset[lensing,linear_correction_bin_window_fix_cl_file]': '../base_omegak_plikHM_TTTEEE_lowl_lowE.minimum.theory_cl'}),
            ('dclpp_p_teb_consext8_CMBmarged', {'redo_theory': True, 'redo_cls': True, 'use_nonlinear_lensing': False}),
            ('dclpp_p_teb_consext8_CMBmarged',
             {'redo_theory': True, 'redo_cls': True, 'accuracy_level': 1.5, 'k_eta_max_scalar': 50000}),
            ('dclpp_p_teb_agr2_CMBmarged',
             {'redo_theory': True, 'redo_cls': True, 'accuracy_level': 1.5, 'k_eta_max_scalar': 50000}),
            ('dclpp_p_teb_consext8_CMBmarged',
             {'redo_theory': True, 'redo_cls': True, 'halofit_version': 4}),
            ('dclpp_p_teb_agr2_CMBmarged',
             {'redo_theory': True, 'redo_cls': True, 'halofit_version': 4}),
            ('smicadx12_Apr6_ndclpp_p_teb_consext8_CMBmarged', {}),
            ]

    for name, var in zip(varnames, vars):
        tag, opt = var
        dic = {'redo_likelihood': True, 'redo_add': False,
               'cmb_dataset[lensing]': '%%DATASETDIR%%planck_lensing_2017/%s%s.dataset' % (
                   base if 'smica' not in tag else '', tag)}
        dic.update(opt)
        if 'ptt' in name and not 'CMBmarged' in name:
            dic['cmb_dataset[lensing,linear_correction_bin_window_fix_cl]'] = 'TT'
        glens.importanceRuns.append([[name], [dic]])
glens.extra_opts = {'sampling_method': 1}  # no fast params
groups.append(glens)

# Things mainly for the lensing paper
glensTT = batchjob.jobGroup('lensTT')
base = 'smicadx12_Dec5_ftl_mv2_n'
tag = 'dclpttptt_p_teb_agr2_CMBmarged'
lensdata = batchjob.dataSet(['lensing', 'lenspriors', 'pttagr2'],
                            [lensonly, 'lensonly_priors',
                             {'cmb_dataset[lensing]': '%%DATASETDIR%%planck_lensing_2017/%s%s.dataset' % (base, tag)}])
glensTT.datasets = [lensdata]
glensTT.datasets += [lensdata.copy().add('theta', theta_prior)]
lensdata = lensdata.copy().add(BAO, BAOdata)
glensTT.datasets.append(lensdata)
lensdata = lensdata.copy().add('theta', theta_prior)
glensTT.datasets.append(lensdata)
glensTT.params = [[], ['mnu']]
groups.append(glensTT)

glens = batchjob.jobGroup('lensonlyastro')
lensdata = batchjob.dataSet(['lensing', 'DESpriors'], [lensonly, 'DES_astro_priors'])
glens.datasets = [lensdata]
glens.datasets += [lensdata.copy().add(BAO, BAOdata)]
lensdata = lensdata.copy().add('CookeDH', 'baryon_density.ini')
glens.datasets += [lensdata]
lensdata = lensdata.copy().add(BAO, BAOdata)
glens.datasets += [lensdata]
glens.param_extra_opts = {'mnu': {'param[mnu]': '0.07 0.05 1 0.1 0.03'}}
glens.params = [[], ['mnu']]
glens.extra_opts = {'sampling_method': 1}  # no fast params
groups.append(glens)

gphi = batchjob.jobGroup('Aphiphi')
gphi.params = [['Aphiphi']]
gphi.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    gphi.datasets.append(d)
gphi.importanceRuns = []
groups.append(gphi)

gphi = batchjob.jobGroup('Alens')
gphi.params = [[], ['Alens']]
gphi.datasets = []
for d in copy.deepcopy(planck_highL_sets):
    gphi.datasets.append(d)
    dlow = d.copy()
    dlow.add(lowl)
    gphi.datasets.append(dlow)
    dtau = d.copy()
    dtau.add(lowE, lowEdata)
    gphi.datasets.append(dtau)
gphi.importanceRuns = [post_BAO]
groups.append(gphi)

gWMAP = batchjob.jobGroup('WMAP')
gWMAP.params = [[]]
gWMAP.datasets = [WMAP9]

gWMAP.importanceRuns = [post_BAO]
groups.append(gWMAP)

for bk in ['BK14']:
    gbkp = batchjob.jobGroup(bk)
    gbkp.datasets = []
    for d in copy.deepcopy(planck_highL_sets):
        d.add(lowl)
        d.add(lowE, lowEdata)
        d.add(bk)
        gbkp.datasets.append(d)
    for d in [copy.deepcopy(baseTTTEEE)]:
        d.add(lowl)
        d.add(lowE, lowEdata)
        d.add(bk)
        d.add(lensing)
        gbkp.datasets.append(d)

    gbkp.params = [['r'], ['nrun', 'r']]
    gbkp.importanceRuns = [post_BAO, post_lensing, post_lensingBAO]
    groups.append(gbkp)

DESdatapriors = [batchjob.dataSet(['DES', 'lenspriors'], ['DES', 'lensonly_priors']),
                 batchjob.dataSet(['DESlens', 'lenspriors'], ['DES_lensing', 'lensonly_priors'])]
gWL = batchjob.jobGroup('DES')
gWL.datasets = copy.deepcopy(DESdatapriors)
for d in copy.deepcopy(DESdatapriors):
    d.add('lensing', lensonly)
    gWL.datasets.append(d)
for d in copy.deepcopy(DESdatapriors):
    d.add(BAO, BAOdata)
    gWL.datasets.append(d)
for d in copy.deepcopy(DESdatapriors):
    d.add('lensing', lensonly)
    d.add(BAO, BAOdata)
    gWL.datasets.append(d)
gWL.params = [[], ['mnu']]
gWL.importanceRuns = []
groups.append(gWL)

DESdatapriors = [batchjob.dataSet(['DES', 'DESpriors'], ['DES', 'DES_astro_priors']),
                 batchjob.dataSet(['DESlens', 'DESpriors'], ['DES_lensing', 'DES_astro_priors']),
                 batchjob.dataSet(['DESwt', 'DESpriors'], ['DES_wt', 'DES_astro_priors'])
                 ]

gWL = batchjob.jobGroup('DESastro')
gWL.datasets = copy.deepcopy(DESdatapriors)
for d in copy.deepcopy(DESdatapriors):
    d.add('lensing', 'lensonly')
    gWL.datasets.append(d)
for d in copy.deepcopy(DESdatapriors):
    d.add(BAO, BAOdata)
    d.add('CookeDH', 'baryon_density.ini')
    gWL.datasets.append(d)
for d in copy.deepcopy(DESdatapriors):
    d.add('lensing', lensonly)
    d.add(BAO, BAOdata)
    d.add('CookeDH', 'baryon_density.ini')
    gWL.datasets.append(d)
gWL.param_extra_opts = {'mnu': {'param[mnu]': '0.07 0.05 1 0.1 0.03'}}
gWL.params = [[], ['mnu']]
gWL.importanceRuns = []
groups.append(gWL)

gDESPlanck = batchjob.jobGroup('DESPlanck')
gDESPlanck.datasets = []
for d in [copy.deepcopy(baseTTTEEE)]:
    d.add(lowl)
    d.add(lowE, lowEdata)
    d.add('DES')
    gDESPlanck.datasets.append(d)
for d in [copy.deepcopy(baseTTTEEE)]:
    d.add(lowl)
    d.add(lowE, lowEdata)
    d.add('DESlens', 'DES_lensing')
gDESPlanck.datasets.append(d)
# for d in copy.deepcopy(gDESPlanck.datasets):
#    d.add(lensing)
#    gDESPlanck.datasets.append(d)
gDESPlanck.params = [[], ['mnu']]
gDESPlanck.importanceRuns = [post_BAO, post_lensing, post_lensingBAO]
groups.append(gDESPlanck)

gext = batchjob.jobGroup('twos')
d = copy.deepcopy(baseTTTEEE)
d.add(lowl)
d.add(lowE, lowEdata)
gext.datasets = [d]
gext.params = [['nrun', 'nrunrun'], ['nnu', 'nrun']]
gext.importanceRuns = [post_BAO, post_lensing, post_lensingBAO]
groups.append(gext)

gext = batchjob.jobGroup('big')
d = copy.deepcopy(baseTTTEEE)
d.add(lowl)
d.add(lowE, lowEdata)
d.add(BAO, BAOdata)
d.add(HST, HSTdata)
d.add(Pantheon)
d.add(lensing)
gext.datasets = [d]
gext.params = [['nrun', 'nnu', 'w', 'mnu']]
gext.importanceRuns = []
groups.append(gext)

gBBN = batchjob.jobGroup('lensingBBN')
BBN = batchjob.dataSet(['lensing', 'lenspriors', BAO, 'Cooke17', 'Aver15'],
                       [lensonly, 'Aver15BBN', 'Cooke17BBN', BAOdata, 'lensonly_priors',
                        {'prior[omegabh2]': ''}])
gBBN.datasets = [BBN]
gBBN.datasets += [BBN.copy().add('theta', theta_prior)]
gBBN.params = [['nnu'], ['nnu', 'mnu']]
gBBN.extra_opts = {'sampling_method': 1}
gBBN.importanceRuns = [post_Pantheon]
groups.append(gBBN)

noCMB = {'lmin_store_all_cmb': 0, 'lmin_computed_cl': 0, 'param[ns]': 0.96, 'param[logA]': 3, 'param[tau]': 0.055,
         'get_sigma8': False}

gBBN = batchjob.jobGroup('BBN')
BBN = batchjob.dataSet([BAO, 'Cooke17'], [BAOdata, 'Cooke17BBN', noCMB])
gBBN.datasets = [BBN]
gBBN.datasets += [BBN.copy().add('Pantheon')]
gBBN.datasets += [BBN.copy().add('JLA')]
gBBN.datasets += [BBN.copy().add('Pantheon').add('theta', theta_prior)]
gBBN.datasets += [BBN.copy().add('theta', theta_prior)]
gBBN.params = [[], ['mnu']]
gBBN.extra_opts = {'sampling_method': 1}
groups.append(gBBN)

gBBN = batchjob.jobGroup('BBNnnu')
BBN1 = batchjob.dataSet([BAO, 'Cooke17', 'Aver15'], ['Aver15BBN', 'Cooke17BBN', BAOdata, noCMB])
BBN2 = batchjob.dataSet([BAO, 'Cooke17Marc', 'Aver15'],
                        ['Aver15BBN', {'abundance_dataset[Cooke17Marc]': '%DATASETDIR%D_Cooke2017_marcucci.dataset'},
                         BAOdata, noCMB])
BBN3 = batchjob.dataSet([BAO, 'Cooke17Adel', 'Aver15'],
                        ['Aver15BBN', {'abundance_dataset[Cooke17Adel]': '%DATASETDIR%D_Cooke2017_adelberger.dataset'},
                         BAOdata, noCMB])
gBBN.datasets = []
for BBN in [BBN1, BBN2, BBN3]:
    gBBN.datasets += [BBN]
    gBBN.datasets += [BBN.copy().add('Pantheon')]
    gBBN.datasets += [BBN.copy().add('theta', theta_prior)]
    gBBN.datasets += [BBN.copy().add('Pantheon').add('theta', theta_prior)]
gBBN.params = [['nnu'], ['nnu', 'mnu']]
gBBN.extra_opts = {'sampling_method': 1}
groups.append(gBBN)

# add zre prior for all runs
importance_filters = [zre_importance(['zre6p5'])]

skip = []

covWithoutNameOrder = ['lenspriors', 'CookeDH', 'pttagr2', HST, 'JLA', Pantheon, BAORSD, 'lensing', 'DESpriors', 'DES',
                       'DESlens', 'BAO', 'reion', 'abundances', 'theta', 'Aver15'] + varnames
covNameMappings = {HST: 'HST', 'CleanedCamSpecHM': 'CamSpecHM', 'Cooke17Adel': 'Cooke17', 'Cooke17Marc': 'Cooke17'}

# try to match run to exisitng covmat
covrenames = []
for planck in planck_vars:
    covrenames.append([planck, 'planck'])

covrenames.append(['CamSpecHM', 'plikHM'])
covrenames.append(['lensing_lenspriors', 'lensonly'])
covrenames.append(['lensing', 'lensonly'])
covrenames.append(['Alensf', 'Alens'])
covrenames.append(['_Aphiphi', ''])
covrenames.append(['Pantheon', 'JLA'])
covrenames.append(['_Aver15', ''])

covrenames.append(['_r', ''])
covrenames.append(['_w', ''])
covrenames.append(['_alpha1', ''])

covrenames.append(['DES_BAO', 'BAO_lensonly'])
covrenames.append(['DESlens_BAO', 'BAO_lensonly'])
covrenames.append(['DES_lensonly', 'lensonly'])
covrenames.append(['DES', 'lensonly'])
covrenames.append(['DESwt', 'DES'])
covrenames.append(['mnu_DES', 'DES'])
covrenames.append(['Riess18', 'HST'])
covrenames.append(['DESwt_DESpriors_lensing', 'DES_DESpriors_lensonly'])
covrenames.append(['DESwt_DESpriors_BAO_CookeDH', 'DES_DESpriors_BAO'])
covrenames.append(['DESwt_DESpriors_lensing_BAO_CookeDH', 'DES_DESpriors_BAO'])
covrenames.append(['DESpriors_BAO', 'DESpriors_CookeDH_BAO'])
covrenames.append(['base_mnu_plikHM_TTTEEE_lowl_lowE_DES', 'base_mnu_plikHM_TTTEEE_lowl_lowE'])
covrenames.append(['base_mnu_plikHM_TTTEEE_lowl_lowE_DESlens', 'base_mnu_plikHM_TTTEEE_lowl_lowE'])
covrenames.append(['lensing_lenspriors_theta', 'lensing_lenspriors_theta_BAO'])
covrenames.append(['lowl_lensing', 'lowl_lowE'])
covrenames.append(['lowl', 'lowl_lowE'])
covrenames.append(['TT', 'TT_lowl_lowE'])
covrenames.append(['TTTEEE', 'TTTEEE_lowl_lowE'])
covrenames.append(['TTTEEE_lowE', 'TTTEEE_lowl_lowE'])
covrenames.append(['TT_lowE', 'TT_lowl_lowE'])

covrenames.append(
    ['nrun_nnu_w_mnu_plikHM_TTTEEE_lowl_lowE_BAO_Riess18_Pantheon_lensing', 'w_plikHM_TTTEEE_lowl_lowE_BAO_Pantheon'])
