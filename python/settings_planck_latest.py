import batchJob, copy, re

ini_dir = 'batch2/'
cov_dir = 'planck_covmats/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
lowl = 'lowl'
lensing = 'lensing'
lensonly = 'lensonly'
WLonly = 'WLonly'
WLonlyHeymans = 'WLonlyHeymans'  # just used as a check against less conservative cuts

highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'H070p6'
JLA = 'JLA'

BAOdata = 'BAO.ini'
HSTdata = 'HST_GPE70p6.ini'

RSDdata = 'BAO_RSD.ini'
BAORSD = 'BAORSD';
WL = 'WL'

Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
lowTEB = 'lowTEB'
# dataset names
tauprior = {'prior[tau]':'0.07 0.02'}
tauname = 'tau07'
WMAPtau = {'prior[tau]':'0.09 0.013'}


camspec_detsets = ['nonclik_detsets.ini']
camspec_CS = ['nonclik.ini']


variant_tag = ['TTTEEE', 'TT']
variant_pol_tag = ['TE', 'EE']
variants = variant_tag

planck_highL_sets = []
planck_pol_sets = []
planck_vars = ['CamSpecHM']
planck_ini = ['CAMspec_%s.ini']
planck_base = [camspec_CS]

if True:
    planck_vars += ['plikHM']
    planck_ini += ['plik_dx11dr2_HM_v15_%s.ini']
    planck_base += [[]]
for planck, ini, base in zip(planck_vars, planck_ini, planck_base):
    for name, var in zip(variant_tag, variants):
        planck_highL_sets.append(batchJob.dataSet([planck , name], base + [ ini % (var)]))
    for var in variant_pol_tag:
        planck_pol_sets.append(batchJob.dataSet([planck , var], base + [ ini % (var)]))


WMAP9 = [[WMAP], ['WMAP.ini']]

likechecks = []
likechecks.append(batchJob.dataSet(['CamSpecDS', 'TT'], camspec_detsets + ['CAMspec_TT.ini']))
likechecks.append(batchJob.dataSet(['Mspec', 'TT'], ['mspec_dx11d_HM_v1_TT.ini']))


start_at_bestfit = False
newCovmats = False

# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return [p for planck in planck_vars if planck in jobItem.data_set.names] and (not'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))

class importanceFilterZre:
    def wantImportance(self, jobItem):
        return [p for planck in planck_vars if planck in jobItem.data_set.names] and not 'reion' in jobItem.data_set.names

class importanceFilterNotOmegak:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set)

class importanceFilterAbundance:
    def wantImportance(self, jobItem):
        return 'nnu' in jobItem.param_set and not 'yhe' in jobItem.param_set


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_BAO = [[BAO], [BAOdata], importanceFilterNotOmegak()]
post_HST = [[HST], [HSTdata], importanceFilterNotOmegak()]
post_JLA = [[JLA], ['JLA_marge.ini'], importanceFilterNotOmegak()]
post_nonBAO = [[HST, JLA], [HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_nonCMB = [[BAO, HST, JLA], [BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_all = [[lensing, BAO, HST, JLA], [lensing, BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_allnonBAO = [[lensing, HST, JLA], [lensing, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]

post_WP = [[ 'WMAPtau'], [WMAPtau, {'redo_no_new_data':'T'}]]
post_abundance = [['abundances'], ['abundances.ini', {'redo_likes':'F', 'redo_add':'T'}], importanceFilterAbundance()]
post_zre = [['zre6p5'], ['zre_prior.ini'], importanceFilterZre()]

# set up groups of parameters and data sets

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(planck_highL_sets)
for d in g.datasets:
    d.add(lowTEB)

g.params = [[], ['omegak'], ['mnu'], ['r'], ['nrun', 'r'], ['nnu'], ['nrun'], ['Alens'], ['yhe'], ['w'], ['alpha1']]
g.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST, post_all]
groups.append(g)


gpol = batchJob.jobGroup('mainpol')
gpol.datasets = copy.deepcopy(planck_pol_sets)
for d in gpol.datasets:
    d.add(lowTEB)

gpol.params = [[], ['mnu'], ['nnu'], ['nrun'], ['Alens'], ['yhe'], ['alpha1']]
gpol.importanceRuns = []
groups.append(gpol)


if False:
    ghigh = batchJob.jobGroup('highL')
    ghigh.datasets = copy.deepcopy(g.datasets)
    for d in ghigh.datasets:
        d.add(highL)

    ghigh.params = [[], ['omegak'], ['mnu'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]
    ghigh.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST, post_all]
    groups.append(ghigh)


g2 = batchJob.jobGroup('ext')
g2.datasets = copy.deepcopy(g.datasets)
g2.params = [ ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['nnu', 'yhe']]
g2.importanceRuns = [post_BAO, post_JLA, post_HST, post_nonCMB]
groups.append(g2)

g3 = batchJob.jobGroup('geom')
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

g3.importanceRuns = [post_lensing]
groups.append(g3)



g5 = batchJob.jobGroup('nopoltau')
g5.params = [[]]
g5.datasets = copy.deepcopy(planck_highL_sets)
for d in g5.datasets:
    d.add(lowl)
for d in copy.deepcopy(planck_highL_sets):
    d.add(lowl)
    d.add(lensing)
    g5.datasets.append(d)
for d in copy.deepcopy(planck_highL_sets):
    d.add(lowl)
    d.add('reion', 'reion_tau.ini')
    g5.datasets.append(d)
g5.importanceRuns = [post_BAO, post_nonCMB, post_zre]
groups.append(g5)


gpolnopoltau = batchJob.jobGroup('polnopoltau')
gpolnopoltau.params = [[]]
gpolnopoltau.datasets = copy.deepcopy(planck_pol_sets)
for d in copy.deepcopy(planck_pol_sets):
    d.add(lensing)
    gpolnopoltau.datasets.append(d)

gpolnopoltau.importanceRuns = [post_BAO, post_nonCMB]
groups.append(gpolnopoltau)

glowllens = batchJob.jobGroup('lowllensing')
glowllens.params = [['mnu']]
glowllens.datasets = copy.deepcopy(planck_highL_sets)
for d in glowllens.datasets:
    d.add(lowl)
    d.add(lensing)

glowllens.importanceRuns = [post_BAO, post_nonCMB]
groups.append(glowllens)

g6 = batchJob.jobGroup('lensing')
g6.datasets = copy.deepcopy(g.datasets)
for d in g6.datasets:
    d.add(lensing)
    d.add(None, {'redo_theory':'F'})

g6.params = [['omegak'], ['mnu'], ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['Alens'], ['nnu', 'meffsterile', 'r']]
g6.importanceRuns = []
groups.append(g6)

gbest = batchJob.jobGroup('basebest')
gbest.datasets = copy.deepcopy(g.datasets)
for d in gbest.datasets:
    d.add(lensing)
    d.add(None, {'redo_theory':'F'})

gbest.params = [[]]
gbest.importanceRuns = [post_BAO, post_JLA, post_HST, post_nonCMB]
groups.append(gbest)


g7 = batchJob.jobGroup('mnu')
g7.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    g7.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    d.add(BAO, BAOdata)
    d.add(None, {'redo_theory':'F'})
    g7.datasets.append(d)

g7.params = [['mnu'], ['nnu', 'meffsterile']]
g7.importanceRuns = [post_JLA, post_HST, post_nonBAO]
groups.append(g7)

gnnu = batchJob.jobGroup('nnu')
gnnu.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    gnnu.datasets.append(d)

gnnu.params = [['nnu']]
gnnu.importanceRuns = [post_nonBAO, post_allnonBAO]
groups.append(gnnu)

if False:
    gmulti = batchJob.jobGroup('multi')
    gmulti.params = [['nnu', 'w'], ['mnu', 'w']]
    gmulti.datasets = []
    for d in copy.deepcopy(g.datasets):
    #    d.covmat = 'planck_covmats/base_w_BAO_HST_JLA_TTTEEE_lensing_lowTEB_plik.covmat'
        d.add(lensing)
        d.add(BAO, BAOdata)
        d.add(JLA)
        gmulti.datasets.append(d)
    gmulti.importanceRuns = [post_HST, post_abundance]
    groups.append(gmulti)


# Things mainly for the lensing paper

glens = batchJob.jobGroup('lensonly')
lensdata = [batchJob.dataSet(lensonly)]
glens.datasets = copy.deepcopy(lensdata)
for d in copy.deepcopy(lensdata):
    d.add(BAO, BAOdata)
    glens.datasets.append(d)
for d in copy.deepcopy(lensdata):
    d.add(HST, HSTdata)
    glens.datasets.append(d)
# for d in copy.deepcopy(lensdata):
#    d.add(HST, HSTdata)
#    d.add('widerns', {'prior[ns]': '0.98 0.05'})
#    glens.datasets.append(d)
for d in copy.deepcopy(lensdata):
    d.add('theta', {'param[theta]':'1.0408'})
    glens.datasets.append(d)
for d in copy.deepcopy(lensdata):
    d.add('theta', {'param[theta]':'1.0408'})
    d.add(BAO, BAOdata)
    glens.datasets.append(d)
glens.params = [[], ['mnu']]
glens.importanceRuns = []
groups.append(glens)

glens = batchJob.jobGroup('lensonlyext')
glens.datasets = copy.deepcopy(lensdata)
glens.params = [['nnu'], ['nnu', 'meffsterile']]
glens.importanceRuns = []
groups.append(glens)

gphi = batchJob.jobGroup('Aphiphi')
gphi.params = [['Aphiphi']]
gphi.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    gphi.datasets.append(d)
gphi.importanceRuns = []
groups.append(gphi)

gphi = batchJob.jobGroup('altAlens')
gphi.params = [['Alensf']]
gphi.datasets = copy.deepcopy(g.datasets)
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    gphi.datasets.append(d)
gphi.importanceRuns = []

groups.append(gphi)


extdata = batchJob.jobGroup('extdata')
extdata.params = [[], ['mnu'], ['nnu', 'meffsterile']]
extdata.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(WL)
    d.add(lensing)
    extdata.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(BAORSD, RSDdata)
    extdata.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(BAORSD, RSDdata)
    d.add(lensing)
    d.add(WL)
    d.add(JLA)
    d.add(HST, HSTdata)
    extdata.datasets.append(d)

extdata.importanceRuns = []
groups.append(extdata)


gWL = batchJob.jobGroup('WLonly')
WLdata = [batchJob.dataSet(WLonly), batchJob.dataSet(WLonlyHeymans)]
gWL.datasets = copy.deepcopy(WLdata)
for d in copy.deepcopy(WLdata):
    d.add(BAO, BAOdata)
    gWL.datasets.append(d)
for d in copy.deepcopy(WLdata):
    d.add(HST, HSTdata)
    gWL.datasets.append(d)
gWL.params = [[], ['mnu']]
gWL.importanceRuns = []
groups.append(gWL)


for g in groups:
    for p in g.params:
        if 'nnu' in p:
            g.importanceRuns.append(post_abundance)
            break


gWMAP = batchJob.jobGroup('WMAP')
gWMAP.params = [[]]
gWMAP.datasets = [WMAP9]
groups.append(gWMAP)

gchecks = batchJob.jobGroup('checks')
gchecks.datasets = likechecks
for d in gchecks.datasets:
    d.add(lowTEB)

gchecks.params = [[], ['mnu'], ['nnu'], ['Alens'], ['yhe']]
gchecks.importanceRuns = []
groups.append(gchecks)

skip = []

covWithoutNameOrder = [HST, 'JLA', BAORSD, 'WL', 'lensing', 'BAO', 'reion']
covNameMappings = {HSTdata:'HST', 'CamSpecHM':'CamSpec', 'CamSpecDS':'CamSpec', 'plikHM':'plik', 'plikDS':'plik', 'Mspec':'CamSpec',
                   WLonlyHeymans: WLonly}

# try to match run to exisitng covmat
covrenames = []
for planck in planck_vars:
    covrenames.append([planck, 'planck'])
covrenames.append(['planck', 'planck_CAMspec'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['Alensf', 'Alens'])
covrenames.append(['_Aphiphi', ''])
covrenames.append(['_r', ''])
covrenames.append(['_w', ''])
covrenames.append(['_alpha1', ''])
covrenames.append(['_WLonly', '_lensonly'])
covrenames.append(['lowl', 'lowTEB'])

def covRenamer(name):
    renamed = re.sub(r'_v.*_highL', '_planck_lowl_lowLike_highL', name, re.I)
    if 'wa_' in name:
        renamed = re.sub(r'_CamSpec.*', '_planck_lowl_lowLike_BAO', renamed, re.I)
    else:
        renamed = re.sub(r'_CamSpec.*', '_planck_lowl_lowLike', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]

