# provisional updated (reduced) grid settings for full mission
# New BBN fitting built in
# New BAO -> DR11
# To add: w/w0, nrun r
import batchJob, copy, re

ini_dir = 'batch2/'
cov_dir = 'planck_covmats/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
lensonly = 'lensonly'

highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST70p6'
JLA = 'JLA'

BAOdata = 'BAODR11.ini'
HSTdata = 'HST_GPE70p6.ini'

Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
lowTEB = 'lowTEB'
# dataset names
tauprior = {'prior[tau]':'0.07 0.02'}
tauname = 'tau07'
WMAPtau = {'prior[tau]':'0.09 0.013'}


camspec_detsets = ['nonclik_v97F.ini']
camspec_CS = ['nonclik_v97CS.ini']


variant_tag = ['TTTEEE', 'TT']
variant_pol_tag = ['TE', 'EE']
variants = variant_tag

planck_highL_sets = []
planck_pol_sets = []
planck_vars = ['v97CS']
planck_ini = ['CAMspec_%s.ini']
planck_base = [camspec_CS]
planck_covmats = [None]

if True:
    planck_vars += ['plik']
    planck_ini += ['plik_dx11c_%s_v13.ini']
    planck_base += [[]]
    planck_covmats += ['planck_covmats/plik_dx11c_%s_v13.covmat']
for planck, ini, base, cov in zip(planck_vars, planck_ini, planck_base, planck_covmats):
    for name, var in zip(variant_tag, variants):
        if cov is not None:
            covM = cov % (var)
        else: covM = None
        planck_highL_sets.append(batchJob.dataSet([planck , name], base + [ ini % (var)], covmat=covM))
    for var in variant_pol_tag:
        if cov is not None:
            covM = cov % (var)
        else: covM = None
        planck_pol_sets.append(batchJob.dataSet([planck , var], base + [ ini % (var)], covmat=covM))


WMAP9 = [[WMAP], ['WMAP.ini']]

start_at_bestfit = False
newCovmats = False

# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return planck in jobItem.data_set.names and (not'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))

class importanceFilterNotOmegak:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set)


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_BAO = [[BAO], [BAOdata], importanceFilterNotOmegak()]
post_HST = [[HST], [HSTdata], importanceFilterNotOmegak()]
post_JLA = [[JLA], ['JLA_marge.ini'], importanceFilterNotOmegak()]
post_nonBAO = [[HST, JLA], [HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_nonCMB = [[BAO, HST, JLA], [BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_all = [[lensing, BAO, HST, JLA], [lensing, BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegak()]
post_WP = [[ 'WMAPtau'], [WMAPtau, {'redo_no_new_data':'T'}]]

# set up groups of parameters and data sets

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(planck_highL_sets)
for d in g.datasets:
    d.add(lowTEB)

g.params = [[], ['omegak'], ['mnu'], ['r'], ['nrun', 'r'], ['nnu'], ['nrun'], ['Alens'], ['yhe'], ['w']]
g.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST, post_all]
groups.append(g)


gpol = batchJob.jobGroup('mainpol')
gpol.datasets = copy.deepcopy(planck_pol_sets)
for d in gpol.datasets:
    d.add(lowTEB)

gpol.params = [[], ['mnu'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]
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
for d in copy.deepcopy(g5.datasets):
    d.add(lensing)
    g5.datasets.append(d)

g5.importanceRuns = [post_BAO, post_nonCMB, post_WP]
groups.append(g5)


gpolnopoltau = batchJob.jobGroup('polnopoltau')
gpolnopoltau.params = [[]]
gpolnopoltau.datasets = copy.deepcopy(planck_pol_sets)
for d in copy.deepcopy(planck_pol_sets):
    d.add(lensing)
    gpolnopoltau.datasets.append(d)

gpolnopoltau.importanceRuns = [post_BAO, post_nonCMB, post_WP]
groups.append(gpolnopoltau)


g6 = batchJob.jobGroup('lensing')
g6.datasets = copy.deepcopy(g.datasets)
for d in g6.datasets:
    d.add(lensing)
    d.add(None, {'redo_theory':'F'})

g6.params = [['omegak'], ['mnu'], ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['Alens']]
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
glens.params = [['nnu'], ['omegak'], ['nnu', 'meffsterile']]
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


skip = []

covNameMappings = {HSTdata:'HST', 'v97CS':'CamSpec'}

# try to match run to exisitng covmat
covrenames = []
for planck in planck_vars:
    covrenames.append([planck, 'planck'])
covrenames.append(['planck', 'planck_CAMspec'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['_lensing', '_post_lensing'])
covrenames.append(['_BAO', '_post_BAO'])
covrenames.append(['_HST', '_post_HST'])
covrenames.append(['_JLA', '_post_SNLS'])

def covRenamer(name):
    renamed = re.sub(r'_v.*_highL', '_planck_lowl_lowLike_highL', name, re.I)
    renamed = re.sub(r'_v.*', '_planck_lowl_lowLike', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]

