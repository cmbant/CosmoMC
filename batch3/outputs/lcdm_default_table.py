import planckStyle as s
from getdist import types
import six

g = s.getSinglePlotter()

pars = ['omegabh2', 'omegach2', 'theta', 'tau', 'logA', 'ns', 'omegamh2', 'H0', 'omegam', 'age', 'sigma8', 'S8', 'zrei',
        'thetastar', 'rdrag']

lines = []
heading = ''
formatter = types.NoLineTableFormatter()


class col(object):
    def __init__(self, datatag, title, samples=None, bestfit=False):
        if datatag is not None:
            if isinstance(datatag,six.string_types): datatag = [datatag]
            samples=[]
            for tag in datatag:
                root = g.getRoot('', tag)
                samples += [g.sampleAnalyser.samplesForRoot(root)]
                samples[-1].paramNames.setLabelsAndDerivedFromParamNames(g.settings.param_names_for_labels)
        self.samples = samples
        self.title = title
        if bestfit:
            self.bestfit = samples[0].getBestFit()
            self.marge = samples[0].getMargeStats()
            self.marge.addBestFit(self.bestfit)
        else:
            self.bestfit = None


items = []

items += [col('plikHM_TTTEEE_lowl_lowE_lensing', 'Best fit', bestfit=True)]

plik = col('plikHM_TTTEEE_lowl_lowE_lensing', 'Marginalized')
items += [plik]

camspec = col('CamSpecHM_TTTEEE_lowl_lowE_lensing', '\\camsepc')
items += [camspec]

diff= col(['CamSpecHM_TTTEEE_lowl_lowE_lensing','plikHM_TTTEEE_lowl_lowE_lensing'], r'(\camspec-\plik)/$\sigma_{\rm \plik}$')
items += [diff]

joint = plik.samples[0].getCombinedSamplesWithSamples(camspec.samples[0])
items += [col(None, 'Combined', samples=[joint])]

for i, par in enumerate(pars):
    param = plik.samples[0].paramNames.parWithName(par)
    line = '$' + param.getLabel() + '$ &'
    # line = '\\hline\n' + line
    for j, item in enumerate(items):
        param = item.samples[0].paramNames.parWithName(par)
        if item.bestfit:
            latex = item.marge.texValues(formatter, param, limit=1)[1]
        else:
            if len(item.samples)>1:
                diff = (item.samples[0].mean(param)-item.samples[1].mean(param))/item.samples[1].std(param)
                latex = r"%+.1f"%diff
            else:
                latex = item.samples[0].getLatex(params=[par], limit=1)[1][0]
        if j: line += ' & '
        line += ' $ ' + latex + ' $ '

    lines.append(line)
print heading + '\\\\\n\\hline\n' + '\\\\\n'.join(lines) + '\\\\\n'
