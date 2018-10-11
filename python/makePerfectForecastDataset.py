# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/MyForecast/test_lensedCls_exactsim.dataset

from __future__ import absolute_import
import shutil
import os
import numpy as np
from getdist import IniFile
from CMBlikes import lastTopComment, DatasetLikelihood, ClsArray


def make_forecast_cmb_dataset(input_cl_file, output_root, output_dir=None, noise_muK_arcmin_T=None,
                              noise_muK_arcmin_P=None, NoiseVar=None, ENoiseFac=2, fwhm_arcmin=None,
                              lmin=2, lmax=None, fsky=1, fields_use=None,
                              lens_recon_noise=None, cl_data_cols=''):
    """
    Make a simulated .dataset and associated files with 'data' set at the input fiducial model.

    :param input_cl_file: input fiducial CL
    :param output_root: root name for output files, e.g. 'my_sim1'
    :param output_dir: output directory
    :param noise_muK_arcmin_T: temperature noise in muK-arcmin
    :param noise_muK_arcmin_P: polarization noise in muK-arcmin
    :param NoiseVar: effective isotropic noise variance for the temperature (N_L=NoiseVar with no beam)
    :param ENoiseFac: factor by which polarization noise variance is higher (usually 2, for Planck about 4
                        as only half the detectors polarized)
    :param fwhm_arcmin: beam fwhm in arcminutes
    :param lmin: l_min
    :param lmax: l_max
    :param fsky: sky fraction
    :param fields_use: optional list of fields to restict to (e.g. 'T E')
    :param lens_recon_noise: optional array, starting at L=0, for the PP lensing reconstruction noise, in [L(L+1)]^2C_L^phi/2pi units
    :param cl_data_cols: if not specified in file header, order of columns in input CL file (e.g. 'TT TE EE BB PP')
    :return:
    """

    use_lensing = lens_recon_noise
    use_CMB = noise_muK_arcmin_T or NoiseVar is not None

    ini = IniFile()
    dataset = ini.params

    if not cl_data_cols:
        cl_data_cols = lastTopComment(input_cl_file)
        if not cl_data_cols:
            raise Exception('input CL file must specific names of columns (TT TE EE..)')
    else:
        dataset['cl_hat_order'] = cl_data_cols

    if use_CMB:
        if NoiseVar is None:
            if noise_muK_arcmin_T is None:
                raise ValueError('Must specify noise')
            NoiseVar = (noise_muK_arcmin_T * np.pi / 180 / 60.) ** 2
            if noise_muK_arcmin_P is not None:
                ENoiseFac = (noise_muK_arcmin_P / noise_muK_arcmin_T) ** 2
        elif noise_muK_arcmin_T is not None or noise_muK_arcmin_P is not None:
            raise ValueError('Specific either noise_muK_arcmin or NoiseVar')
        if not fields_use:
            fields_use = ''
            if 'TT' or 'TE' in cl_data_cols: fields_use = 'T'
            if 'EE' or 'TE' in cl_data_cols: fields_use += ' E'
            if 'BB' in cl_data_cols: fields_use += ' B'
            if 'PP' in cl_data_cols and use_lensing: fields_use += ' P'
    else:
        fields_use = fields_use or 'P'

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', output_root)
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    dataset['fields_use'] = fields_use

    if use_CMB:
        fwhm = fwhm_arcmin / 60
        xlc = 180 * np.sqrt(8. * np.log(2.)) / np.pi
        sigma2 = (fwhm / xlc) ** 2
        noise_cols = 'TT           EE          BB'
        if use_lensing: noise_cols += '          PP'
    elif use_lensing:
        noise_cols = 'PP'
    noise_file = output_root + '_Noise.dat'
    with open(os.path.join(output_dir, noise_file), 'w') as f:
        f.write('#L %s\n' % noise_cols)

        for l in range(lmin, lmax + 1):
            NoiseCl = l * (l + 1.) / 2 / np.pi * NoiseVar * np.exp(l * (l + 1) * sigma2)
            noises = []
            if use_CMB: noises += [NoiseCl, ENoiseFac * NoiseCl, ENoiseFac * NoiseCl]
            if use_lensing: noises += [lens_recon_noise[l]]
            f.write("%d " % l + " ".join("%E" % elem for elem in noises) + "\n")

    dataset['fullsky_exact_fksy'] = fsky
    dataset['dataset_format'] = 'CMBLike2'
    dataset['like_approx'] = 'exact'

    dataset['cl_lmin'] = lmin
    dataset['cl_lmax'] = lmax

    dataset['binned'] = False

    dataset['cl_hat_includes_noise'] = False

    shutil.copy(input_cl_file, os.path.join(output_dir, output_root + '.dat'))
    dataset['cl_hat_file'] = output_root + '.dat'
    dataset['cl_noise_file '] = noise_file

    ini.saveFile(os.path.join(output_dir, output_root + '.dataset'))


if __name__ == "__main__":
    import tempfile

    # Edit parameters you want to change here
    lensedTotClFileRoot = os.path.join(os.path.dirname(__file__), '..', 'data',
                                       'base_plikHM_TT_lowTEB.minimum.theory_cl')
    # these numbers are Planck-like
    # Noise var is N_l in muK^2 for white noise
    # note  NoiseVar = (muKArcmin * np.pi / 180 / 60.) ** 2
    # Pol noise var = ENoiseFac * NoiseVar
    # 2 normally, but for Planck only half detectors are polarized
    output_dir = tempfile.gettempdir()
    output_root = 'test_sim'
    make_forecast_cmb_dataset(input_cl_file=lensedTotClFileRoot, output_root=output_root,
                              output_dir=output_dir, lmin=2, lmax=2500,
                              fwhm_arcmin=5, fsky=0.7, NoiseVar=4e-5, ENoiseFac=4)
    print('Made ' + os.path.join(output_dir, output_root + '.dataset'))

    # The rest is just a test on files produced above
    like = DatasetLikelihood(os.path.join(output_dir, output_root + '.dataset'))
    cls = ClsArray(lensedTotClFileRoot)
    cls.cls_array[0, 0] *= 1.004
    cls.cls_array[1, 1] *= 0.991
    import time

    start = time.time()
    chi2 = like.chi_squared(cls)
    end = time.time() - start
    print('Test chi2 = %s' % chi2)
    print('Time: %s' % end)
    if not np.allclose(49.055, chi2, rtol=1e-5): raise Exception('likelihood test failed')
