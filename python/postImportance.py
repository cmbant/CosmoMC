import os, sys, batchJob, iniFile
import planckStyle as s

datafiles = ['plenslike_dx9nom_143_217_857m_g30_co_simnfg_simbeam_simlens_mono_pp_lmin_10.dat',
      'plenslike_dx9nom_143_217_857m_g30_co_simnfg_simbeam_simlens_mono_pp_lmax_2020.dat',
      'plenslike_dx9nom_143_217_857m_g30_co_simnfg_simbeam_simlens_mono_pp_lmin_10_lmax_2020.dat']

offsets = [0, 0, 20]

tags = ['lowbin', 'highbin', 'lowhighbin']

names = ['base_planck_lowl_lowLike_highL', 'base_omegak_planck_lowl_lowLike_highL', 'base_mnu_planck_lowl_lowLike_highL']

impjob = [True, False, False]
nums = [2, 2, 1]

def runJobs():
    for name, imp in zip(names, impjob):
        if imp: root = 'main/postIniFiles/' + name + '_post_lensing.ini'
        else: root = 'main/iniFiles/' + name + '_lensing.ini'
        for tag, data, offset, num in zip(tags, datafiles, offsets, nums):
            ini = iniFile.iniFile(root, keep_includes=True)
            if imp:
                improot = ini.params['redo_outroot']
            ini = iniFile.iniFile()
            ini.defaults.append(root)
            if imp: ini.params['file_root'] = improot
            outroot = name + '_' + tag
            ini.params['redo_outroot'] = 'chains/' + outroot
            ini.params['redo_theory'] = False
            ini.params['redo_no_new_data'] = True
            ini.params['redo_add'] = True
            ini.params['feedback'] = 2
            ini.params['action'] = 1
            ini.params['clik_data_lensing'] = '%DATASETDIR%clik/' + data
            if imp: ini.params['redo_skip'] = 0
            else: ini.params['redo_skip'] = 0.3
            ini.params['redo_likeoffset'] = offset
            ini.saveFile('lensing_jobs/' + outroot + '.ini')
            command = 'perl runMPI_HPCS.pl ' + 'lensing_jobs/' + outroot + ' ' + str(num)
            print 'Submitting...' + command
            os.system(command)


def runGetDist():
    for name in names:
        for tag in tags:
            ini = iniFile.iniFile()
            outroot = name + '_' + tag
            ini.params['file_root'] = 'chains/' + outroot
            ini.defaults.append('batch1/getdist_common.ini')
            ini.params['out_dir'] = 'lensing_jobs/dist/'
            ini.params['plot_data_dir'] = 'main/plot_data/'
            ini.params['ignore_rows'] = 0
            fname = 'lensing_jobs/' + outroot + '_getdist.ini'
            ini.saveFile(fname)
            os.system('./getdist ' + fname)

def makePlots():
    g = s.plotter
    for name, imp in zip(names, impjob):
        if imp:root = name + '_post_lensing'
        else: root = name + '_lensing'
        roots = [root] + [ name + '_' + tag for tag in tags]
        g.plots_1d(roots, legend_labels=['reference'] + tags, nx=5, legend_ncol=4)
        g.export('plots/lens_compare_' + name + '.pdf')
        g.newPlot()


# runGetDist()
makePlots()
