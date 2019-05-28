import os

#--------------------------------
# Functions for E3SM resolution comparsion
#--------------------------------

def load_paths_ocn(climo_ys=41, climo_ye=50, ts_ys=1, ts_ye=50, runname=None):
    """Load paths of ocean data
    :climo_ys: starting year of climatology
    :climo_ye: ending year of climatology
    :ts_ys: starting year of timeseries
    :ts_ye: ending year of timeseries
    :returns: paths for restart file, monthly mean file, processed climotology data,
              timeseries data and output figures, saved in a dictionary

    """
    # get hostname
    hostname = os.uname()[1]
    print('Running on machine {}'.format(hostname))
    # set paths
    if 'theta' in hostname:
        yshift = 0
        data_root = '/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG'
        climo_root = '/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_climo/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
        ts_root = '/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_ts/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/high_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/run'
    elif 'edison' in hostname:
        yshift = 0
        data_root = '/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG'
        climo_root = '/global/project/projectdirs/acme/qingli/e3sm_climo/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
        ts_root = '/global/project/projectdirs/acme/qingli/e3sm_ts/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/archive/ocn/hist'
    elif 'blues' in hostname:
        yshift = 0
        if runname == 'low-res-noSI':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res-noSI/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'low-res-g':
            data_root = '/lcrc/group/acme/qingli/acme_scratch/anvil/GMPAS-IAF_T62_oEC60to30v3_CTRL'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/GMPAS-IAF_T62_oEC60to30v3_CTRL/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/GMPAS-IAF_T62_oEC60to30v3_CTRL/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res-g/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'gl-mesh':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/gl-mesh/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'low-res-gm6h':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo//20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil/{:04d}-{:04d}/ocn'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low-res-gm6h/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        else:
            raise ValueError('Run \'{}\' not supported'.format(runname))
    elif 'pn1803144' in hostname:
        if runname is None:
            # for testing
            data_root = os.environ['HOME']+'/data/mpas/test'
            climo_root = data_root+'/climo'
            ts_root = data_root+'/ts'
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/test'
            rst_root = data_root
            mon_root = data_root
        elif runname == 'gl-mesh':
            data_root = os.environ['HOME']+'/data/mpas/test/gl-mesh'
            climo_root = data_root+'/climo'
            ts_root = data_root+'/ts'
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/test/gl-mesh'
            rst_root = data_root
            mon_root = data_root
        else:
            raise ValueError('Run \'{}\' not supported'.format(runname))
    else:
        raise EnvironmentError('This script should be executed on edison, theta, blues or pn1803144')
    os.makedirs(fig_root, exist_ok=True)
    path = {'rst_root': rst_root,
            'mon_root': mon_root,
            'climo_root': climo_root,
            'ts_root': ts_root,
            'fig_root': fig_root}
    return path

def load_paths_ice(climo_ys=41, climo_ye=50, ts_ys=1, ts_ye=50, runname=None):
    """Load paths of sea ice data
    :climo_ys: starting year of climatology
    :climo_ye: ending year of climatology
    :ts_ys: starting year of timeseries
    :ts_ye: ending year of timeseries
    :returns: paths for restart file, monthly mean file, processed climotology data,
              timeseries data and output figures, saved in a dictionary

    """
    # get hostname
    hostname = os.uname()[1]
    print('Running on machine {}'.format(hostname))
    # set paths
    if 'theta' in hostname:
        yshift = 0
        data_root = '/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG'
        climo_root = '/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_climo/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
        ts_root = '/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_ts/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/high_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/run'
    elif 'edison' in hostname:
        yshift = 0
        data_root = '/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG'
        climo_root = '/global/project/projectdirs/acme/qingli/e3sm_climo/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
        ts_root = '/global/project/projectdirs/acme/qingli/e3sm_ts/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/archive/ice/hist'
    elif 'blues' in hostname:
        yshift = 0
        if runname == 'low-res-noSI':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res-noSI/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'low-res-g':
            data_root = '/lcrc/group/acme/qingli/acme_scratch/anvil/GMPAS-IAF_T62_oEC60to30v3_CTRL'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/GMPAS-IAF_T62_oEC60to30v3_CTRL/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/GMPAS-IAF_T62_oEC60to30v3_CTRL/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res-g/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'gl-mesh':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/gl-mesh/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif runname == 'low-res-gm6h':
            data_root = '/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil'
            climo_root = '/lcrc/group/acme/qingli/e3sm_climo//20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil/{:04d}-{:04d}/ice'.format(climo_ys+yshift, climo_ye+yshift)
            ts_root = '/lcrc/group/acme/qingli/e3sm_ts/20190326.GM600.A_WCYCL1850S.ne30_oECv3.anvil/{:04d}-{:04d}/ice'.format(ts_ys+yshift, ts_ye+yshift)
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low-res-gm6h/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        else:
            raise ValueError('Run \'{}\' not supported'.format(runname))
    elif 'pn1803144' in hostname:
        if runname is None:
            # for testing
            data_root = os.environ['HOME']+'/data/mpas/test'
            climo_root = data_root+'/climo'
            ts_root = data_root+'/ts'
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/test'
            rst_root = data_root
            mon_root = data_root
        elif runname == 'gl-mesh':
            data_root = os.environ['HOME']+'/data/mpas/test/gl-mesh'
            climo_root = data_root+'/climo'
            ts_root = data_root+'/ts'
            fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/test/gl-mesh'
            rst_root = data_root
            mon_root = data_root
        else:
            raise ValueError('Run \'{}\' not supported'.format(runname))
    else:
        raise EnvironmentError('This script should be executed on edison, theta, blues or pn1803144')
    os.makedirs(fig_root, exist_ok=True)
    path = {'rst_root': rst_root,
            'mon_root': mon_root,
            'climo_root': climo_root,
            'ts_root': ts_root,
            'fig_root': fig_root}
    return path
