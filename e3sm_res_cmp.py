import os

#--------------------------------
# Functions for E3SM resolution comparsion
#--------------------------------

def load_paths_ocn(climo_ys=41, climo_ye=50):
    """Load paths of ocean data
    :climo_ys: starting year of climatology
    :climo_ye: ending year of climatology
    :returns: paths for restart file, monthly mean file, processed climotology data and output figures

    """
    # get hostname
    hostname = os.uname()[1]
    print('Running on machine {}'.format(hostname))
    # set paths
    if 'theta' in hostname:
        yshift = 5
        data_root = '/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG'
        climo_root = '/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_climo/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/high_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/run'
    elif 'edison' in hostname:
        yshift = 0
        data_root = '/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG'
        climo_root = '/global/project/projectdirs/acme/qingli/e3sm_climo/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/{:04d}-{:04d}/ocn'.format(climo_ys+yshift, climo_ye+yshift)
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res/{:04d}-{:04d}'.format(climo_ys+yshift, climo_ye+yshift)
        rst_root = data_root+'/run'
        mon_root = data_root+'/archive/ocn/hist'
    elif 'pn1803144' in hostname:
        # for testing
        data_root = os.environ['HOME']+'/data/mpas/test'
        climo_root = data_root+'/climo'
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/low_res'
        rst_root = data_root
        mon_root = data_root
    else:
        raise EnvironmentError('This script should be executed on edison, theta or pn1803144')
    os.makedirs(fig_root, exist_ok=True)
    return rst_root, mon_root, climo_root, fig_root
