#!/usr/bin/env python

from mpasanalysis import *
import e3sm_res_cmp

def main():
    global fig_dir, mpasmesh
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth

    # get paths of restart files, monthly mean output files, processed climatology files and output figures
    ts_ys = 46
    ts_ye = 55
    plt_ys = 46
    plt_ye = 55
    nmon = 12 # 12 for production and 1 for testing
    data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye,
                                            ts_ys=ts_ys, ts_ye=ts_ye)
    rst_root = data_root['rst_root']
    mon_root = data_root['mon_root']
    fig_root = data_root['fig_root']
    rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)

    # load dataset
    mpasmesh = MPASMesh(filepath=rst_file)
    f_rst = Dataset(rst_file, 'r')

    # read grid information
    lon = np.degrees(f_rst.variables['lonCell'][:])
    lat = np.degrees(f_rst.variables['latCell'][:])
    cellArea = f_rst.variables['areaCell'][:]
    bottomDepth = f_rst.variables['bottomDepth'][:]

    refBottomDepth = f_rst.variables['refBottomDepth'][:]
    nVertLevels = len(refBottomDepth)
    refTopDepth = np.zeros(nVertLevels)
    refTopDepth[1:nVertLevels] = refBottomDepth[0:nVertLevels-1]
    refLayerThickness = refTopDepth-refBottomDepth
    refMidDepth = 0.5*(refTopDepth+refBottomDepth)

    # Temperature (degC)

    # varname = 'temperatureAtSurface'
    # units = 'degC'
    # levels = np.linspace(-2, 26, 57)

    # varname = 'temperatureAt250m'
    # units = 'degC'
    # levels = np.linspace(2, 8, 61)

    # Salinity (psu)

    # varname = 'salinityAtSurface'
    # units = 'psu'
    # levels = np.linspace(28, 36, 81)

    varname = 'salinityAt250m'
    units = 'psu'
    levels = np.linspace(34.7, 35.5, 41)

    # MLD (m)

    # varname = 'tThreshMLD'
    # units = 'm'
    # levels = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
    #                    110, 130, 150, 180, 210, 240, 280, 320, 360,
    #                    407, 454, 500, 1000, 1500, 2000])

    # relative vorticity at 250 m (s^{-1})
    # varname = 'relativeVorticityAt250m'
    # units = 's^{-1}'
    # levels = np.linspace(-5e-5, 5e-5, 51)

    fig_dir = fig_root+'/Animation/highfreq/'+varname
    os.makedirs(fig_dir, exist_ok=True)
    for y in np.arange(plt_ye-plt_ys+1)+plt_ys:
        for m in np.arange(nmon)+1:
            print('{:04d}-{:02d}'.format(y, m))
            mon_file = mon_root+'/mpaso.hist.am.highFrequencyOutput.{:04d}-{:02d}-01_00.00.00.nc'.format(y, m)
            print(mon_file)
            f_mon = Dataset(mon_file, 'r')
            plot_labsea_highfreq(f_mon, varname, units, levels, y, m)
            f_mon.close()

def plot_labsea_highfreq(f_in, vname, units, levels, iyear, imon):
    """Plot map.

    :f_in: (netcdf4 Dataset) input file
    :vname: (str) variable name
    :units: (str) variable units
    :levels: (list) levels for contours
    :iyear: (int) year index
    :imon: (int) month index
    :returns: TODO

    """
    global fig_dir, mpasmesh
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth

    # year and month
    yyyy = '{:04d}'.format(iyear)
    mm = '{:02d}'.format(imon)

    # read monthly mean ocean data
    data_in = f_in.variables[vname][:]
    nt = data_in.shape[0]
    for i in np.arange(nt):
        ii = '{:02d}'.format(i)
        data = data_in[i,:]
        mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                             mesh=mpasmesh, name=vname, units=units)
        # plot figure: map
        fig = plt.figure(figsize=[6, 5.5])
        # m,tmp = mpaso_obj.plot(region='LabSea', levels=levels, ptype='contourf', cmap='RdBu_r')
        m,tmp = mpaso_obj.plot(region='LabSea', levels=levels, ptype='contourf')
        axis = plt.gca()
        axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        figname = fig_dir+'/LabSea_climo_Map_'+yyyy+'-'+mm+'-'+ii+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)

if __name__ == "__main__":
    main()


