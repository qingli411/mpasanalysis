#!/usr/bin/env python

from mpasanalysis import *
import e3sm_res_cmp

def main():
    global fig_dir
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth

    # get paths of restart files, monthly mean output files, processed climatology files and output figures
    ts_ys = 1
    ts_ye = 20
    plt_ys = 1
    plt_ye = 20
    nmon = 1 # 12 for production and 1 for testing
    data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye,
                                            ts_ys=ts_ys, ts_ye=ts_ye, runname='gl-mesh')
    rst_root = data_root['rst_root']
    mon_root = data_root['mon_root']
    fig_root = data_root['fig_root']
    rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)

    # load dataset
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

    # Salinity (psu)

    varname = 'salinityAtSurface'
    units = 'psu'
    levels = np.linspace(28, 36, 41)

    fig_dir = fig_root+'/Animation/highfreq/'+varname
    os.makedirs(fig_dir, exist_ok=True)
    for y in np.arange(plt_ye-plt_ys+1)+1:
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
    global fig_dir
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
                            name=vname, units=units)
        # plot figure: map
        fig = plt.figure(figsize=[6, 5.5])
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


