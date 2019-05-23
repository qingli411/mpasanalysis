#!/usr/bin/env python
from mpasanalysis import *
import e3sm_res_cmp

def main():
    global fig_dir
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth
    global s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat
    global s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat

    # get paths of restart files, monthly mean output files, processed climatology files and output figures
    ts_ys = 1
    ts_ye = 20
    plt_ys = 1
    plt_ye = 20
    nmon = 12 # 12 for production and 1 for testing
    runname = 'gl-mesh'
    data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye, runname=runname)
    rst_root = data_root['rst_root']
    mon_root = data_root['mon_root']
    fig_root = data_root['fig_root']
    rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)
    data_root_ice = e3sm_res_cmp.load_paths_ice(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye, runname=runname)
    mon_root_ice = data_root_ice['mon_root']

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

    # transect
    trnsct = transect('AR7W')

    # velocity
    varname_x = 'timeMonthly_avg_velocityZonal'
    varname_y = 'timeMonthly_avg_velocityMeridional'
    varname = 'timeMonthly_avg_velocity'
    units = 'm/s'
    levels = np.linspace(-0.3, 0.3, 21)

    fig_dir = fig_root+'/Animation/'+varname
    os.makedirs(fig_dir, exist_ok=True)
    for y in np.arange(plt_ye-plt_ys)+1:
        for m in np.arange(nmon)+1:
            print('{:04d}-{:02d}'.format(y, m))
            mon_file = mon_root+'/mpaso.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(y, m)
            print(mon_file)
            f_mon = Dataset(mon_file, 'r')
            plot_LabSea_velocity(f_mon, varname_x, varname_y, varname, units, levels, y, m, trnsct)
            f_mon.close()

def plot_LabSea_velocity(f_in, vname_x, vname_y, vname, units, levels, iyear, imon, trnsct):
    """Plot map.

    :f_in: (netcdf4 Dataset) ocean input file
    :vname_x: (str) variable name x-component
    :vname_y: (str) variable name y-component
    :vname: (str) variable name for display
    :units: (str) variable units for display
    :levels: (list) levels for contours
    :iyear: (int) year index
    :imon: (int) month index
    :trnsct: (VerticalTransect object) transect
    :returns: none

    """
    global fig_dir
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth

    # read monthly mean ocean data
    data_x = f_in.variables[vname_x][0,:,:]
    data_y = f_in.variables[vname_y][0,:,:]
    mpasovol_obj_x = MPASOVolume(data=data_x, lat=lat, lon=lon, depth=refMidDepth, cellarea=cellArea,
                                 layerthickness=refLayerThickness, bottomdepth=bottomDepth,
                                 name=vname_x, units=units)
    mpasovol_obj_y = MPASOVolume(data=data_y, lat=lat, lon=lon, depth=refMidDepth, cellarea=cellArea,
                                 layerthickness=refLayerThickness, bottomdepth=bottomDepth,
                                 name=vname_y, units=units)
    mpaso_obj_x = mpasovol_obj_x.get_map(depth=0.0)
    mpaso_obj_y = mpasovol_obj_y.get_map(depth=0.0)
    data = np.sqrt(mpaso_obj_x.data**2+mpaso_obj_y.data**2)
    mpaso_obj = MPASOMap(data=data, lon=mpaso_obj_x.lon, lat=mpaso_obj_x.lat, cellarea=mpaso_obj_x.cellarea,
                         name=vname, units=units)
    # year and month
    yyyy = '{:04d}'.format(iyear)
    mm = '{:02d}'.format(imon)

    # plot figure 1: map
    fig = plt.figure(figsize=[6, 5.5])
    levels1 = np.linspace(0, levels[-1], 31)
    m,tmp = mpaso_obj.plot(region='LabSea', levels=levels1, ptype='contourf')
    m.drawgreatcircle(trnsct.lon0, trnsct.lat0, trnsct.lon1, trnsct.lat1, color='gray')
    axis = plt.gca()
    axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                 fontsize=12, color='k', va='top',
                 bbox=dict(boxstyle='square',ec='k',fc='w'))
    figname = fig_dir+'/LabSea_climo_Map_'+yyyy+'-'+mm+'.png'
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

    # plot figure 2: vertical transect
    fig = plt.figure(figsize=[6, 4])
    plot_transect_normal(mpasovol_obj_x, mpasovol_obj_y, trnsct, name='Normal velocity', levels=levels,
                         depth_mode='native', cmap='RdBu_r' )
    axis = plt.gca()
    axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                 fontsize=12, color='k', va='top',
                 bbox=dict(boxstyle='square',ec='k',fc='w'))
    figname = fig_dir+'/LabSea_climo_VCSec1_'+yyyy+'-'+mm+'.png'
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


if __name__ == "__main__":
    main()

