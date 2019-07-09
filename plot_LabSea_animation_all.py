#!/usr/bin/env python

from mpasanalysis import *
import e3sm_res_cmp

def main():
    global fig_dir
    global lon, lat, refMidDepth, cellArea, refLayerThickness, bottomDepth
    global s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat
    global s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat

    # get paths of restart files, monthly mean output files, processed climatology files and output figures
    ts_ys = 21
    ts_ye = 50
    plt_ys = 21
    plt_ye = 50
    nmon = 12 # 12 for production and 1 for testing
    runname = 'gl-mesh-gm1800'
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

    # ### Cross sections

    # cross section 1
    s1_s_lat = 63
    s1_s_lon = 296
    s1_e_lat = 50
    s1_e_lon = 320
    # cross section 2 (WOCE AR7W)
    s2_s_lat = 53.5
    s2_s_lon = 304.5
    s2_e_lat = 61
    s2_e_lon = 312

    # ### Salinity (psu) with ice fraction (unitless)

    varname = 'timeMonthly_avg_activeTracers_salinity'
    units = 'psu'
    levels = np.linspace(28, 36, 41)
    varname1 = 'timeMonthly_avg_iceAreaCell'
    levels1 = [0.15, 0.85]
    units1 = 'unitless'

    fig_dir = fig_root+'/Animation/'+varname
    os.makedirs(fig_dir, exist_ok=True)
    for y in np.arange(plt_ys, plt_ye):
        for m in np.arange(nmon)+1:
            print('{:04d}-{:02d}'.format(y, m))
            mon_file = mon_root+'/mpaso.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(y, m)
            f_mon = Dataset(mon_file, 'r')
            mon_file1 = mon_root_ice+'/mpascice.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(y, m)
            f_mon1 = Dataset(mon_file1, 'r')
            plot_labsea_ocn_ice(f_mon, f_mon1, varname, varname1, units, units1, levels, levels1, y, m)
            f_mon.close()
            f_mon1.close()

    # ### Temperature (degC) with ice thickness (m)

    varname = 'timeMonthly_avg_activeTracers_temperature'
    units = 'degC'
    levels = np.linspace(-2, 26, 57)
    varname1 = 'timeMonthly_avg_iceVolumeCell'
    levels1 = np.linspace(0.5, 5, 10)
    units1 = 'm'

    fig_dir = fig_root+'/Animation/'+varname
    os.makedirs(fig_dir, exist_ok=True)
    for y in np.arange(plt_ys, plt_ye):
        for m in np.arange(nmon)+1:
            print('{:04d}-{:02d}'.format(y, m))
            mon_file = mon_root+'/mpaso.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(y, m)
            f_mon = Dataset(mon_file, 'r')
            mon_file1 = mon_root_ice+'/mpascice.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(y, m)
            f_mon1 = Dataset(mon_file1, 'r')
            plot_labsea_ocn_ice(f_mon, f_mon1, varname, varname1, units, units1, levels, levels1, y, m)
            f_mon.close()
            f_mon1.close()

def plot_labsea_ocn_ice(f_ocn, f_ice, vname_ocn, vname_ice, units_ocn, units_ice, levels_ocn, levels_ice, iyear, imon):
    """Plot map.

    :f_ocn: (netcdf4 Dataset) ocean input file
    :f_ice: (netcdf4 Dataset) sea ice input file
    :vname_ocn: (str) ocean variable name
    :vname_ice: (str) sea ice variable name
    :units_ocn: (str) ocean variable units
    :units_ice: (str) sea ice variable units
    :levels_ocn: (list) levels for ocean contour
    :levels_ice: (list) levels for sea ice contour
    :iyear: (int) year index
    :imon: (int) month index
    :returns: TODO

    """
    global fig_dir
    global lon, lat, refMidDepth, cellArea
    global s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat
    global s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat

    # read monthly mean ocean data
    data_ocn = f_ocn.variables[vname_ocn][0,:,:]
    mpasovol_obj = MPASOVolume(data=data_ocn, lat=lat, lon=lon, depth=refMidDepth, cellarea=cellArea,
                               layerthickness=refLayerThickness, bottomdepth=bottomDepth,
                               name=vname_ocn, units=units_ocn)
    mpaso_obj = mpasovol_obj.get_map(depth=0.0)
    # read monthly mean sea ice data
    data_ice = f_ice.variables[vname_ice][0,:]
    mpascice_obj = MPASCICEMap(data=data_ice, lat=lat, lon=lon, cellarea=cellArea,
                               name=vname_ice, units=units_ice)
    # year and month
    yyyy = '{:04d}'.format(iyear)
    mm = '{:02d}'.format(imon)

    # plot figure 1: map
    fig = plt.figure(figsize=[6, 5.5])
    m,tmp = mpaso_obj.plot(region='LabSea', levels=levels_ocn, ptype='contourf')
    mpascice_obj.overlay(m, levels=levels_ice, cmap='bone_r')
    m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='gray')
    m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='gray')
    axis = plt.gca()
    axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                 fontsize=12, color='k', va='top',
                 bbox=dict(boxstyle='square',ec='k',fc='w'))
    figname = fig_dir+'/LabSea_climo_Map_'+yyyy+'-'+mm+'.png'
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

    # plot figure 2: vertical cross section 1
    fig = plt.figure(figsize=[6, 4])
    mpaso_vcsec1 = mpasovol_obj.get_vertical_cross_section(lon0=s1_s_lon, lat0=s1_s_lat,
                                                           lon1=s1_e_lon, lat1=s1_e_lat)
    mpaso_vcsec1.plot(levels=levels_ocn, depth_mode='native')
    axis = plt.gca()
    axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                 fontsize=12, color='k', va='top',
                 bbox=dict(boxstyle='square',ec='k',fc='w'))
    figname = fig_dir+'/LabSea_climo_VCSec1_'+yyyy+'-'+mm+'.png'
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

    # figure 3: vertical cross section 2
    fig = plt.figure(figsize=[6, 4])
    mpaso_vcsec2 = mpasovol_obj.get_vertical_cross_section(lon0=s2_s_lon, lat0=s2_s_lat,
                                                         lon1=s2_e_lon, lat1=s2_e_lat, depth_bottom=4500)
    mpaso_vcsec2.plot(levels=levels_ocn, depth_mode='native')
    axis = plt.gca()
    axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                 fontsize=12, color='k', va='top',
                 bbox=dict(boxstyle='square',ec='k',fc='w'))
    figname = fig_dir+'/LabSea_climo_VCSec2_'+yyyy+'-'+mm+'.png'
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

if __name__ == "__main__":
    main()

