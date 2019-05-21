# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from mpasanalysis import *
import e3sm_res_cmp
# %matplotlib inline

# get paths of restart files, monthly mean output files, processed climatology files and output figures
ts_ys = 1
ts_ye = 25
data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye)
rst_root = data_root['rst_root']
ts_root = data_root['ts_root']
climo_root = data_root['climo_root']
fig_root = data_root['fig_root']
rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)
data_root_ice = e3sm_res_cmp.load_paths_ice(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye)
mon_root_ice = data_root_ice['mon_root']

# load dataset
f_rst = Dataset(rst_file, 'r')

# read grid information
lon = np.degrees(f_rst.variables['lonCell'][:])
lat = np.degrees(f_rst.variables['latCell'][:])
cellArea = f_rst.variables['areaCell'][:]


# ### March MLD

# +
climo_file = climo_root+'/mixedLayerDepth/mpaso_03_climo.nc'
f_climo = Dataset(climo_file, 'r')
levels = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                   110, 130, 150, 180, 210, 240, 280, 320, 360,
                   407, 454, 500, 1000, 1500, 2000])

fig = plt.figure(figsize=[6,6])
mld_d = f_climo.variables['timeMonthly_avg_dThreshMLD'][0,:]
mpaso_mld_d = MPASOMap(data=mld_d, lat=lat, lon=lon, cellarea=cellArea,
                       name='MLD (density threshold)', units='m')
m = mpaso_mld_d.plot(region='LabSea', ptype='contourf', levels=levels, label=calendar.month_abbr[3])


climo_file_ice = mon_root_ice+'/mpascice.hist.am.timeSeriesStatsMonthly.0041-01-01.nc'
f_climo_ice = Dataset(climo_file_ice, 'r')
icearea_d = f_climo_ice.variables['timeMonthly_avg_iceAreaCell'][0,:]
mpascice_icearea_d = MPASCICEMap(data=icearea_d, lat=lat, lon=lon, cellarea=cellArea,
                                 name='Ice Area', units='none')
mpascice_icearea_d.overlay(m, levels=[0.15, 0.85], cmap='bone_r')

robj = region_latlon('LabSea_SD1')
m.drawgreatcircle(robj.lon_ll, robj.lat_ll, robj.lon_ll, robj.lat_ur, color='y')
m.drawgreatcircle(robj.lon_ll, robj.lat_ur, robj.lon_ur, robj.lat_ur, color='y')
m.drawgreatcircle(robj.lon_ur, robj.lat_ur, robj.lon_ur, robj.lat_ll, color='y')
m.drawgreatcircle(robj.lon_ur, robj.lat_ll, robj.lon_ll, robj.lat_ll, color='y')
plt.show()

# save data
path_mld = fig_root+'/data_LabSea_climo_03_dThreshMLD.npz'
mpaso_mld_d.save(path_mld)
path_ice = fig_root+'/data_LabSea_climo_03_iceArea.npz'
mpascice_icearea_d.save(path_ice)

figname = fig_root+'/LabSea_climo_03_dThreshMLD_ice.png'
# fig.savefig(figname, dpi = 300)
# plt.close(fig)
# -

path_mld = fig_root+'/data_LabSea_climo_03_dThreshMLD.npz'
path_ice = fig_root+'/data_LabSea_climo_03_iceArea.npz'
mpaso_mld_d1 = MPASOMap().load(path_mld)
m1 = mpaso_mld_d1.plot(region='LabSea', ptype='contourf', levels=levels, label=calendar.month_abbr[3])
mpascice_icearea_d1 = MPASCICEMap().load(path_ice)
mpascice_icearea_d.overlay(m1, levels=[0.15, 0.85], cmap='bone_r')

# +
grpname = 'mixedLayerDepth'
region = 'LabSea_SD1'
varlist = ['timeMonthly_avg_dThreshMLD',
           'timeMonthly_avg_tThreshMLD']
nvar = len(varlist)
for j, varname in enumerate(varlist):
    ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(grpname, varname, ts_ys, ts_ye)
    f_ts = Dataset(ts_file, 'r')
    ncvar = f_ts.variables[varname]
    data = ncvar[:]
    if j == 0:
        nt = data.shape[0]
        time = np.linspace(1,nt,nt)/12.
        mdata = np.zeros([nvar, nt])
    for i in np.arange(nt):
        mpaso_obj = MPASOMap(data=data[i,:], lat=lat, lon=lon, cellarea=cellArea, name=varname, units=ncvar.units)
        mdata[j, i] = mpaso_obj.mean(region=region)

print(mdata.shape)


# +
fig, axarr = plt.subplots(nvar, sharex='col')
fig.set_size_inches(8, 8)
for i in np.arange(nvar):
    axarr[i].plot(time, mdata[i, :])
    axarr[i].set_title(varlist[i][16:])
    axarr[i].set_ylabel('MLD (m)')
    if i == nvar-1:
        axarr[i].set_xlabel('Time (Year)')

# plt.show()
# figname = fig_root+'/LabSea_ts_MLD.png'
# fig.savefig(figname, dpi = 300)
# plt.close(fig)
# -

# ### Heat Flux (W m$^{-2}$)

# +
grpname = 'heatFlux'
region = 'LabSea_SD1'
varlist = ['timeMonthly_avg_sensibleHeatFlux',
           'timeMonthly_avg_latentHeatFlux',
           'timeMonthly_avg_shortWaveHeatFlux',
           'timeMonthly_avg_longWaveHeatFluxDown',
           'timeMonthly_avg_longWaveHeatFluxUp']
nvar = len(varlist)
for j, varname in enumerate(varlist):
    ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(grpname, varname, ts_ys, ts_ye)
    f_ts = Dataset(ts_file, 'r')
    ncvar = f_ts.variables[varname]
    data = ncvar[:]
    if j == 0:
        nt = data.shape[0]
        time = np.linspace(1,nt,nt)/12.
        mdata = np.zeros([nvar, nt])
    for i in np.arange(nt):
        mpaso_obj = MPASOMap(data=data[i,:], lat=lat, lon=lon, cellarea=cellArea, name=varname, units=ncvar.units)
        mdata[j, i] = mpaso_obj.mean(region=region)

print(mdata.shape)


# +
fig, axarr = plt.subplots(5, sharex='col')
fig.set_size_inches(8, 12)
for i in np.arange(5):
    if i < 3:
        axarr[i].plot(time, mdata[i, :])
        axarr[i].set_title(varlist[i][16:])
    elif i == 3:
        axarr[i].plot(time, np.sum(mdata[3:5, :], axis=0))
        axarr[i].set_title('longWaveHeatFlux')
    else:
        axarr[i].plot(time, np.sum(mdata, axis=0))
        axarr[i].set_title('totalHeatFlux')
    axarr[i].set_ylabel('Heat Flux (W m$^{-2}$)')
    if i == nvar-1:
        axarr[i].set_xlabel('Time (Year)')


# plt.show()
# figname = fig_root+'/LabSea_ts_heatFlux.png'
# fig.savefig(figname, dpi = 300)
# plt.close(fig)
# -

# ### Fresh Water Flux (kg m$^{-2}$ s$^{-1}$)

# +
grpname = 'freshWaterFlux'
region = 'LabSea_SD1'
varlist = ['timeMonthly_avg_evaporationFlux',
           'timeMonthly_avg_rainFlux',
           'timeMonthly_avg_snowFlux',
           'timeMonthly_avg_seaIceSalinityFlux',
           'timeMonthly_avg_seaIceFreshWaterFlux',
           'timeMonthly_avg_riverRunoffFlux',
           'timeMonthly_avg_iceRunoffFlux']
nvar = len(varlist)
for j, varname in enumerate(varlist):
    ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(grpname, varname, ts_ys, ts_ye)
    f_ts = Dataset(ts_file, 'r')
    ncvar = f_ts.variables[varname]
    data = ncvar[:]
    if j == 0:
        nt = data.shape[0]
        time = np.linspace(1,nt,nt)/12.
        mdata = np.zeros([nvar, nt])
    for i in np.arange(nt):
        mpaso_obj = MPASOMap(data=data[i,:], lat=lat, lon=lon, cellarea=cellArea, name=varname, units=ncvar.units)
        mdata[j, i] = mpaso_obj.mean(region=region)

print(mdata.shape)


# +
fig, axarr = plt.subplots(5, sharex='col')
fig.set_size_inches(8, 12)
for i in np.arange(5):
    if i == 0:
        axarr[i].plot(time, mdata[i, :])
        axarr[i].set_title(varlist[0][16:])
    elif i == 1:
        axarr[i].plot(time, np.sum(mdata[1:3, :], axis=0))
        axarr[i].set_title('precipitationFlux')
    elif i == 2:
        axarr[i].plot(time, np.sum(mdata[3:5, :], axis=0))
        axarr[i].set_title('seaIceFlux')
    elif i == 3:
        axarr[i].plot(time, np.sum(mdata[5:7, :], axis=0))
        axarr[i].set_title('runoffFlux')
    else:
        axarr[i].plot(time, np.sum(mdata,axis=0))
        axarr[i].set_title('totalFlux')
    axarr[i].set_ylabel('Heat Flux (kg m$^{-2}$ s$^{-1}$)')
    if i == nvar-1:
        axarr[i].set_xlabel('Time (Year)')

# plt.show()
# figname = fig_root+'/LabSea_ts_heatFlux.png'
# fig.savefig(figname, dpi = 300)
# plt.close(fig)
# -


