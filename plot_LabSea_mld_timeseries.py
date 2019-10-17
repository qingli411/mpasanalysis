# get time series of mean mld
from mpasanalysis import *
import e3sm_res_cmp

# get paths of restart files, monthly mean output files, processed climatology files and output figures
ts_ys = 1
ts_ye = 25
data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye)
rst_root = data_root['rst_root']
ts_root = data_root['ts_root']
fig_root = data_root['fig_root']
rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)

# load dataset
f_rst = Dataset(rst_file, 'r')

# +
# read grid information
lon = np.degrees(f_rst.variables['lonCell'][:])
lat = np.degrees(f_rst.variables['latCell'][:])
cellArea = f_rst.variables['areaCell'][:]

# Mixed Layer Depth (m)

grpname = 'mixedLayerDepth'
region = 'LabSea_SD1'
varlist = ['timeMonthly_avg_dThreshMLD',
           'timeMonthly_avg_tThreshMLD']
nvar = len(varlist)
for j, varname in enumerate(varlist):
    print('Saving time series of \'{}\'...'.format(varname))
    ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(grpname, varname, ts_ys, ts_ye)
    f_ts = Dataset(ts_file, 'r')
    ncvar = f_ts.variables[varname]
    data = ncvar[:]
    if j == 0:
        nt = data.shape[0]
        time = np.linspace(1,nt,nt)/12+ts_ys
        mdata = np.zeros([nvar, nt])
        npcount = np.floor(nt/10)
    for i in np.arange(nt):
        if np.mod(i, npcount) == 0:
            print('{:6.1f} %'.format(i/nt*100.0))
        mpaso_obj = MPASOMap(data=data[i,:], lat=lat, lon=lon, cellarea=cellArea, name=varname, units=ncvar.units)
        mdata[j, i] = mpaso_obj.mean(region=region)
    print('{:6.1f} %'.format(100.0))

print('Time dimension: {}'.format(mdata.shape[1]))
data_name = fig_root+'/ts_mld.npz'
np.savez(data_name, time=time, dThreshMLD=mdata[0,:], tThreshMLD=mdata[1,:])

