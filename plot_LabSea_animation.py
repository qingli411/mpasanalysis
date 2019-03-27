
# coding: utf-8

# In[ ]:


from mpasanalysis import *
import e3sm_res_cmp
# get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


# get paths of restart files, monthly mean output files, processed climatology files and output figures
ts_ys = 1
ts_ye = 25
plt_ys = 1
plt_ye = 25
nmon = 12 # 12 for production and 1 for testing
data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye)
rst_root = data_root['rst_root']
ts_root = data_root['ts_root']
fig_root = data_root['fig_root']
rst_file = rst_root+'/mpaso.rst.{:04d}-01-01_00000.nc'.format(ts_ye+1)
data_root_ice = e3sm_res_cmp.load_paths_ice(climo_ys=ts_ys, climo_ye=ts_ye, ts_ys=ts_ys, ts_ye=ts_ye)
ts_root_ice = data_root_ice['ts_root']


# In[ ]:


# load dataset
f_rst = Dataset(rst_file, 'r')


# In[ ]:


# read grid information
lon = np.degrees(f_rst.variables['lonCell'][:])
lat = np.degrees(f_rst.variables['latCell'][:])
cellArea = f_rst.variables['areaCell'][:]

refBottomDepth = f_rst.variables['refBottomDepth'][:]
nVertLevels = len(refBottomDepth)
refTopDepth = np.zeros(nVertLevels)
refTopDepth[1:nVertLevels] = refBottomDepth[0:nVertLevels-1]
refLayerThickness = refTopDepth-refBottomDepth
refMidDepth = 0.5*(refTopDepth+refBottomDepth)

# ### Cross section

# In[ ]:


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


# ### Mixed Layer Depth (m)

# In[ ]:


fldname = 'mixedLayerDepth'
varname = 'timeMonthly_avg_dThreshMLD'
units = 'm'
levels = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                   110, 130, 150, 180, 210, 240, 280, 320, 360,
                   407, 454, 500, 1000, 1500, 2000])
ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(fldname, varname, ts_ys, ts_ye)
f_ts = Dataset(ts_file, 'r')
data = f_ts.variables[varname][:]
fig_dir = fig_root+'/Animation/'+varname
os.makedirs(fig_dir, exist_ok=True)
for y in np.arange(plt_ye-plt_ys+1):
    yyyy = '{:04d}'.format(y+1)
    for m in np.arange(nmon):
        fig = plt.figure(figsize=[6,5.5])
        mm = '{:02d}'.format(m+1)
        i = 12*y+m
        mpaso_obj = MPASOMap(data=data[i,:], lat=lat, lon=lon, cellarea=cellArea, name=fldname, units=units)
        mpaso_obj.plot(region='LabSea', levels=levels, ptype='contourf')
        axis = plt.gca()
        axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)



# ### Salinity (psu) with ice fraction (unitless)

# In[ ]:


fldname = 'salinity'
varname = 'timeMonthly_avg_activeTracers_salinity'
units = 'psu'
levels = np.linspace(28, 36, 41)
ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(fldname, varname, ts_ys, ts_ye)
f_ts = Dataset(ts_file, 'r')
fldname1 = 'iceFraction'
varname1 = 'timeMonthly_avg_iceAreaCell'
levels1 = [0.15, 0.85]
ts_file1 = ts_root_ice+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(fldname1, varname1, ts_ys, ts_ye)
f_ts1 = Dataset(ts_file1, 'r')

fig_dir = fig_root+'/Animation/'+varname
os.makedirs(fig_dir, exist_ok=True)
for y in np.arange(plt_ye-plt_ys+1):
    yyyy = '{:04d}'.format(y+1)
    for m in np.arange(nmon):
        # figure 1: map
        fig = plt.figure(figsize=[6,5.5])
        mm = '{:02d}'.format(m+1)
        i = 12*y+m
        data = f_ts.variables[varname][i,:,:]
        mpasovol_obj = MPASOVolume(data=data, lat=lat, lon=lon, depth=refMidDepth, cellarea=cellArea, name=fldname, units=units)
        mpaso_obj = mpasovol_obj.get_map(depth=0.0)
        m = mpaso_obj.plot(region='LabSea', levels=levels, ptype='contourf')
        data1 = f_ts1.variables[varname1][i,:]
        mpascice_obj = MPASCICEMap(data=data1, lat=lat, lon=lon, cellarea=cellArea, name=fldname1, units=None)
        mpascice_obj.overlay(m, levels=levels1, colors='w')
        m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='gray')
        m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='gray')
        axis = plt.gca()
        axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_Map_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        # figure 2: vertical cross section 1
        fig = plt.figure(figsize=[6,4])
        mpaso_vcsec1 = mpasovol_obj.get_vertical_cross_section(lon0=s1_s_lon, lat0=s1_s_lat,
                                                             lon1=s1_e_lon, lat1=s1_e_lat)
        mpaso_vcsec1.plot(levels=levels, depth_mode='native')
        axis = plt.gca()
        axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_VCSec1_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        # figure 3: vertical cross section 2
        fig = plt.figure(figsize=[6,4])
        mpaso_vcsec2 = mpasovol_obj.get_vertical_cross_section(lon0=s2_s_lon, lat0=s2_s_lat,
                                                             lon1=s2_e_lon, lat1=s2_e_lat, depth_bottom=4500)
        mpaso_vcsec2.plot(levels=levels, depth_mode='native')
        axis = plt.gca()
        axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_VCSec2_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)

f_ts.close()
f_ts1.close()
del data, data1

# ### Temperature (degC) with ice thickness (m)

# In[ ]:


fldname = 'temperature'
varname = 'timeMonthly_avg_activeTracers_temperature'
units = 'degC'
levels = np.linspace(-2, 26, 57)
ts_file = ts_root+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(fldname, varname, ts_ys, ts_ye)
f_ts = Dataset(ts_file, 'r')

fldname1 = 'iceThickness'
varname1 = 'timeMonthly_avg_iceVolumeCell'
levels1 = np.linspace(0.5, 5, 10)
ts_file1 = ts_root_ice+'/{:s}/{:s}_{:04d}01_{:04d}12.nc'.format(fldname1, varname1, ts_ys, ts_ye)
f_ts1 = Dataset(ts_file1, 'r')

fig_dir = fig_root+'/Animation/'+varname
os.makedirs(fig_dir, exist_ok=True)
for y in np.arange(plt_ye-plt_ys+1):
    yyyy = '{:04d}'.format(y+1)
    for m in np.arange(nmon):
        # figure 1: map
        fig = plt.figure(figsize=[6,5.5])
        mm = '{:02d}'.format(m+1)
        i = 12*y+m
        data = f_ts.variables[varname][i,:,:]
        mpasovol_obj = MPASOVolume(data=data, lat=lat, lon=lon, depth=refMidDepth, cellarea=cellArea, name=fldname, units=units)
        mpaso_obj = mpasovol_obj.get_map(depth=0.0)
        m = mpaso_obj.plot(region='LabSea', levels=levels, ptype='contourf')
        data1 = f_ts1.variables[varname1][i,:]
        mpascice_obj = MPASCICEMap(data=data1, lat=lat, lon=lon, cellarea=cellArea, name=fldname1, units=None)
        mpascice_obj.overlay(m, levels=levels1, cmap='bone_r')
        m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='gray')
        m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='gray')
        axis = plt.gca()
        axis.text(0.06, 0.62, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_Map_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        # figure 2: vertical cross section 1
        fig = plt.figure(figsize=[6,4])
        mpaso_vcsec1 = mpasovol_obj.get_vertical_cross_section(lon0=s1_s_lon, lat0=s1_s_lat,
                                                             lon1=s1_e_lon, lat1=s1_e_lat)
        mpaso_vcsec1.plot(levels=levels, depth_mode='native')
        axis = plt.gca()
        axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_VCSec1_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        # figure 3: vertical cross section 2
        fig = plt.figure(figsize=[6,4])
        mpaso_vcsec2 = mpasovol_obj.get_vertical_cross_section(lon0=s2_s_lon, lat0=s2_s_lat,
                                                             lon1=s2_e_lon, lat1=s2_e_lat, depth_bottom=4500)
        mpaso_vcsec2.plot(levels=levels, depth_mode='native')
        axis = plt.gca()
        axis.text(0.06, 0.12, yyyy+'-'+mm, transform=axis.transAxes,
                     fontsize=12, color='k', va='top',
                     bbox=dict(boxstyle='square',ec='k',fc='w'))
        plt.show()
        figname = fig_dir+'/LabSea_climo_VCSec2_'+yyyy+'-'+mm+'.png'
        fig.savefig(figname, dpi = 300)
        plt.close(fig)

f_ts.close()
f_ts1.close()
del data, data1
