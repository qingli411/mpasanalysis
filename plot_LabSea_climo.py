
# coding: utf-8

# In[ ]:


from mpasanalysis import *
import e3sm_res_cmp
get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


# get paths of restart files, monthly mean output files, processed climatology files and output figures
data_root = e3sm_res_cmp.load_paths_ocn(climo_ys=41, climo_ye=50, ts_ys=1, ts_ye=50)
rst_root = data_root['rst_root']
climo_root = data_root['climo_root']
fig_root = data_root['fig_root']
rst_file = rst_root+'/mpaso.rst.0051-01-01_00000.nc'


# In[ ]:


# load dataset
f_rst = Dataset(rst_file, 'r')


# In[ ]:


# read grid information
lon = np.degrees(f_rst.variables['lonCell'][:])
lat = np.degrees(f_rst.variables['latCell'][:])
cellArea = f_rst.variables['areaCell'][:]


# ### SST (degC)

# In[ ]:


levels = np.linspace(-2, 26, 57)
for i in np.arange(12):
    climo_file = climo_root+'/temperature/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_activeTracers_temperature'][0,:,0]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea, name='SST', units='degC')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1])
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_SST.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### SSS (psu)

# In[ ]:


levels = np.linspace(16, 36, 41)
for i in np.arange(12):
    climo_file = climo_root+'/salinity/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_activeTracers_salinity'][0,:,0]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea, name='SSS', units='psu')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1])
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_SSS.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Sea Surface Potential Density (kg m^{-2}

# In[ ]:


levels = np.linspace(1012, 1028, 33)
for i in np.arange(12):
    climo_file = climo_root+'/potentialDensity/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_potentialDensity'][0,:,0]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                         name='Potential Density', units='kg m^{-3}')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1])
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_potentialDensity.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Sea Surface Current (m s^{-1})

# In[ ]:


levels = np.linspace(-0.4, 0.4, 21)
for i in np.arange(12):
    climo_file = climo_root+'/velocity/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    # zonal
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_velocityZonal'][0,:,0]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                         name='Zonal Surface Current', units='m s^{-1}')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_velocityZonal.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # meridional
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_velocityMeridional'][0,:,0]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                         name='Meridional Surface Current', units='m s^{-1}')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_velocityMeridional.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Mixed Layer Depth (m)

# In[ ]:


levels = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                   110, 130, 150, 180, 210, 240, 280, 320, 360,
                   407, 454, 500, 1000, 1500, 2000])
levels_diff = np.linspace(-50, 50, 21)
for i in np.arange(12):
    # read data
    climo_file = climo_root+'/mixedLayerDepth/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    # mld_d
    fig = plt.figure(figsize=[6,6])
    mld_d = f_climo.variables['timeMonthly_avg_dThreshMLD'][0,:]
    mpaso_mld_d = MPASOMap(data=mld_d, lat=lat, lon=lon, cellarea=cellArea,
                           name='MLD (density threshold)', units='m')
    mpaso_mld_d.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1])
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_dThreshMLD.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # mld_t
    fig = plt.figure(figsize=[6,6])
    mld_t = f_climo.variables['timeMonthly_avg_tThreshMLD'][0,:]
    mpaso_mld_t = MPASOMap(data=mld_t, lat=lat, lon=lon, cellarea=cellArea,
                         name='MLD (temperature threshold)', units='m')
    mpaso_mld_t.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1])
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_tThreshMLD.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # diff_mld
    fig = plt.figure(figsize=[6,6])
    mpaso_diff = MPASOMap(data=mld_t-mld_d, lat=lat, lon=lon, cellarea=cellArea,
                         name='MLD (temperature) - MLD (density)', units='m')
    mpaso_diff.plot(region='LabSea', levels=levels_diff, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_tThreshMLD-dThreshMLD.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Surface Heat Flux (W m$^{-2}$)

# In[ ]:


levels = np.linspace(-400, 400, 81)
for i in np.arange(12):
    # read file
    climo_file = climo_root+'/heatFlux/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    # shf
    fig = plt.figure(figsize=[6,6])
    shf = f_climo.variables['timeMonthly_avg_sensibleHeatFlux'][0,:]
    mpaso_shf = MPASOMap(data=shf, lat=lat, lon=lon, cellarea=cellArea, name='SHF', units=r'W m^{-2}')
    mpaso_shf.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_sensibleHeatFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # lhf
    fig = plt.figure(figsize=[6,6])
    lhf = f_climo.variables['timeMonthly_avg_latentHeatFlux'][0,:]
    mpaso_lhf = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea, name='LHF', units=r'W m^{-2}')
    mpaso_lhf.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_latentHeatFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # swf
    fig = plt.figure(figsize=[6,6])
    swf = f_climo.variables['timeMonthly_avg_shortWaveHeatFlux'][0,:]
    mpaso_shf = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea, name='SWF', units=r'W m^{-2}')
    mpaso_shf.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_shortWaveHeatFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # lwf
    fig = plt.figure(figsize=[6,6])
    lwfu = f_climo.variables['timeMonthly_avg_longWaveHeatFluxUp'][0,:]
    lwfd = f_climo.variables['timeMonthly_avg_longWaveHeatFluxDown'][0,:]
    lwf = lwfu+lwfd
    mpaso_lwf = MPASOMap(data=lwf, lat=lat, lon=lon, cellarea=cellArea, name='LWF', units=r'W m^{-2}')
    mpaso_lwf.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_longWaveHeatFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # total
    fig = plt.figure(figsize=[6,6])
    tot = shf+lhf+swf+lwf
    mpaso_tot = MPASOMap(data=tot, lat=lat, lon=lon, cellarea=cellArea, name='Total Heat Flux', units=r'W m^{-2}')
    mpaso_tot.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_totalHeatFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Surface Freshwater  Flux (kg m$^{-2}$ s$^{-1}$)

# In[ ]:


levels = np.linspace(-2e-4, 2e-4, 41)
for i in np.arange(12):
    # read file
    climo_file = climo_root+'/freshWaterFlux/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    # evap
    fig = plt.figure(figsize=[6,6])
    evap = f_climo.variables['timeMonthly_avg_evaporationFlux'][0,:]
    mpaso_evap = MPASOMap(data=evap, lat=lat, lon=lon, cellarea=cellArea,
                          name='Evaporation', units=r'kg m^{-2} s^{-1}')
    mpaso_evap.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_evaporationFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # prec
    fig = plt.figure(figsize=[6,6])
    rain = f_climo.variables['timeMonthly_avg_rainFlux'][0,:]
    snow = f_climo.variables['timeMonthly_avg_snowFlux'][0,:]
    prec = rain+snow
    mpaso_prec = MPASOMap(data=prec, lat=lat, lon=lon, cellarea=cellArea,
                          name='Precipitation', units=r'kg m^{-2} s^{-1}')
    mpaso_prec.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_precipitationFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # sea ice
    fig = plt.figure(figsize=[6,6])
    si_s = f_climo.variables['timeMonthly_avg_seaIceSalinityFlux'][0,:]
    si_f = f_climo.variables['timeMonthly_avg_seaIceFreshWaterFlux'][0,:]
    sice = si_s+si_f
    mpaso_sice = MPASOMap(data=sice, lat=lat, lon=lon, cellarea=cellArea,
                          name='Sea Ice', units=r'kg m^{-2} s^{-1}')
    mpaso_sice.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_seaIceFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # run off
    fig = plt.figure(figsize=[6,6])
    ro_r = f_climo.variables['timeMonthly_avg_riverRunoffFlux'][0,:]
    ro_i = f_climo.variables['timeMonthly_avg_iceRunoffFlux'][0,:]
    roff = ro_r+ro_i
    mpaso_roff = MPASOMap(data=roff, lat=lat, lon=lon, cellarea=cellArea,
                          name='Runoff', units=r'kg m^{-2} s^{-1}')
    mpaso_roff.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_runOffFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # total
    fig = plt.figure(figsize=[6,6])
    tot = evap+prec+sice+roff
    mpaso_tot = MPASOMap(data=tot, lat=lat, lon=lon, cellarea=cellArea,
                         name='Total Fresh Water Flux', units=r'kg m^{-2} s^{-1}')
    mpaso_tot.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_totalFreshWaterFlux.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ### Wind stress (N m^{-2})

# In[ ]:


levels = np.linspace(-0.4, 0.4, 41)
for i in np.arange(12):
    climo_file = climo_root+'/windStress/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    # zonal
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_windStressZonal'][0,:]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                         name='Zonal Wind Stress', units='N m^{-2}')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_windStressZonal.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # meridional
    fig = plt.figure(figsize=[6,6])
    data = f_climo.variables['timeMonthly_avg_windStressMeridional'][0,:]
    mpaso_obj = MPASOMap(data=data, lat=lat, lon=lon, cellarea=cellArea,
                         name='Meridional Wind Stress', units='N m^{-2}')
    mpaso_obj.plot(region='LabSea', levels=levels, label=calendar.month_abbr[i+1], cmap='RdBu_r')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_windStressMeridional.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

