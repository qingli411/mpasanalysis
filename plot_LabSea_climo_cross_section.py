
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


# flags


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


# ## Cross sections

# In[ ]:


# cross section 1
s1_s_lat = 62
s1_s_lon = 300
s1_e_lat = 50
s1_e_lon = 320
# cross section 2
s2_s_lat = 54
s2_s_lon = 305
s2_e_lat = 60.5
s2_e_lon = 313
# list of depths
depth_list = [0, 100, 200, 500, 1000]


# ## Bathymetry

# In[ ]:


fig = plt.figure(figsize=[6,6])
levels = np.linspace(0, 6000, 13)
bottomDepth = f_rst.variables['bottomDepth'][:]
mpaso_bath = MPASOMap(data=bottomDepth, lat=lat, lon=lon, cellarea=cellArea, name='Bottom Depth', units='m')
m = mpaso_bath.plot(region='LabSea', ptype='contourf', cmap='bone_r', levels=levels)
m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='y')
m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='y')
plt.show()
figname = fig_root+'/LabSea_bathymetry.png'
fig.savefig(figname, dpi = 300)
plt.close(fig)


# ## Temperature

# In[ ]:


for i in np.arange(12):
    climo_file = climo_root+'/temperature/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    ncvar_temp = f_climo.variables['timeMonthly_avg_activeTracers_temperature']
    temp = ncvar_temp[0,:,:]
    mpaso_temp = MPASOVolume(data=temp, lon=lon, lat=lat, depth=refMidDepth, cellarea=cellArea,
                             name='Temperature', units='degC')
    levels = np.linspace(-2, 26, 57)
    
    # Temperature map at different depths
    for depth in depth_list:
        fig = plt.figure(figsize=[6,6])
        mpaso_dat = mpaso_temp.get_map(depth=depth)
        m = mpaso_dat.plot(region='LabSea', levels=levels)
        m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='y')
        m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='y')
        plt.show()
        figname = fig_root+'/LabSea_climo_{:02d}_Map_temperature_D{:d}.png'.format(i+1, depth)
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        
    # Vertical cross sections of temeprature
    # section 1
    fig = plt.figure(figsize=[6,4])
    mpaso_vcsec1 = mpaso_temp.get_vertical_cross_section(lon0=s1_s_lon, lat0=s1_s_lat,
                                                         lon1=s1_e_lon, lat1=s1_e_lat)
    mpaso_vcsec1.plot(levels=levels, depth_mode='native')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_VCSec1_temperature.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # section 2
    fig = plt.figure(figsize=[6,4])
    mpaso_vcsec2 = mpaso_temp.get_vertical_cross_section(lon0=s2_s_lon, lat0=s2_s_lat,
                                                         lon1=s2_e_lon, lat1=s2_e_lat, depth_bottom=4500)
    mpaso_vcsec2.plot(levels=levels, depth_mode='native')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_VCSec2_temperature.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)


# ## Salinity

# In[ ]:


for i in np.arange(12):
    climo_file = climo_root+'/salinity/mpaso_{:02d}_climo.nc'.format(i+1)
    f_climo = Dataset(climo_file, 'r')
    ncvar_salt = f_climo.variables['timeMonthly_avg_activeTracers_salinity']
    salt = ncvar_salt[0,:,:]
    mpaso_salt = MPASOVolume(data=salt, lon=lon, lat=lat, depth=refMidDepth, cellarea=cellArea,
                             name='Salinity', units='psu')
    levels = np.linspace(16, 40, 49)
    
    # Salinity map at different depths
    for depth in depth_list:
        fig = plt.figure(figsize=[6,6])
        mpaso_dat = mpaso_salt.get_map(depth=depth)
        m = mpaso_dat.plot(region='LabSea', levels=levels)
        m.drawgreatcircle(s1_s_lon, s1_s_lat, s1_e_lon, s1_e_lat, color='y')
        m.drawgreatcircle(s2_s_lon, s2_s_lat, s2_e_lon, s2_e_lat, color='y')
        plt.show()
        figname = fig_root+'/LabSea_climo_{:02d}_Map_salinity_D{:d}.png'.format(i+1, depth)
        fig.savefig(figname, dpi = 300)
        plt.close(fig)
        
    # Vertical cross sections of temeprature
    # section 1
    fig = plt.figure(figsize=[6,4])
    mpaso_vcsec1 = mpaso_salt.get_vertical_cross_section(lon0=s1_s_lon, lat0=s1_s_lat,
                                                         lon1=s1_e_lon, lat1=s1_e_lat)
    mpaso_vcsec1.plot(levels=levels, depth_mode='native')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_VCSec1_salinity.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)
    # section 2
    fig = plt.figure(figsize=[6,4])
    mpaso_vcsec2 = mpaso_salt.get_vertical_cross_section(lon0=s2_s_lon, lat0=s2_s_lat,
                                                         lon1=s2_e_lon, lat1=s2_e_lat, depth_bottom=4500)
    mpaso_vcsec2.plot(levels=levels, depth_mode='native')
    plt.show()
    figname = fig_root+'/LabSea_climo_{:02d}_VCSec2_salinity.png'.format(i+1)
    fig.savefig(figname, dpi = 300)
    plt.close(fig)

