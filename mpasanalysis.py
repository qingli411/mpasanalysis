import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

#--------------------------------
# GOMTMap
#--------------------------------

class MPASOMap(object):

    """MPASOMap object"""

    def __init__(self, data=None, lon=None, lat=None, name=None, units=None):
        """Initialize MPASOMap

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.data = data
        self.lon = lon
        self.lat = lat
        self.name = name
        self.units = units

    def save(self, path):
        """Save MPASOMap object

        :path: (str) path of file to save
        :returns: none

        """
        np.savez(path, data=self.data, lon=self.lon, lat=self.lat, name=self.name, units=self.units)

    def load(self, path):
        """Load data to MPASOMap object

        :path: (str) path of file to load
        :returns: (MPASOMap object)

        """
        dat = np.load(path)
        self.__init__(data=dat['data'], lon=dat['lon'], lat=dat['lat'],
                name=str(dat['name']), units=str(dat['units']))
        return self

    def masked(self, mask, mask_data=np.nan):
        """Apply mask to MPASOMap object. The mask should also be a MPASOMap object,
           with 1 for valid and 0 for invalid.

        :mask: (MPASOMap object) mask, 1 for valid, 0 for invalid
        :mask_data: (optional) values to be filled in maked points
        :return: (MPASOMap object) masked MPASOMap

        """
        if mask.data.size != self.data.size:
            raise ValueError('The dimension of mask does not match.')
        dat = self.data
        self.data = np.where(mask.data==0, mask_data, dat)


    def plot(self, axis=None, region='Global', levels=None, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :leveles: (list, optional) list of levels
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) to be passed to mpl_toolkits.basemap.scatter()
        :return: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot map
        if region == 'Global':
            # global map
            m = Basemap(projection='cyl', llcrnrlat=-72, urcrnrlat=72, llcrnrlon=20, urcrnrlon=380, ax=axis)
            # plot coastlines, draw label meridians and parallels.
            m.drawcoastlines()
            m.drawmapboundary(fill_color='lightgray')
            m.fillcontinents(color='gray',lake_color='lightgray')
            m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,0,1])
            data = self.data
            lat = self.lat
            lon = self.lon
            # shift longitude
            lon = np.where(lon < 20., lon+360., lon)
            x, y = m(lon, lat)
            # marker size
            markersize = 0.1
        elif region == 'LabSea':
            # regional map for Labrador sea
            lon_ll = 296.0
            lat_ll = 36.0
            lon_ur = 356.0
            lat_ur = 70.0
            lon_c = 0.5*(lon_ll+lon_ur)
            lat_c = 0.5*(lat_ll+lat_ur)
            m = Basemap(projection='cass', llcrnrlon=lon_ll, llcrnrlat=lat_ll, urcrnrlon=lon_ur, urcrnrlat=lat_ur,
                        resolution='l', lon_0=lon_c, lat_0=lat_c, ax=axis)
            m.drawcoastlines()
            m.drawmapboundary(fill_color='lightgray')
            m.fillcontinents(color='gray',lake_color='lightgray')
            m.drawparallels(np.arange(-90.,91.,10.), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-180.,181.,10.), labels=[1,0,0,1])
            lon_mask = (self.lon >= lon_ll-26.0) & (self.lon <= lon_ur)
            lat_mask = (self.lat >= lat_ll) & (self.lat <= lat_ur+4.0)
            region_mask = lon_mask & lat_mask
            data = self.data[region_mask]
            lat = self.lat[region_mask]
            lon = self.lon[region_mask]
            x, y = m(lon, lat)
            # marker size
            markersize = 1
        elif region == 'test':
            # regional map for test
            lon_ll = 310.0
            lat_ll = 55.0
            lon_ur = 320.0
            lat_ur = 65.0
            lon_c = 0.5*(lon_ll+lon_ur)
            lat_c = 0.5*(lat_ll+lat_ur)
            m = Basemap(projection='cass', llcrnrlon=lon_ll, llcrnrlat=lat_ll, urcrnrlon=lon_ur, urcrnrlat=lat_ur,
                        resolution='l', lon_0=lon_c, lat_0=lat_c, ax=axis)
            m.drawcoastlines()
            m.drawmapboundary(fill_color='lightgray')
            m.fillcontinents(color='gray',lake_color='lightgray')
            m.drawparallels(np.arange(-90.,91.,10.), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-180.,181.,10.), labels=[1,0,0,1])
            lon_mask = (self.lon >= lon_ll-26.0) & (self.lon <= lon_ur)
            lat_mask = (self.lat >= lat_ll) & (self.lat <= lat_ur+4.0)
            region_mask = lon_mask & lat_mask
            data = self.data[region_mask]
            lat = self.lat[region_mask]
            lon = self.lon[region_mask]
            x, y = m(lon, lat)
            # marker size
            markersize = 2
        else:
            raise ValueError('Region {} not supported.'.format(region))
        # manually mapping levels to the colormap if levels is passed in,
        # otherwise linear mapping
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            fig = m.scatter(x, y, marker='.', s=markersize, c=data, norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            fig = m.scatter(x, y, marker='.', s=markersize, c=data, cmap=plt.cm.get_cmap(cmap), **kwargs)
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-2, 2))
            cb.update_ticks()
        return fig

def plot_map(xx, yy, data, axis=None, levels=None, add_colorbar=True, draw_eq=False, cmap='rainbow', **kwargs):
    # use curret axis if not specified
    if not axis:
        axis = plt.gca()
    # plot map
    m = Basemap(projection='cyl', llcrnrlat=-85, urcrnrlat=85, llcrnrlon=20, urcrnrlon=380, ax=axis)
    # plot coastlines, draw label meridians and parallels.
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='gray',lake_color='lightgray')
    m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1])
    m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,0,1])
    # manually mapping levels to the colormap if levels is passed in,
    # otherwise linear mapping
    if levels:
        bounds = np.array(levels)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        fig = m.pcolormesh(x=xx, y=yy, data=data, latlon=True, norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
    else:
        fig = m.pcolormesh(x=xx, y=yy, data=data, latlon=True, cmap=plt.cm.get_cmap(cmap), **kwargs)
    # draw equator if requested
    if draw_eq:
        m.drawparallels([0], labels=[0,0,0,0], linewidth=1.5, dashes=[1, 0])
    # add colorbar
    if add_colorbar:
        cb = m.colorbar(fig, ax=axis)
        # cb.formatter.set_powerlimits((-2, 2))
        cb.update_ticks()
    return fig

def rmse_lat(data0, data1, lat, lat_bnd, wgt=None):
    # weight
    if wgt is None:
        wgt = np.ones(data0.shape)
    # number of lat bin = number of lat boundary + 1
    nl = len(lat_bnd)
    rmse = np.zeros(nl+2)
    # overall rmse
    rmse[0] = np.sqrt(((data1-data0)**2*wgt).mean()/wgt.mean())
    # find rmse in each lat bin
    lidx0 = 0
    for i in np.arange(nl):
        j = i+1
        # find latitude index
        lidx1 = np.argmin(np.abs(lat-lat_bnd[i]))+1
        # print('lat0 = {}, lat1 = {}'.format(lat[lidx0], lat[lidx1]))
        rmse[j] = np.sqrt(((data1[lidx0:lidx1,:]-data0[lidx0:lidx1,:])**2*wgt[lidx0:lidx1,:]).mean()
                /wgt[lidx0:lidx1,:].mean())
        lidx0 = lidx1
    # rmse in the last lat bin
    rmse[-1] = np.sqrt(((data1[lidx0:,:]-data0[lidx0:,:])**2*wgt[lidx0:,:]).mean()/wgt[lidx0:,:].mean())
    return rmse

def me_lat(data0, data1, lat, lat_bnd, wgt=None):
    # weight
    if wgt is None:
        wgt = np.ones(data0.shape)
    # number of lat bin = number of lat boundary + 1
    nl = len(lat_bnd)
    me = np.zeros(nl+2)
    # overall rmse
    me[0] = ((data1-data0)*wgt).mean()/wgt.mean()
    # find rmse in each lat bin
    lidx0 = 0
    for i in np.arange(nl):
        j = i+1
        # find latitude index
        lidx1 = np.argmin(np.abs(lat-lat_bnd[i]))+1
        # print('lat0 = {}, lat1 = {}'.format(lat[lidx0], lat[lidx1]))
        me[j] = ((data1[lidx0:lidx1,:]-data0[lidx0:lidx1,:])*wgt[lidx0:lidx1,:]).mean()/wgt[lidx0:lidx1,:].mean()
        lidx0 = lidx1
    # rmse in the last lat bin
    me[-1] = ((data1[lidx0:,:]-data0[lidx0:,:])*wgt[lidx0:,:]).mean()/wgt[lidx0:,:].mean()
    return me
