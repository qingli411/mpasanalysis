import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
from scipy import spatial
from mpl_toolkits.basemap import Basemap

#--------------------------------
# MPASOVolume
#--------------------------------

class MPASOVolume(object):

    """MPASOVolume object"""

    def __init__(self, data, lon, lat, depth, cellarea, name, units):
        """Iniitalize MPASOVolume

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :depth: (1D numpy array) depth
        :cellarea: (1D numpy array) area of cells
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.fillvalue = -1.e30
        self.data = np.where(data<=self.fillvalue, np.nan, data)
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.cellarea = cellarea
        self.name = name
        self.units = units

    def get_map(self, depth=0.0):
        """ Return a map at depth

        :depth: (float) depth of the cross section
        :returns: (MPASOMap object) map

        """
        # get z index
        zidx = np.argmin(np.abs(self.depth-depth))
        # MPASOMap object
        obj = MPASOMap(data=self.data[:,zidx], lon=self.lon, lat=self.lat, cellarea=self.cellarea, name=self.name, units=self.units)
        return obj

    def get_vertical_cross_section(self, lon0, lat0, lon1, lat1, depth_bottom=5000.0, depth_top=0.0):
        """Return cross section defined by two points [lon0, lat0] and [lon1, lat1] and
           the depth range [depth_bottom, depth_top]

        :lon0: (float) Longitude of starting point
        :lat0: (float) Latitude of starting point
        :lon1: (float) Longitude of ending point
        :lat1: (float) Latitude of ending point
        :depth_bottom: (float) depth of the bottom
        :depth_top: (float) depth of the top
        :returns: (MPASOVertCrossSection object) vertical cross section

        """
        # automatically adjust number of points
        d = distance(lon0, lat0, lon1, lat1)
        d_cell = np.sqrt(self.cellarea.mean())
        npoints = int(np.ceil(d/d_cell*1000))
        print('Nearest neighbor interpolation to {} points.'.format(npoints))
        xx = np.linspace(lon0, lon1, npoints)
        yy = np.linspace(lat0, lat1, npoints)
        # select nearest neighbor
        pts = np.array(list(zip(xx,yy)))
        tree = spatial.KDTree(list(zip(self.lon, self.lat)))
        p = tree.query(pts)
        cidx = p[1]
        zidx0 = np.argmin(np.abs(self.depth-depth_top))
        zidx1 = np.argmin(np.abs(self.depth-depth_bottom))
        data = self.data[cidx,zidx0:zidx1]
        depth = self.depth[zidx0:zidx1]
        # calculate distance from [lon0, lat0]
        dist = np.zeros(npoints)
        for i in np.arange(npoints-1):
            dist[i+1] = distance(lon0, lat0, xx[i+1], yy[i+1])
        obj = MPASOVertCrossSection(data=data, lon=xx, lat=yy, dist=dist, depth=depth, name=self.name, units=self.units)
        return obj

#--------------------------------
# MPASOMap
#--------------------------------

class MPASOMap(object):

    """MPASOMap object"""

    def __init__(self, data, lon, lat, cellarea, name, units):
        """Initialize MPASOMap

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :cellarea: (1D numpy array) area of cells
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.data = data
        self.lon = lon
        self.lat = lat
        self.cellarea = cellarea
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


    def plot(self, axis=None, region='Global', levels=None, add_title=True, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) to be passed to mpl_toolkits.basemap.scatter()
        :return: (basemap) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot map
        if region == 'Global':
            # global map
            lon_ll = 20.0
            lat_ll = -80.0
            lon_ur = 380.0
            lat_ur = 80.0
            m = Basemap(projection='cyl', llcrnrlat=lat_ll, urcrnrlat=lat_ur, llcrnrlon=lon_ll, urcrnrlon=lon_ur, ax=axis)
            # parallels and meridians
            mdlat = 30.0
            mdlon = 60.0
            # markersize
            markersize = 1
            # data
            data = self.data
            lat = self.lat
            lon = self.lon
            cellarea = self.cellarea
            # shift longitude
            lon = np.where(lon < 20., lon+360., lon)
        else:
            # regional map
            region_obj = region_latlon(region)
            lon_ll, lat_ll, lon_ur, lat_ur = region_obj.lon_ll, region_obj.lat_ll, region_obj.lon_ur, region_obj.lat_ur
            lon_c = 0.5*(lon_ll+lon_ur)
            lat_c = 0.5*(lat_ll+lat_ur)
            m = Basemap(projection='cass', llcrnrlon=lon_ll, llcrnrlat=lat_ll, urcrnrlon=lon_ur, urcrnrlat=lat_ur,
                        resolution='l', lon_0=lon_c, lat_0=lat_c, ax=axis)
            # parallels and meridians
            mdlat = 10.0
            mdlon = 10.0
            # region mask
            lon_mask = (self.lon >= lon_ll-26.0) & (self.lon <= lon_ur)
            lat_mask = (self.lat >= lat_ll) & (self.lat <= lat_ur+4.0)
            region_mask = lon_mask & lat_mask
            # apply region mask to data
            data = self.data[region_mask]
            lat = self.lat[region_mask]
            lon = self.lon[region_mask]
            cellarea = self.cellarea[region_mask]
        # print message
        print('Plotting map of {} at region \'{}\''.format(self.name+' ('+self.units+')', region))
        # plot coastlines, draw label meridians and parallels.
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.fillcontinents(color='gray',lake_color='lightgray')
        m.drawparallels(np.arange(-90.,91.,mdlat), labels=[1,0,0,1])
        m.drawmeridians(np.arange(-180.,181.,mdlon), labels=[1,0,0,1])
        # automatically adjust the marker size for regional plot
        if region != 'Global':
            plt.gcf().canvas.draw()
            axwidth = m.ax.get_window_extent().width/72.
            axheight = m.ax.get_window_extent().height/72.
            axarea = axwidth*axheight
            cellarea_mean = np.nanmean(cellarea)
            cellarea_norm = cellarea/cellarea_mean
            area = m.xmax * m.ymax
            markersize = cellarea_mean*18000/area*axarea*cellarea_norm
            markersize[markersize<1.0] = 1.0
            print('Minimum and maximum markersizes: {:4.2f} and {:4.2f}'.format(np.min(markersize), np.max(markersize)))
            if np.max(markersize) < 1.0:
                print('Set markersize to 1')
                markersize = 1
        # plot data on map
        x, y = m(lon, lat)
        if levels is not None:
            # manually mapping levels to the colormap if levels is passed in,
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            fig = m.scatter(x, y, marker='.', s=markersize, c=data, norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            # otherwise linear mapping
            fig = m.scatter(x, y, marker='.', s=markersize, c=data, cmap=plt.cm.get_cmap(cmap), **kwargs)
        # add title
        if add_title:
            axis.set_title('{} ({})'.format(self.name, self.units))
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        return m

#--------------------------------
# MPASOVertCrossSection
#--------------------------------

class MPASOVertCrossSection(object):

    """MPASOVertCrossSection object"""

    def __init__(self, data, lon, lat, dist, depth, name, units):
        """Initialize MPASOCrossSection

        :data: (1D numpy array) data array
        :lon: (1D numpy array) longitude array
        :lat: (1D numpy array) latitude array
        :dist: (1D numpy array) distance array
        :depth: (1D numpy array) depth array
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.data = data
        self.lon = lon
        self.lat = lat
        self.dist = dist
        self.depth = depth
        self.name = name
        self.units = units

    def plot(self, axis=None, ptype='contourf', levels=None, add_title=True, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) to be passed to mpl_toolkits.basemap.scatter()
        :return: (basemap) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        if levels:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        if ptype == 'pcolor':
            fig = axis.pcolor(self.dist, self.depth, np.transpose(self.data), norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        elif ptype == 'contourf':
            fig = axis.contourf(self.dist, self.depth, np.transpose(self.data), norm=norm, cmap=plt.cm.get_cmap(cmap),  **kwargs)

        axis.set_xlabel('Distance (km)')
        axis.set_ylabel('Depth (m)')
        axis.invert_yaxis()
        # add title
        if add_title:
            axis.set_title('{} ({})'.format(self.name, self.units))
        # add colorbar
        if add_colorbar:
            cb = plt.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        return fig

#--------------------------------
# Region object
#--------------------------------

class regionMap(object):
    """ Region object"""
    def __init__(self, lon_ll, lat_ll, lon_ur, lat_ur):
        self.lon_ll = lon_ll
        self.lat_ll = lat_ll
        self.lon_ur = lon_ur
        self.lat_ur = lat_ur

#--------------------------------
# Share functions
#--------------------------------

def region_latlon(region_name):
    """Return longitude and latitude of lower left an upper right of the region.

    :region_name: (string) region name
    :return: (region object) predefined region object

    """

    if region_name == 'LabSea':
        # regional map for Labrador sea
        region = regionMap(lon_ll=296.0, lat_ll=36.0, lon_ur=356.0, lat_ur=70.0)
    elif region_name == 'test':
        region = regionMap(lon_ll=310.0, lat_ll=55.0, lon_ur=320.0, lat_ur=65.0)
    else:
        raise ValueError('Region {} not supported.'.format(region_name))
    return region

def distance(lon0, lat0, lon1, lat1):
    """Calculate the distance between two points on Earth.
    https://stackoverflow.com/questions/19412462/getting-distance-between-two-points-based-on-latitude-longitude

    """
    radius = 6371  # km
    dlat = np.radians(lat1 - lat0)
    dlon = np.radians(lon1 - lon0)
    a = (np.sin(dlat / 2) * np.sin(dlat / 2) +
         np.cos(np.radians(lat0)) * np.cos(np.radians(lat1)) *
         np.sin(dlon / 2) * np.sin(dlon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = radius * c
    return d

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
