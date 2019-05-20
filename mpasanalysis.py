import os
import sys
import calendar
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

    def __init__(self, data, lon, lat, depth, cellarea, layerthickness, name, units):
        """Iniitalize MPASOVolume

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :depth: (1D numpy array) depth
        :cellarea: (1D numpy array) area of cells
        :layerthickness: (1D numpy array) layer thickness
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.fillvalue = -9.99999979021476795361e+33
        self.data = np.where(data<=self.fillvalue, np.nan, data)
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.cellarea = cellarea
        self.layerthickness = layerthickness
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
        name = self.name+' at {:6.2f} m'.format(self.depth[zidx])
        obj = MPASOMap(data=self.data[:,zidx], lon=self.lon, lat=self.lat, cellarea=self.cellarea, name=name, units=self.units)
        return obj

    def get_map_vertical_sum(self, depth_bottom=6000.0, depth_top=0.0):
        """Return a map of vertically integrated field
        :depth_bottom: (float) depth of the bottom
        :depth_top: (float) depth of the top
        :returns: (MPASOMap object) map of vertically integrated field

        """
        # get z indices
        zidx0 = np.argmin(np.abs(self.depth-depth_top))
        zidx1 = np.argmin(np.abs(self.depth-depth_bottom))
        z_top = self.depth[zidx0]-0.5*self.layerthickness[zidx0]
        z_bottom = self.depth[zidx1]+0.5*self.layerthickness[zidx1]
        data = np.nansum(self.data[:,zidx0:zidx1+1] * self.layerthickness[zidx0:zidx1+1].reshape((1,zidx1+1-zidx0)), axis=1)
        # MPASOMap object
        name = 'Vertical sum of '+self.name+' between {:6.2f} m and {:6.2f} m'.format(z_top, z_bottom)
        units = self.units+' m'
        obj = MPASOMap(data=data, lon=self.lon, lat=self.lat, cellarea=self.cellarea, name=name, units=units)
        return obj

    def get_map_vertical_mean(self, depth_bottom=6000.0, depth_top=0.0):
        """Return a map of vertically averaged field
        :depth_bottom: (float) depth of the bottom
        :depth_top: (float) depth of the top
        :returns: (MPASOMap object) map of vertically integrated field

        """
        # get z indices
        zidx0 = np.argmin(np.abs(self.depth-depth_top))
        zidx1 = np.argmin(np.abs(self.depth-depth_bottom))
        print(zidx0)
        print(zidx1)
        z_top = self.depth[zidx0]-0.5*self.layerthickness[zidx0]
        z_bottom = self.depth[zidx1]+0.5*self.layerthickness[zidx1]
        data = np.nansum(self.data[:,zidx0:zidx1+1] * self.layerthickness[zidx0:zidx1+1].reshape((1,zidx1+1-zidx0)), axis=1)/(z_bottom-z_top)
        # MPASOMap object
        name = 'Vertical mean of '+self.name+' between {:6.2f} m and {:6.2f} m'.format(z_top, z_bottom)
        units = self.units
        obj = MPASOMap(data=data, lon=self.lon, lat=self.lat, cellarea=self.cellarea, name=name, units=units)
        return obj

    def get_vertical_cross_section(self, lon0, lat0, lon1, lat1, depth_bottom=6000.0, depth_top=0.0):
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
        d = gc_distance(lon0, lat0, lon1, lat1)
        d_cell = np.sqrt(self.cellarea.mean())
        npoints = int(np.ceil(d/d_cell*1000))
        print('Nearest neighbor interpolation to {} points.'.format(npoints))
        loni, lati = gc_interpolate(lon0, lat0, lon1, lat1, npoints)
        # adjust lon in range [0, 360)
        loni = np.where(loni<0.0, loni+360.0, loni)
        # select nearest neighbor
        pts = np.array(list(zip(loni,lati)))
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
            dist[i+1] = gc_distance(lon0, lat0, loni[i+1], lati[i+1])
        obj = MPASOVertCrossSection(data=data, lon=loni, lat=lati, dist=dist, depth=depth, name=self.name, units=self.units)
        return obj


#--------------------------------
# MPASOMap
#--------------------------------

class MPASOMap(object):

    """MPASOMap object"""

    def __init__(self, data=None, lon=None, lat=None, cellarea=None, name=None, units=None):
        """Initialize MPASOMap

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :cellarea: (1D numpy array) area of cells
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.fillvalue = -9.99999979021476795361e+33
        self.lon = lon
        self.lat = lat
        self.cellarea = cellarea
        self.name = name
        self.units = units
        if data is not None:
            self.data = np.where(data<=self.fillvalue, np.nan, data)

    def save(self, path):
        """Save MPASOMap object

        :path: (str) path of file to save
        :returns: none

        """
        np.savez(path, data=self.data, lon=self.lon, lat=self.lat, cellarea=self.cellarea,
                 name=self.name, units=self.units)

    def load(self, path):
        """Load data to MPASOMap object

        :path: (str) path of file to load
        :returns: (MPASOMap object)

        """
        dat = np.load(path)
        self.__init__(data=dat['data'], lon=dat['lon'], lat=dat['lat'],
                cellarea=dat['cellarea'], name=str(dat['name']), units=str(dat['units']))
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

    def mean(self, region='Global', lon_min=None, lat_min=None, lon_max=None, lat_max=None):
        """Area average over a region

        :region: (str) region name
        :returns: (float) mean value

        """
        data_mask = ~np.isnan(self.data)
        data = self.data[data_mask]
        lon = self.lon[data_mask]
        lat = self.lat[data_mask]
        cellarea = self.cellarea[data_mask]
        if region == 'Global':
            mean = np.sum(data*cellarea)/np.sum(cellarea)
        else:
            if region == 'Custom':
                assert all(x is not None for x in [lon_min, lat_min, lon_max, lat_max]),\
                    "Provide lon_min, lat_min, lon_max, and lat_max to customize region"
                lon_ll, lat_ll, lon_ur, lat_ur = lon_min, lat_min, lon_max, lat_max
            else:
                # region mask
                region_obj = region_latlon(region)
                lon_ll, lat_ll, lon_ur, lat_ur = region_obj.lon_ll, region_obj.lat_ll, region_obj.lon_ur, region_obj.lat_ur
            lon_mask = (lon >= lon_ll) & (lon <= lon_ur)
            lat_mask = (lat >= lat_ll) & (lat <= lat_ur)
            region_mask = lon_mask & lat_mask
            # apply region mask to data
            data = data[region_mask]
            lat = lat[region_mask]
            lon = lon[region_mask]
            cellarea = cellarea[region_mask]
            mean = np.sum(data*cellarea)/np.sum(cellarea)
        return mean

    def plot(self, axis=None, region='Global', ptype='scatter', levels=None,
             label=None, add_title=True, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :region: (str) region name
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :label: (str) label
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
            m = Basemap(projection='cyl', llcrnrlat=lat_ll, urcrnrlat=lat_ur,
                    llcrnrlon=lon_ll, urcrnrlon=lon_ur, ax=axis)
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
            m = Basemap(projection='cass', llcrnrlon=lon_ll, llcrnrlat=lat_ll,
                    urcrnrlon=lon_ur, urcrnrlat=lat_ur, resolution='l', lon_0=lon_c, lat_0=lat_c, ax=axis)
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
        m.drawcoastlines(zorder=3)
        m.drawmapboundary(fill_color='lightgray')
        m.fillcontinents(color='gray',lake_color='lightgray', zorder=2)
        m.drawparallels(np.arange(-90.,91.,mdlat), labels=[1,0,0,1])
        m.drawmeridians(np.arange(-180.,181.,mdlon), labels=[1,0,0,1])
        if levels is not None:
            # manually mapping levels to the colormap if levels is passed in,
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # plot type
        if ptype == 'scatter':
            # marker size
            if region != 'Global':
                # automatically adjust the marker size for regional plot
                plt.gcf().canvas.draw()
                axwidth = m.ax.get_window_extent().width/72.
                axheight = m.ax.get_window_extent().height/72.
                axarea = axwidth*axheight
                cellarea_mean = np.nanmean(cellarea)
                cellarea_norm = cellarea/cellarea_mean
                area = m.xmax * m.ymax
                markersize = cellarea_mean*18000/area*axarea*cellarea_norm
                markersize[markersize<1.0] = 1.0
                print('Minimum and maximum markersizes: {:4.2f} and {:4.2f}'.
                        format(np.min(markersize), np.max(markersize)))
                if np.max(markersize) < 1.0:
                    print('Set markersize to 1')
                    markersize = 1
            # plot scatter on map
            x, y = m(lon, lat)
            fig = m.scatter(x, y, marker='.', s=markersize, c=data,
                        norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        elif ptype == 'contourf':
            x, y = m(lon, lat)
            fig = m.contourf(x, y, data, tri=True, levels=levels, extend='both',
                        norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            raise ValueError('Plot type {} not supported.'.format(ptype))
        # add label
        if label is not None:
            if region == 'LabSea':
                axis.text(0.1, 0.67, label, transform=axis.transAxes,
                    fontsize=12, color='k', fontweight='bold', va='top',
                    bbox=dict(boxstyle='square',ec='k',fc='w'))
        # add title
        if add_title:
            axis.set_title('{} ({})'.format(self.name, self.units))
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        return m, fig

    def overlay(self, m, axis=None, **kwargs):
        """Overlay contours on a map.

        :m: (Basemap) map
        :returns: none

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        data = self.data
        lat = self.lat
        lon = self.lon
        x, y = m(lon, lat)
        fig = axis.tricontour(x, y, data, **kwargs)

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

    def plot(self, axis=None, ptype='contourf', depth_mode='linear', levels=None, add_title=True, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, contourf, pcolor etc.
        :depth_mode: (str) 'linear', 'native' (native grid) or 'symlog' (linear above 100 m and log below)
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
        # levels
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # depth mode
        if depth_mode == 'native':
            depth = np.arange(self.depth.size)
        elif depth_mode == 'symlog':
            depth = self.depth
        elif depth_mode == 'linear':
            depth = self.depth
        else:
            print('Depth mode \'{}\' not supported, using \'linear\' instead.'.format(depth_mode))
            depth = self.depth
        # plot type
        if ptype == 'pcolor':
            fig = axis.pcolor(self.dist, depth, np.transpose(self.data),
                    norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        elif ptype == 'contourf':
            fig = axis.contourf(self.dist, depth, np.transpose(self.data), levels=levels, extend='both',
                    norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            raise ValueError('Plot type {} not supported.'.format(ptype))
        # update depth scale
        if depth_mode == 'native':
            ndepth = self.depth.size
            plt.gcf().canvas.draw()
            ticks = axis.get_yticks()
            idx = [int(item) for item in ticks]
            depth_label = [str(int(round(self.depth[item]))) for item in idx if item < ndepth]
            axis.set_yticklabels(depth_label)
        elif depth_mode == 'symlog':
            axis.set_yscale('symlog', linthreshy=100)
            axis.axhline(y=100, linewidth=0.5, color='k')
            ylims = axis.get_ylim()
            axis.set_ylim([0, ylims[1]])
            ticks = [20, 40, 60, 80, 100, 500, 1000, 2000]
            ticks_new = [item for item in ticks if item < self.depth[-1]]
            axis.set_yticks(ticks_new)
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
# MPASCICEMap
#--------------------------------

class MPASCICEMap(MPASOMap):

    """MPASCICEMap object"""
    pass

#--------------------------------
# Region object
#--------------------------------

class region(object):
    """ Region object"""
    def __init__(self, lon_ll, lat_ll, lon_ur, lat_ur):
        self.lon_ll = lon_ll
        self.lat_ll = lat_ll
        self.lon_ur = lon_ur
        self.lat_ur = lat_ur

#--------------------------------
# Functions
#--------------------------------

def region_latlon(region_name):
    """Return longitude and latitude of lower left an upper right of the region.

    :region_name: (string) region name
    :return: (region object) predefined region object

    """

    if region_name == 'LabSea':
        # regional map for Labrador sea
        rg = region(lon_ll=296.0, lat_ll=36.0, lon_ur=356.0, lat_ur=70.0)
    elif region_name == 'LabSea_SD1':
        rg = region(lon_ll=304.0, lat_ll=56.0, lon_ur=312.0, lat_ur=60.0)
    elif region_name == 'test':
        rg = region(lon_ll=310.0, lat_ll=55.0, lon_ur=320.0, lat_ur=65.0)
    else:
        raise ValueError('Region {} not supported.'.format(region_name))
    return rg

def select_path(lonP0, latP0, lonP1, latP1,
                lonVertex, latVertex, lonEdge, latEdge,
                indexToEdgeID, indexToVertexID,
                edgesOnVertex, verticesOnEdge,
                debug_info=False):
    """ Select the edges and vertices on a path given by the two endpoints.
    """
    idxP0 = get_index_latlon([lonP0], [latP0], lonVertex, latVertex)
    idxP1 = get_index_latlon([lonP1], [latP1], lonVertex, latVertex)
    print('Vertex closest to P0: {:4.1f} {:4.1f}'.format(lonVertex[idxP0], latVertex[idxP0]))
    print('Vertex closest to P1: {:4.1f} {:4.1f}'.format(lonVertex[idxP1], latVertex[idxP1]))
    # initialize arrays
    edges_on_path        = []
    idx_edges_on_path    = []
    vertices_on_path     = []
    idx_vertices_on_path = []
    # start from vertex P0
    idx_vertex_now = idxP0
    # record vortices on path and the indices
    vertices_on_path.append(indexToVertexID[idx_vertex_now])
    idx_vertices_on_path.append(idx_vertex_now)
    if debug_info:
        print('Vertex on path: {:4.1f} {:4.1f}'.format(lonVertex[idx_vertex_now], latVertex[idx_vertex_now]))

    # continue if not reached P1
    while idx_vertex_now != idxP1:

        # find the indices of the three edges on vertex
        edge_arr     = edgesOnVertex[idx_vertex_now,:]
        idx_edge_arr = [np.where(indexToEdgeID==val)[0][0] for val in edge_arr]
        # print the location of the three edges
        if debug_info:
            for i in np.arange(len(idx_edge_arr)):
                print('   Edge {:d}: {:4.1f} {:4.1f}'.\
                      format(i, lonEdge[idx_edge_arr[i]], latEdge[idx_edge_arr[i]]))
        # choose the edge from the three that is closest to vertex P1
        dist = [gc_distance(loni, lati, lonP1, latP1) \
                for (loni, lati) in zip(lonEdge[idx_edge_arr], latEdge[idx_edge_arr])]
        idx3_next     = np.argmin(dist)
        edge_next     = edge_arr[idx3_next]
        idx_edge_next = np.where(indexToEdgeID==edge_next)[0][0]
        # print the edge on path
        if debug_info:
            print('Edge on path: [Edge {:d}] {:4.1f} {:4.1f}'.\
                  format(idx3_next, lonEdge[idx_edge_arr[idx3_next]], latEdge[idx_edge_arr[idx3_next]]))
        # record edges on path and the indices
        edges_on_path.append(edge_next)
        idx_edges_on_path.append(idx_edge_next)

        # find the other vertex on this edge
        vertex_arr      = verticesOnEdge[idx_edge_next,:]
        vertex_next     = vertex_arr[vertex_arr!=indexToVertexID[idx_vertex_now]][0]
        idx_vertex_next = np.where(indexToVertexID==vertex_next)[0][0]
        # record vortices on path and the indices
        vertices_on_path.append(vertex_next)
        idx_vertices_on_path.append(idx_vertex_next)
        if debug_info:
            print('Vertex on path: {:4.1f} {:4.1f}'.\
                  format(lonVertex[idx_vertex_next], latVertex[idx_vertex_next]))
        # move to next vertex
        idx_vertex_now  = idx_vertex_next

    out = {'edge': edges_on_path,
           'edge_idx': idx_edges_on_path,
           'vertex': vertices_on_path,
           'vertex_idx': idx_vertices_on_path}
    return out

def get_index_latlon(loni, lati, lon_arr, lat_arr):
    pts = np.array(list(zip(loni,lati)))
    tree = spatial.KDTree(list(zip(lon_arr, lat_arr)))
    p = tree.query(pts)
    cidx = p[1]
    return cidx[0]

def gc_radius():
    """Return the radius of Earth
    :returns: (float) radius of Earth in km

    """
    return 6371.0

def gc_distance(lon0, lat0, lon1, lat1):
    """Calculate the great circle distance (km) between two points [lon0, lat0] and [lon1, lat1]
    http://www.movable-type.co.uk/scripts/latlong.html

    :lon0: (float) longitude of point 1 in degrees
    :lat0: (float) longitude of point 1 in degrees
    :lon1: (float) longitude of point 2 in degrees
    :lat1: (float) longitude of point 2 in degrees
    :returns: (numpy array) longitude and latitude

    """
    radius = gc_radius() # km
    dlat_r = np.radians(lat1 - lat0)
    dlon_r = np.radians(lon1 - lon0)
    lat0_r = np.radians(lat0)
    lat1_r = np.radians(lat1)
    a = (np.sin(dlat_r / 2) * np.sin(dlat_r / 2) +
         np.cos(lat0_r) * np.cos(lat1_r) *
         np.sin(dlon_r / 2) * np.sin(dlon_r / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = radius * c
    return d

def gc_interpolate(lon0, lat0, lon1, lat1, npoints):
    """Interpolate on a great circle between two points [lon0, lat0] and [lon1, lat1]
    http://www.movable-type.co.uk/scripts/latlong.html

    :lon0: (float) longitude of point 1 in degrees
    :lat0: (float) longitude of point 1 in degrees
    :lon1: (float) longitude of point 2 in degrees
    :lat1: (float) longitude of point 2 in degrees
    :npoints: (int) number of points for interpolation
    :returns: (numpy array) longitude and latitude

    """
    radius = gc_radius() # km
    frac = np.linspace(0, 1, npoints)
    lon0_r = np.radians(lon0)
    lat0_r = np.radians(lat0)
    lon1_r = np.radians(lon1)
    lat1_r = np.radians(lat1)
    delta = gc_distance(lon0, lat0, lon1, lat1) / radius
    a = np.sin((1 - frac) * delta) / np.sin(delta)
    b = np.sin(frac * delta) / np.sin(delta)
    x = a * np.cos(lat0_r) * np.cos(lon0_r) + b * np.cos(lat1_r) * np.cos(lon1_r)
    y = a * np.cos(lat0_r) * np.sin(lon0_r) + b * np.cos(lat1_r) * np.sin(lon1_r)
    z = a * np.sin(lat0_r) + b * np.sin(lat1_r)
    lat_out = np.arctan2(z, np.sqrt(x**2 + y**2))
    lon_out = np.arctan2(y, x)
    return np.degrees(lon_out), np.degrees(lat_out)

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
