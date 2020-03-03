import os
import sys
import calendar
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from datetime import datetime
from netCDF4 import Dataset, chartostring
from scipy import spatial
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon

#--------------------------------
# MPASMesh
#--------------------------------

class MPASMesh(object):

    """MPASMesh object"""

    def __init__(self, filepath=None):
        assert filepath is not None, 'Please set the path of MPAS mesh file.'
        self.filepath = filepath

    def load(self):
        """read dataset
        """
        out = Dataset(self.filepath, 'r')
        return out

    def get_edge_sign_on_cell(self, mask=None):
        """ Get the sign of edges on cells

        :mask: (numpy array) mask
        :return: (numpy array) sign of edges

        """
        fmesh = self.load()
        if mask is None:
            edgesOnCell = fmesh.variables['edgesOnCell'][:]
            nEdgesOnCell = fmesh.variables['nEdgesOnCell'][:]
            indexToCellID = fmesh.variables['indexToCellID'][:]
        else:
            edgesOnCell = fmesh.variables['edgesOnCell'][mask,:]
            nEdgesOnCell = fmesh.variables['nEdgesOnCell'][mask]
            indexToCellID = fmesh.variables['indexToCellID'][mask]
        cellsOnEdge = fmesh.variables['cellsOnEdge'][:]
        ncell, nec = edgesOnCell.shape
        edge_sign_on_cell = np.zeros([ncell,nec])
        for i, in np.arange(ncell):
            for j in np.arange(nEdgesOnCell[i]):
                idx_e = edgesOnCell[i,j]-1
                if indexToCellID[i] == cellsOnEdge[idx_e, 0]:
                    edge_sign_on_cell[i,j] = -1
                else:
                    edge_sign_on_cell[i,j] = 1
        return edge_sign_on_cell

    def get_edge_sign_on_vertex(self, mask=None):
        """ Get the sign of edges on vertices

        :mask: (numpy array) mask
        :return: (numpy array) sign of edges

        """
        fmesh = self.load()
        if mask is None:
            edgesOnVertex = fmesh.variables['edgesOnVertex'][:]
            indexToVertexID = fmesh.variables['indexToVertexID'][:]
        else:
            edgesOnVertex = fmesh.variables['edgesOnVertex'][mask,:]
            indexToVertexID = fmesh.variables['indexToVertexID'][mask]
        verticesOnEdge = fmesh.variables['verticesOnEdge'][:]
        nvertex, nev = edgesOnVertex.shape
        edge_sign_on_vertex = np.zeros([nvertex, nev])
        for i in np.arange(nvertex):
            for j in np.arange(3):
                idx_e = edgesOnVertex[i,j]-1
                if indexToVertexID[i] == verticesOnEdge[idx_e, 0]:
                    edge_sign_on_vertex[i,j] = -1
                else:
                    edge_sign_on_vertex[i,j] = 1
        return edge_sign_on_vertex

    def get_shortest_path(self, lonP0, latP0, lonP1, latP1, npoint_ref=1, debug_info=False):
        """ Get the shorted path that connects two endpoints.

        :lonP0: (float) Longitude of endpoint 0
        :latP0: (float) Latitude of endpoint 0
        :lonP1: (float) Longitude of endpoint 1
        :latP1: (float) Latitude of endpoint 1
        :npoint_ref: (int) number of reference points along the great circle
        :weight_ref: (float) weight of reference points
        :debug_info: (bool) print out additional debug information if True

        """
        fmesh           = self.load()
        lonVertex       = np.degrees(fmesh.variables['lonVertex'][:])
        latVertex       = np.degrees(fmesh.variables['latVertex'][:])
        lonEdge         = np.degrees(fmesh.variables['lonEdge'][:])
        latEdge         = np.degrees(fmesh.variables['latEdge'][:])
        indexToVertexID = fmesh.variables['indexToVertexID'][:]
        edgesOnVertex   = fmesh.variables['edgesOnVertex'][:]
        verticesOnEdge  = fmesh.variables['verticesOnEdge'][:]

        # find indices of endpoints
        idxP0 = get_index_latlon(lonP0, latP0, lonVertex, latVertex)
        idxP1 = get_index_latlon(lonP1, latP1, lonVertex, latVertex)
        print('Vertex closest to P0: {:8.5f} {:8.5f}'.format(lonVertex[idxP0], latVertex[idxP0]))
        print('Vertex closest to P1: {:8.5f} {:8.5f}'.format(lonVertex[idxP1], latVertex[idxP1]))
        # find reference points
        lon_ref, lat_ref = gc_interpolate(lonVertex[idxP0], latVertex[idxP0], \
                                          lonVertex[idxP1], latVertex[idxP1], npoint_ref+2)
        lon_ref = np.mod(lon_ref[1:-1], 360)
        lat_ref = np.mod(lat_ref[1:-1], 360)
        # initialize path
        out = self.Path(mesh=self,
                        idx_edge=[],
                        idx_vertex=[],
                        sign_edges=[],
                        lon_edge=[],
                        lat_edge=[],
                        lon_vertex=[],
                        lat_vertex=[])
        # loop over reference points, find the path between these points
        idx_sp0 = idxP0
        for i in np.arange(npoint_ref):
            idx_vertex = np.minimum(i,1)
            idx_sp1 = get_index_latlon(lon_ref[i], lat_ref[i], lonVertex, latVertex)
            print(' - Vertex closest to RefP{:d}: {:8.5f} {:8.5f}'.format(i+1, lonVertex[idx_sp1], latVertex[idx_sp1]))
            out_i = self._get_path(idx_sp0, idx_sp1, lonVertex, latVertex, lonEdge, latEdge, \
                                   indexToVertexID, edgesOnVertex, verticesOnEdge, debug_info)
            out = out + out_i
            idx_sp0 = idx_sp1
        # last path, start from end points P1
        out_n = self._get_path(idxP1, idx_sp1, lonVertex, latVertex, lonEdge, latEdge, \
                               indexToVertexID, edgesOnVertex, verticesOnEdge, debug_info)
        out = out + out_n.reverse()
        return out

    def _get_path(self, idxP0, idxP1, lonVertex, latVertex, lonEdge, latEdge, indexToVertexID, edgesOnVertex, verticesOnEdge, debug_info):
        # initialize arrays
        idx_edges_on_path    = []
        idx_vertices_on_path = []
        sign_edges = []
        # start from vertex P0
        idx_vertex_now = idxP0
        # record vortices on path and the indices
        idx_vertices_on_path.append(idx_vertex_now)
        if debug_info:
            print('\nVertex on path ({:d}): {:8.5f} {:8.5f}'.format(idx_vertex_now, lonVertex[idx_vertex_now], latVertex[idx_vertex_now]))

        # continue if not reached P1
        istep = 0
        while idx_vertex_now != idxP1:
            # print the step
            if debug_info:
                print('\nStep {:d}'.format(istep))
            # find the indices of the three edges on vertex
            edge_arr     = edgesOnVertex[idx_vertex_now,:]
            idx_edge_arr = edge_arr-1
            # compute the distance from P1
            dist = []
            idx_tmp = []
            for idx in idx_edge_arr:
                if idx not in idx_edges_on_path:
                    loni = lonEdge[idx]
                    lati = latEdge[idx]
                    dist.append(gc_distance(loni, lati, lonVertex[idxP1], latVertex[idxP1]))
                    idx_tmp.append(idx)
            # print the location of the three edges
            if debug_info:
                print('\nEdges on vertex:')
                for i, idx in enumerate(idx_tmp):
                    print('   Edge {:d} ({:d}): {:8.5f} {:8.5f} ({:10.4f})'.\
                          format(i, idx, lonEdge[idx], latEdge[idx], dist[i]))
            # choose the edge from the three that is closest to vertex P1
            idx_min = np.argmin(dist)
            idx_edge_next = idx_tmp[idx_min]
            # print the edge on path
            if debug_info:
                print('\nEdge on path : [Edge {:d} ({:d})] {:8.5f} {:8.5f}'.\
                      format(idx_min, idx_edge_next, lonEdge[idx_edge_next], latEdge[idx_edge_next]))
            # record edges on path and the indices
            idx_edges_on_path.append(idx_edge_next)
            # find the other vertex on this edge
            vertex_arr = verticesOnEdge[idx_edge_next,:]
            if vertex_arr[0] == indexToVertexID[idx_vertex_now]:
                vertex_next = vertex_arr[1]
                sign_edges.append(-1)
            else:
                vertex_next = vertex_arr[0]
                sign_edges.append(1)
            idx_vertex_next = vertex_next-1
            # record vortices on path and the indices
            idx_vertices_on_path.append(idx_vertex_next)
            if debug_info:
                print('\nVertex on path ({:d}): {:8.5f} {:8.5f}'.\
                      format(idx_vertex_next, lonVertex[idx_vertex_next], latVertex[idx_vertex_next]))
            # move to next vertex
            idx_vertex_now  = idx_vertex_next
            # count steps
            istep += 1

        # create a path on MPAS mesh
        lon_edge = lonEdge[idx_edges_on_path]
        lat_edge = latEdge[idx_edges_on_path]
        lon_vertex = lonVertex[idx_vertices_on_path]
        lat_vertex = latVertex[idx_vertices_on_path]
        out = self.Path(mesh=self,
                        idx_edge=idx_edges_on_path,
                        idx_vertex=idx_vertices_on_path,
                        sign_edges=sign_edges,
                        lon_edge=lon_edge,
                        lat_edge=lat_edge,
                        lon_vertex=lon_vertex,
                        lat_vertex=lat_vertex)
        return out

    def get_closed_path_cell(self, idx_cell):
        """Get the closed path around a cell defined by all the edges on cell.

        :idx_cell: (int) index for cell

        """

        # load mesh
        fmesh        = self.load()
        nedges       = fmesh.variables['nEdgesOnCell'][idx_cell]
        idx_edges    = fmesh.variables['edgesOnCell'][idx_cell,:nedges]-1
        idx_vertices = fmesh.variables['verticesOnCell'][idx_cell,:nedges]-1
        idx_vertices = np.append(idx_vertices, idx_vertices[0])
        lon_edges    = np.degrees(fmesh.variables['lonEdge'][idx_edges])
        lat_edges    = np.degrees(fmesh.variables['latEdge'][idx_edges])
        lon_vertices = np.degrees(fmesh.variables['lonVertex'][idx_vertices])
        lat_vertices = np.degrees(fmesh.variables['latVertex'][idx_vertices])
        sign_edges = np.ones(nedges)
        for i in np.arange(nedges):
            tmp_idx_vertices = fmesh.variables['verticesOnEdge'][idx_edges[i],:]-1
            if tmp_idx_vertices[1] == idx_vertices[i]:
                sign_edges[i] = -1
        # create a path
        out = self.Path(mesh=self,
                        idx_edge=idx_edges,
                        idx_vertex=idx_vertices,
                        sign_edges=sign_edges,
                        lon_edge=lon_edges,
                        lat_edge=lat_edges,
                        lon_vertex=lon_vertices,
                        lat_vertex=lat_vertices)
        return out

    def get_map(self, varname, position='cell', name=None, units=None):
        """Get map for variable

        :varname: (str) variable name
        :position: (str) data position on grid, 'cell' or 'vertex'
        :name: (str) name of variable, optional
        :units: (str) units of variable, optional

        """

        ncdata = self.load().variables[varname]
        if name is None:
            name = ncdata.long_name
        if units is None:
            units = ncdata.units
        out = MPASOMap(data=ncdata[:], mesh=self, position=position, name=name, units=units)
        return out

    def plot_edges(self, m, **kwargs):
        fmesh = self.load()
        lonEdge         = np.degrees(fmesh.variables['lonEdge'][:])
        latEdge         = np.degrees(fmesh.variables['latEdge'][:])
        lonVertex       = np.degrees(fmesh.variables['lonVertex'][:])
        latVertex       = np.degrees(fmesh.variables['latVertex'][:])
        verticesOnEdge  = fmesh.variables['verticesOnEdge'][:]
        lines = []
        lonmax = np.mod(m.lonmax, 360)
        lonmin = np.mod(m.lonmin, 360)
        latmax = m.latmax
        latmin = m.latmin
        idx = (lonEdge <= lonmax) & (lonEdge >= lonmin) & \
              (latEdge <= latmax) & (latEdge >= latmin)
        verticesOnEdge_arr = verticesOnEdge[idx,:]
        nedges=verticesOnEdge_arr.shape[0]
        for i in np.arange(nedges):
            idx_vertex0 = verticesOnEdge_arr[i,0]-1
            idx_vertex1 = verticesOnEdge_arr[i,1]-1
            lonP0 = lonVertex[idx_vertex0]
            latP0 = latVertex[idx_vertex0]
            lonP1 = lonVertex[idx_vertex1]
            latP1 = latVertex[idx_vertex1]
            x0, y0 = m(lonP0, latP0)
            x1, y1 = m(lonP1, latP1)
            lines.append([(x0, y0), (x1, y1)])
        lc = LineCollection(lines, **kwargs)
        out = m.ax.add_collection(lc)
        return out

    def plot_edges_xy(self, axis=None, **kwargs):
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # load edge data
        fmesh = self.load()
        xCell           = fmesh.variables['xCell'][:]
        yCell           = fmesh.variables['yCell'][:]
        xVertex         = fmesh.variables['xVertex'][:]
        yVertex         = fmesh.variables['yVertex'][:]
        verticesOnEdge  = fmesh.variables['verticesOnEdge'][:]
        lines = []
        nedges = verticesOnEdge.shape[0]
        for i in np.arange(nedges):
            idx_vertex0 = verticesOnEdge[i,0]-1
            idx_vertex1 = verticesOnEdge[i,1]-1
            xP0 = xVertex[idx_vertex0]
            yP0 = yVertex[idx_vertex0]
            xP1 = xVertex[idx_vertex1]
            yP1 = yVertex[idx_vertex1]
            lines.append([(xP0, yP0), (xP1, yP1)])
        lc = LineCollection(lines, **kwargs)
        out = axis.add_collection(lc)
        return out

    def plot_celledges_xy(self, axis=None):
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # load edge data
        fmesh = self.load()
        xCell           = fmesh.variables['xCell'][:]
        yCell           = fmesh.variables['yCell'][:]
        xVertex         = fmesh.variables['xVertex'][:]
        yVertex         = fmesh.variables['yVertex'][:]
        verticesOnCell  = fmesh.variables['verticesOnCell'][:]
        nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
        dvEdge          = fmesh.variables['dvEdge'][:]
        dvEdge_small = 1.0
        dvEdge_max = dvEdge.max() + dvEdge_small
        x_period = fmesh.x_period
        y_period = fmesh.y_period
        patches = []
        ncells = nEdgesOnCell.shape[0]
        for i in np.arange(ncells):
            nedges = nEdgesOnCell[i]
            idx_v = verticesOnCell[i,:nedges]-1
            xp = xVertex[idx_v]
            yp = yVertex[idx_v]
            if any(np.abs(xp[0:-1]-xp[1:]) > dvEdge_max) or \
               any(np.abs(yp[0:-1]-yp[1:]) > dvEdge_max):
                xc = xCell[i]
                yc = yCell[i]
                for j in np.arange(nedges):
                    if xp[j] - xc > dvEdge_max:
                        xp[j] -= x_period
                    elif xp[j] - xc < -dvEdge_max:
                        xp[j] += x_period
                    if yp[j] - yc > dvEdge_max:
                        yp[j] -= y_period
                    elif yp[j] - yc < -dvEdge_max:
                        yp[j] += y_period
            patches.append(Polygon(list(zip(xp,yp))))
        pc = PatchCollection(patches, facecolors='none', edgecolors='black', alpha=1.0)
        pc.set_lw(1)
        out = axis.add_collection(pc)
        return out

    class Path(object):

        """Path on MPASMesh object"""

        def __init__(self, mesh=None, idx_edge=[], idx_vertex=[], sign_edges=[], lon_edge=[], lat_edge=[], lon_vertex=[], lat_vertex=[]):

            self.mesh = mesh
            self.idx_edge = list(idx_edge)
            self.idx_vertex = list(idx_vertex)
            self.sign_edges = list(sign_edges)
            self.lon_edge = list(lon_edge)
            self.lat_edge = list(lat_edge)
            self.lon_vertex = list(lon_vertex)
            self.lat_vertex = list(lat_vertex)

        def __add__(self, other):

            if len(self.idx_vertex) == 0:
                idx_v = 0
            else:
                idx_v = 1
            for attr in self.__dict__.keys():
                attr_val = getattr(self, attr)
                if isinstance(attr_val, list):
                    if 'vertex' in attr:
                        attr_val.extend(getattr(other, attr)[idx_v:])
                    else:
                        attr_val.extend(getattr(other, attr))
                    setattr(self, attr, attr_val)
            return self

        def reverse(self):

            for attr in self.__dict__.keys():
                attr_val = getattr(self, attr)
                if isinstance(attr_val, list):
                    attr_val.reverse()
                    setattr(self, attr, attr_val)
            self.sign_edges = [-1*val for val in self.sign_edges]
            return self


        def plot_edge_center(self, m, s=1, **kwargs):
            x_e, y_e = m(self.lon_edge, self.lat_edge)
            out = m.scatter(x_e, y_e, s=s, **kwargs)
            return out

        def plot_vertex(self, m, s=1, **kwargs):
            x_v, y_v = m(self.lon_vertex, self.lat_vertex)
            out = m.scatter(x_v, y_v, s=s, **kwargs)
            return out

        def plot_edge(self, m, **kwargs):
            x_e, y_e = m(self.lon_vertex, self.lat_vertex)
            out = m.plot(x_e, y_e, **kwargs)
            return out

#--------------------------------
# MPASOData
#--------------------------------

class MPASOData(object):

    """MPASOData object"""

    def __init__(self, filepath=None, filepath_mesh=None):
        """Initialize MPASOData

        :filepath: (str) path of MPASO data file
        :filepath_mesh: (str) path of corresponding mesh file

        """
        assert filepath is not None, 'Please set the path of MPAS data file.'
        assert filepath_mesh is not None, 'Please set the path of corresponding mesh file.'
        self.filepath = filepath
        self.mesh = MPASMesh(filepath=filepath_mesh)

    def load(self):
        """read dataset

        :return: (netCDF4.Dataset) netcdf data

        """
        out = Dataset(self.filepath, 'r')
        return out

    def get_volume(self, varname, name=None, units=None, tidx=0):
        """Get volume for variable

        :varname: (str) variable name
        :name: (str) name of variable, optional
        :units: (str) units of variable, optional
        :return: (MPASOVolume object) volume

        """
        ncdata = self.load().variables[varname]
        if name is None:
            name = ncdata.long_name
        if units is None:
            units = ncdata.units
        out = MPASOVolume(data=ncdata[tidx,:,:], mesh=self.mesh, name=name, units=units)
        return out

    def get_map(self, varname, position='cell', name=None, units=None, tidx=0):
        """Get map for variable

        :varname: (str) variable name
        :position: (str) cell or vertex
        :name: (str) name of variable, optional
        :units: (str) units of variable, optional
        :tidx: (int) time index
        :return: (MPASOMap object) map

        """

        ncdata = self.load().variables[varname]
        if name is None:
            name = ncdata.long_name
        if units is None:
            units = ncdata.units
        out = MPASOMap(data=ncdata[tidx,:], mesh=self.mesh, position=position, name=name, units=units)
        return out

    def get_domain(self, varname, position='cell', name=None, units=None, tidx=0):
        """Get domain for variable

        :varname: (str) variable name
        :position: (str) cell or vertex
        :name: (str) name of variable, optional
        :units: (str) units of variable, optional
        :tidx: (int) time index
        :return: (MPASODomain object) domain

        """

        ncdata = self.load().variables[varname]
        if name is None:
            name = ncdata.long_name
        if units is None:
            units = ncdata.units
        ndim = np.ndim(ncdata)
        if ndim == 1:
            out = MPASODomain(data=ncdata[:], mesh=self.mesh, position=position, name=name, units=units)
        elif ndim == 2:
            out = MPASODomain(data=ncdata[tidx,:], mesh=self.mesh, position=position, name=name, units=units)
        else:
            out = MPASODomain(data=ncdata[tidx,:,:], mesh=self.mesh, position=position, name=name, units=units)
        return out

    def get_profile(self, varname, idx, position='cell', name=None, units=None,
                    tidx_start=None, tidx_end=None):
        """Get profile for variable

        :varname: (str) variable name
        :position: (str) cell or vertex
        :idx: (int) index for cell or vertex
        :name: (str) name of variable, optional
        :units: (str) units of variable, optional
        :tidx_start: (int) starting index
        :tidx_end: (int) ending index
        :return: (MPASOProfile object) profile

        """
        # data
        fdata = self.load()
        ncdata = fdata.variables[varname]
        if varname == 'temperatureLES':
            data = ncdata[tidx_start:tidx_end,idx,:] - 273.15
        else:
            data = ncdata[tidx_start:tidx_end,idx,:]
        if name is None:
            name = ncdata.long_name
        if units is None:
            units = ncdata.units
        ndim = np.ndim(ncdata)
        assert ndim == 3, 'Cannot get profile for {:d}-dimensional variable. Stop.'.format(ndim)
        # time
        xtime = chartostring(fdata.variables['xtime'][tidx_start:tidx_end,:])
        time = [datetime.strptime(x.strip().replace('0000', '0001'), '%Y-%m-%d_%H:%M:%S') for x in xtime]
        # z
        fmesh = self.mesh.load()
        if 'LES' in position:
            depth = fmesh.variables['zLES'][0,idx,:]
        else:
            if 'refZMid' in fmesh.variables.keys():
                depth = fmesh.variables['refZMid'][:]
            elif 'zMid' in fmesh.variables.keys():
                depth = fmesh.variables['zMid'][0,idx,:]
            else:
                raise KeyError('Neither refZMid or zMid is found.' )
        # MPASOProfile
        out = MPASOProfile(time=time, time_name='Time', time_units=None,
                           z=depth, z_name='z', z_units='m',
                           data=data, data_name=name, data_units=units)
        return out


    def get_map_relative_vorticity(self, depth=0.0, mask=None, varname_prefix=''):
        """Get map of relative vorticity

        """
        fdata = self.load()
        fmesh = self.mesh.load()
        data = fdata.variables[varname_prefix+'normalVelocity'][:]
        dcEdge = fmesh.variables['dcEdge'][:]
        verticesOnEdge = fmesh.variables['verticesOnEdge'][:]
        if mask is None:
            edgesOnVertex = fmesh.variables['edgesOnVertex'][:]
            areaTriangle = fmesh.variables['areaTriangle'][:]
        else:
            edgesOnVertex = fmesh.variables['edgesOnVertex'][mask,:]
            areaTriangle = fmesh.variables['areaTriangle'][mask]
        edgeSignOnVertex = self.mesh.get_edge_sign_on_vertex(mask=mask)
        nt, ne, nz = data.shape
        refbottomdepth = fmesh.variables['refBottomDepth'][:]
        nvertlevels = len(refbottomdepth)
        reftopdepth = np.zeros(nvertlevels)
        reftopdepth[1:nvertlevels] = refbottomdepth[0:nvertlevels-1]
        refmiddepth = 0.5*(reftopdepth+refbottomdepth)
        zidx = np.argmin(np.abs(refmiddepth-depth))
        normal_velocity = data[:,:,zidx]
        nv = edgesOnVertex.shape[0]
        vorticity = np.zeros([nt, nv])
        for i in np.arange(nt):
            for idx_v in np.arange(nv):
                idx_ev = edgesOnVertex[idx_v,:]-1
                vorticity[i,idx_v] = np.sum(normal_velocity[i,idx_ev] * dcEdge[idx_ev] * edgeSignOnVertex[idx_v,:])
                vorticity[i,idx_v] = vorticity[i,idx_v]/areaTriangle[idx_v]
        return vorticity


    def get_transport(self, transect=None, path=None, varname=None, varname_prefix='', bolus=False, \
                      refval_diff=None, refval_ratio=None):
        """Compute the transport of variable across transect

        :transect: (VerticalTransect object) transect
        :path: (MPASMesh.Path) path of the transect
        :varname: (str) variable name
        :varname_prefix: (str) variable name prefix (e.g., timeMonthly_avg_)
        :bolus: (bool) Add GM bolus velocity if True
        :refval_diff: (float) Reference value (v-vref)
        :refval_ratio: (float) Reference value (1-v/vref)
        :return: (float) transport

        """

        tran, dist = self.get_transport_cumulative(transect=transect,
                                                   path=path,
                                                   varname=varname,
                                                   varname_prefix=varname_prefix,
                                                   bolus=bolus,
                                                   refval=refval)

        transport = tran[:,-1]
        return transport

    def get_transport_cumulative(self, transect=None, path=None, varname=None, varname_prefix='', bolus=False,\
                                 refval_diff=None, refval_ratio=None):
        """Compute the cumulative transport of variable across transect

        :transect: (VerticalTransect object) transect
        :path: (MPASMesh.Path) path of the transect
        :varname: (str) variable name
        :varname_prefix: (str) variable name prefix (e.g., timeMonthly_avg_)
        :bolus: (bool) Add GM bolus velocity if True
        :refval_diff: (float) Reference value (v-vref)
        :refval_ratio: (float) Reference value (1-v/vref)
        :return: ([float, float]) transport, distance

        """

        fdata = self.load()
        fmesh = self.mesh.load()
        if path is None:
            assert transect is not None, 'Transect \'transect\' required if \'path\' is not set.'
            path = self.mesh.get_shortest_path(transect.lon0, transect.lat0, transect.lon1, transect.lat1)
        idx_edge = path.idx_edge
        sign_edges_1d = path.sign_edges
        dv_edge_1d = fmesh.variables['dvEdge'][idx_edge]
        normal_velocity = fdata.variables[varname_prefix+'normalVelocity'][:,idx_edge,:]
        if bolus:
            normal_gm_bolus_velocity = fdata.variables[varname_prefix+'normalGMBolusVelocity'][:,idx_edge,:]
            normal_velocity = normal_velocity + normal_gm_bolus_velocity
        cells_on_edge = fmesh.variables['cellsOnEdge'][idx_edge,:]
        idx_cells_on_edge = cells_on_edge - 1
        if varname is None:
            var = np.ones(normal_velocity.shape)
        else:
            var_c0 = fdata.variables[varname_prefix+varname][:,idx_cells_on_edge[:,0],:]
            var_c1 = fdata.variables[varname_prefix+varname][:,idx_cells_on_edge[:,1],:]
            var = 0.5*(var_c0+var_c1)
            if refval_diff is not None:
                var = var-refval_diff
            if refval_ratio is not None:
                var = 1.0-var/refval_ratio
        layer_thickness_c0 = fdata.variables[varname_prefix+'layerThickness'][:,idx_cells_on_edge[:,0],:]
        layer_thickness_c1 = fdata.variables[varname_prefix+'layerThickness'][:,idx_cells_on_edge[:,1],:]
        dh_edge = 0.5*(layer_thickness_c0+layer_thickness_c1)
        nt, ne, nz = normal_velocity.shape
        dv_edge = np.transpose(np.tile(dv_edge_1d,[nz,1]))
        sign_edges = np.transpose(np.tile(sign_edges_1d,[nz,1]))
        transport = np.zeros([nt, ne+1])
        for i in np.arange(nt):
            tmp = normal_velocity[i,:,:]*var[i,:,:]*sign_edges*dv_edge*dh_edge[i,:,:]
            transport[i, 1:] = np.cumsum(np.sum(tmp, axis=1))
        dist = np.zeros(ne+1)
        for j in np.arange(ne)+1:
            dist[j] = gc_distance(path.lon_vertex[j], path.lat_vertex[j], path.lon_vertex[0], path.lat_vertex[0])
        return transport, dist


#--------------------------------
# MPASOVolume
#--------------------------------

class MPASOVolume(object):

    """MPASOVolume object"""

    def __init__(self, data=None, lon=None, lat=None, depth=None, cellarea=None, \
                 layerthickness=None, bottomdepth=None, position='cell', name=None, units=None, mesh=None):
        """Iniitalize MPASOVolume

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :depth: (1D numpy array) depth
        :cellarea: (1D numpy array) area of cells
        :layerthickness: (1D numpy array) layer thickness
        :bottomdepth: (1D numpy array) depth of bottom
        :position: (str) data position on grid, 'cell' or 'vertex'
        :name: (str) name of variable
        :units: (str) units of variable
        :mesh: (MPASMesh) mesh object

        """
        assert data is not None, 'Data array \'data\' required.'
        self.fillvalue = -9.99999979021476795361e+33
        self.data = np.where(data<=self.fillvalue, np.nan, data)
        self.name = name
        self.units = units
        self.mesh = mesh
        self.position = position
        if mesh is None:
            assert lon is not None, 'Longitude array \'lon\' required.'
            assert lat is not None, 'Latitude array \'lat\' required.'
            assert depth is not None, 'Depth array \'depth\' required.'
            assert cellarea is not None, 'Cell area array \'cellarea\' required.'
            assert layerthickness is not None, 'Layer thickness array \'layerthickness\' required.'
            assert bottomdepth is not None, 'Bottom depth array \'bottomdepth\' required.'
            self.lon = lon
            self.lat = lat
            self.depth = depth
            self.cellarea = cellarea
            self.layerthickness = layerthickness
            self.bottomdepth = bottomdepth
        else:
            print("Reading mesh data from {}".format(mesh.filepath))
            fmesh = mesh.load()
            self.lon = np.degrees(fmesh.variables['lonCell'][:])
            self.lat = np.degrees(fmesh.variables['latCell'][:])
            self.cellarea = fmesh.variables['areaCell'][:]
            self.bottomdepth = fmesh.variables['bottomDepth'][:]
            refbottomdepth = fmesh.variables['refBottomDepth'][:]
            nvertlevels = len(refbottomdepth)
            reftopdepth = np.zeros(nvertlevels)
            reftopdepth[1:nvertlevels] = refbottomdepth[0:nvertlevels-1]
            reflayerthickness = reftopdepth-refbottomdepth
            refmiddepth = 0.5*(reftopdepth+refbottomdepth)
            self.depth = refmiddepth
            self.layerthickness = reflayerthickness

    def get_map(self, depth=0.0):
        """ Return a map at depth

        :depth: (float) depth of the cross section
        :returns: (MPASOMap object) map

        """
        # get z index
        zidx = np.argmin(np.abs(self.depth-depth))
        # MPASOMap object
        name = self.name+' at {:6.2f} m'.format(self.depth[zidx])
        obj = MPASOMap(data=self.data[:,zidx], mesh=self.mesh, lon=self.lon, lat=self.lat, cellarea=self.cellarea, position=self.position, name=name, units=self.units)
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
        obj = MPASOMap(data=data, mesh=self.mesh, lon=self.lon, lat=self.lat, cellarea=self.cellarea, position=self.position, name=name, units=units)
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
        obj = MPASOMap(data=data, mesh=self.mesh, lon=self.lon, lat=self.lat, cellarea=self.cellarea, position=self.position, name=name, units=units)
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
        bottomdepth = self.bottomdepth[cidx]
        # calculate distance from [lon0, lat0]
        dist = np.zeros(npoints)
        for i in np.arange(npoints-1):
            dist[i+1] = gc_distance(lon0, lat0, loni[i+1], lati[i+1])
        obj = MPASOVertCrossSection(data=data, lon=loni, lat=lati, dist=dist, depth=depth, bottomdepth=bottomdepth, name=self.name, units=self.units)
        return obj

    def get_transect(self, transect):
        """Return the transect defined by the VerticalTransect object.

        :transect: (VerticalTransect object) vertical transect
        :returns: (MPASOVertCrossSection object) vertical cross section

        """
        obj = self.get_vertical_cross_section(lon0=transect.lon0,
                                              lat0=transect.lat0,
                                              lon1=transect.lon1,
                                              lat1=transect.lat1,
                                              depth_bottom=transect.depth,
                                              depth_top=0.0)
        return obj

#--------------------------------
# MPASOMap
#--------------------------------

class MPASOMap(object):

    """MPASOMap object"""

    def __init__(self, data=None, lon=None, lat=None, cellarea=None, position='cell', name=None, units=None, mesh=None):
        """Initialize MPASOMap

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :cellarea: (1D numpy array) area of cells
        :position: (str) data position on grid, 'cell' or 'vertex'
        :name: (str) name of variable
        :units: (str) units of variable
        :mesh: (MPASMesh) mesh object

        """
        if data is not None:
            self.fillvalue = -9.99999979021476795361e+33
            self.data = np.where(data<=self.fillvalue, np.nan, data)
            self.name = name
            self.units = units
            self.mesh = mesh
            self.position = position
            if mesh is None:
                assert lon is not None, 'Longitude array \'lon\' required.'
                assert lat is not None, 'Latitude array \'lat\' required.'
                assert cellarea is not None, 'Cell area array \'cellarea\' required.'
                self.lon = lon
                self.lat = lat
                self.cellarea = cellarea
            else:
                print("Reading mesh data from {}".format(mesh.filepath))
                fmesh = mesh.load()
                if position == 'cell':
                    self.lon = np.degrees(fmesh.variables['lonCell'][:])
                    self.lat = np.degrees(fmesh.variables['latCell'][:])
                    self.cellarea = fmesh.variables['areaCell'][:]
                elif position == 'vertex':
                    self.lon = np.degrees(fmesh.variables['lonVertex'][:])
                    self.lat = np.degrees(fmesh.variables['latVertex'][:])
                    self.cellarea = fmesh.variables['areaTriangle'][:]
                else:
                    raise ValueError('Unsupported position \'{}\''.format(position))

    def __sub__(self, other):
        """Subtract 'other' from an MPASOMap object

        :other: (float, int, or MPASOMap object): object to be subtracted
        :returns: (MPASOMap object) the modified MPASOMap object

        """
        if isinstance(other, float) or isinstance(other, int):
            out = MPASOMap(data=self.data-other, lon=self.lon, lat=self.lat, cellarea=self.cellarea, \
                           position=self.position, name=self.name, units=self.name, mesh=self.mesh)
        elif isinstance(other, MPASOMap):
            assert self.units == other.units, 'MPASOMap has a different unit'
            assert all(self.lon == other.lon), 'MPASOMap has different Longitude'
            assert all(self.lat == other.lat), 'MPASOMap has different Latitude'
            assert self.position == other.position, 'MPASOMap is on different cell position'
            out = MPASOMap(data=self.data-other.data, lon=self.lon, lat=self.lat, cellarea=self.cellarea, \
                           position=self.position, name=self.name+' $-$ '+other.name, \
                           units=self.units, mesh=self.mesh)
        else:
            raise TypeError('Subtraction is not defined between an MPASOMap object and a {} object'.format(type(other)))
        return out

    def save(self, filepath):
        """Save MPASOMap object

        :filepath: (str) path of file to save
        :returns: none

        """
        np.savez(filepath, data=self.data, lon=self.lon, lat=self.lat, cellarea=self.cellarea, \
                 name=self.name, units=self.units, meshfile=self.mesh.filepath)

    def load(self, filepath):
        """Load data to MPASOMap object

        :filepath: (str) path of file to load
        :returns: (MPASOMap object)

        """
        dat = np.load(filepath)
        if 'meshfile' not in dat.keys():
            mesh = None
        else:
            meshfile = str(dat['meshfile'])
            mesh = MPASMesh(filepath=meshfile)
        self.__init__(data=dat['data'], lon=dat['lon'], lat=dat['lat'], cellarea=dat['cellarea'], \
                      name=str(dat['name']), units=str(dat['units']), mesh=mesh)
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
                    'Provide lon_min, lat_min, lon_max, and lat_max to customize region'
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

    def _pcolor(self, m, position='cell', **kwargs):
        assert self.mesh is not None, 'Mesh file required for _pcolor.'
        fmesh = self.mesh.load()
        lonCell   = np.degrees(fmesh.variables['lonCell'][:])
        latCell   = np.degrees(fmesh.variables['latCell'][:])
        lonVertex = np.degrees(fmesh.variables['lonVertex'][:])
        latVertex = np.degrees(fmesh.variables['latVertex'][:])
        # bounds of the domain
        lonmax = np.mod(m.lonmax, 360)
        lonmin = np.mod(m.lonmin, 360)
        latmax = m.latmax
        latmin = m.latmin
        if position == 'cell':
            verticesOnCell  = fmesh.variables['verticesOnCell'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            idx = (lonCell <= lonmax) & (lonCell >= lonmin) & \
                  (latCell <= latmax) & (latCell >= latmin) & \
                  (~ np.isnan(self.data))
            verticesOnCell_arr = verticesOnCell[idx,:]
            nEdgesOnCell_arr = nEdgesOnCell[idx]
            data = self.data[idx]
            # patches
            patches = []
            ncell = verticesOnCell_arr.shape[0]
            for i in np.arange(ncell):
                idx_v = verticesOnCell_arr[i,:nEdgesOnCell_arr[i]]-1
                lonp = lonVertex[idx_v]
                latp = latVertex[idx_v]
                xp, yp = m(lonp, latp)
                patches.append(Polygon(list(zip(xp,yp))))
        elif position == 'vertex':
            cellsOnVertex = fmesh.variables['cellsOnVertex'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            idx = (lonVertex <= lonmax) & (lonVertex >= lonmin) & \
                  (latVertex <= latmax) & (latVertex >= latmin) & \
                  (~ np.isnan(self.data))
            cellsOnVertex_arr = cellsOnVertex[idx,:]
            data = self.data[idx]
            # patches
            patches = []
            idx_mask = []
            nvertex = cellsOnVertex_arr.shape[0]
            for i in np.arange(nvertex):
                idx_c = cellsOnVertex_arr[i,:]-1
                if any(idx_c == -1):
                    idx_mask.append(i)
                    continue
                lonp = lonCell[idx_c]
                latp = latCell[idx_c]
                xp, yp = m(lonp, latp)
                patches.append(Polygon(list(zip(xp,yp))))
            data = np.delete(data, idx_mask)
        else:
            raise ValueError('Unsupported position \'{}\''.format(position))
        # plot patch collection
        pc = PatchCollection(patches, **kwargs)
        pc.set_array(data)
        pc.set_lw(0.1)
        out = m.ax.add_collection(pc)
        return out

    def _pcolor_nan_mask(self, m, position='cell'):
        assert self.mesh is not None, 'Mesh file required for _pcolor_nan_mask.'
        fmesh = self.mesh.load()
        lonCell         = np.degrees(fmesh.variables['lonCell'][:])
        latCell         = np.degrees(fmesh.variables['latCell'][:])
        lonVertex       = np.degrees(fmesh.variables['lonVertex'][:])
        latVertex       = np.degrees(fmesh.variables['latVertex'][:])
        # bounds of the domain
        lonmax = np.mod(m.lonmax, 360)
        lonmin = np.mod(m.lonmin, 360)
        latmax = m.latmax
        latmin = m.latmin
        if position == 'cell':
            verticesOnCell  = fmesh.variables['verticesOnCell'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            idx = (lonCell <= lonmax) & (lonCell >= lonmin) & \
                  (latCell <= latmax) & (latCell >= latmin) & \
                  (np.isnan(self.data))
            verticesOnCell_arr = verticesOnCell[idx,:]
            nEdgesOnCell_arr = nEdgesOnCell[idx]
            data = self.data[idx]
            # patches
            patches = []
            ncell = verticesOnCell_arr.shape[0]
            for i in np.arange(ncell):
                idx_v = verticesOnCell_arr[i,:nEdgesOnCell_arr[i]]-1
                lonp = lonVertex[idx_v]
                latp = latVertex[idx_v]
                xp, yp = m(lonp, latp)
                patches.append(Polygon(list(zip(xp,yp))))
        elif position == 'vertex':
            cellsOnVertex = fmesh.variables['cellsOnVertex'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            idx = (lonVertex <= lonmax) & (lonVertex >= lonmin) & \
                  (latVertex <= latmax) & (latVertex >= latmin) & \
                  (np.isnan(self.data))
            cellsOnVertex_arr = cellsOnVertex[idx,:]
            data = self.data[idx]
            # patches
            patches = []
            idx_mask = []
            nvertex = cellsOnVertex_arr.shape[0]
            for i in np.arange(nvertex):
                idx_c = cellsOnVertex_arr[i,:]-1
                if any(idx_c == -1):
                    idx_mask.append(i)
                    continue
                lonp = lonCell[idx_c]
                latp = latCell[idx_c]
                xp, yp = m(lonp, latp)
                patches.append(Polygon(list(zip(xp,yp))))
            data = np.delete(data, idx_mask)
        else:
            raise ValueError('Unsupported position \'{}\''.format(position))
        # plot patch collection
        pc = PatchCollection(patches, facecolors='lightgray', edgecolors='lightgray', alpha=1.0)
        pc.set_lw(0.1)
        out = m.ax.add_collection(pc)
        return out

    def plot(self, axis=None, region='Global', ptype='contourf', mask_nan=False, levels=None,
             label=None, add_title=True, title=None, add_colorbar=True, cmap='rainbow', projection=None, **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :region: (str) region name
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :label: (str) label
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :projection: (str, optional) projection type
        :**kwargs: (keyword arguments) other arguments
        :return: (basemap) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # print message
        print('Plotting map of {} at region \'{}\''.format(self.name+' ('+self.units+')', region))
        m = plot_basemap(region=region, axis=axis, projection=projection)
        # mask out nan
        nan_mask = (~ np.isnan(self.data))
        # process data
        if region == 'Global':
            # markersize
            markersize = 1
            # data
            data = self.data[nan_mask]
            lat = self.lat[nan_mask]
            lon = self.lon[nan_mask]
            cellarea = self.cellarea[nan_mask]
            # shift longitude
            lon = np.where(lon < 20., lon+360., lon)
        else:
            # regional map
            region_obj = region_latlon(region)
            lon_ll, lat_ll, lon_ur, lat_ur = region_obj.lon_ll, region_obj.lat_ll, region_obj.lon_ur, region_obj.lat_ur
            # region mask
            lonmax = np.mod(m.lonmax, 360)
            lonmin = np.mod(m.lonmin, 360)
            latmax = m.latmax
            latmin = m.latmin
            lon_mask = (self.lon <= lonmax) & (self.lon >= lonmin)
            lat_mask = (self.lat <= latmax) & (self.lat >= latmin)
            region_mask = lon_mask & lat_mask & nan_mask
            # apply region mask to data
            data = self.data[region_mask]
            lat = self.lat[region_mask]
            lon = self.lon[region_mask]
            cellarea = self.cellarea[region_mask]
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
            if mask_nan:
                self._pcolor_nan_mask(m, position=self.position)
        elif ptype == 'pcolor':
            fig = self._pcolor(m, position=self.position, norm=norm, cmap=plt.cm.get_cmap(cmap), alpha=1.0, **kwargs)
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
            if title is None:
                axis.set_title('{} ({})'.format(self.name, self.units))
            else:
                axis.set_title(title)
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        return m, fig

    def overlay(self, m, axis=None, label=False, label_fmt='%1.2f', **kwargs):
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
        if label:
            axis.clabel(fig, fig.levels, fmt=label_fmt)

#--------------------------------
# MPASOVertCrossSection
#--------------------------------

class MPASOVertCrossSection(object):

    """MPASOVertCrossSection object"""

    def __init__(self, data, lon, lat, dist, depth, name, units, bottomdepth=None):
        """Initialize MPASOCrossSection

        :data: (1D numpy array) data array
        :lon: (1D numpy array) longitude array
        :lat: (1D numpy array) latitude array
        :dist: (1D numpy array) distance array
        :depth: (1D numpy array) depth array
        :bottomdepth: (1D numpy array) depth of bottom
        :name: (str) name of variable
        :units: (str) units of variable

        """
        if bottomdepth is not None:
            for i in np.arange(bottomdepth.size):
                idxd = depth>bottomdepth[i]
                data[i,idxd] = np.nan
        self.data = data
        self.lon = lon
        self.lat = lat
        self.dist = dist
        self.depth = depth
        self.bottomdepth = bottomdepth
        self.name = name
        self.units = units

    def plot(self, axis=None, ptype='contourf', depth_mode='linear', levels=None, add_title=True, add_colorbar=True, invert_yaxis=True, cmap='rainbow', **kwargs):
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
        elif ptype == 'contour':
            fig = axis.contour(self.dist, depth, np.transpose(self.data), levels=levels,
                    norm=norm, **kwargs)
            axis.clabel(fig, fig.levels, fmt='%1.2f')
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
        if invert_yaxis:
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
# MPASODomain
#--------------------------------

class MPASODomain(object):

    """MPASODomain object"""

    def __init__(self, data=None, mesh=None, position='cell', name=None, units=None):
        """Initialize MPASODomain

        :data: (1D numpy array) data at each location
        :mesh: (MPASMesh) mesh object
        :position: (str) data position on grid, 'cell' or 'vertex'
        :name: (str) name of variable
        :units: (str) units of variable

        """
        if data is not None:
            self.data = data
            self.ndim = np.ndim(data)
            self.name = name
            self.units = units
            self.mesh = mesh
            self.position = position
            print("Reading mesh data from {}".format(mesh.filepath))
            fmesh = mesh.load()
            if 'cell' in position:
                self.x = fmesh.variables['xCell'][:]
                self.y = fmesh.variables['yCell'][:]
            elif 'vertex' in position:
                self.x = fmesh.variables['xVertex'][:]
                self.y = fmesh.variables['yVertex'][:]
            else:
                raise ValueError('Unsupported position \'{}\''.format(position))
            if self.ndim == 2:
                if 'LES' in position:
                    self.z = fmesh.variables['zLES'][0,0,:]
                else:
                    if 'refZMid' in fmesh.variables.keys():
                        self.z = fmesh.variables['refZMid'][:]
                    elif 'zMid' in fmesh.variables.keys():
                        self.z = fmesh.variables['zMid'][0,0,:]
                    else:
                        raise KeyError('Neither refZMid or zMid is found.' )

    def _pcolor(self, data, axis=None, position='cell', **kwargs):
        assert self.mesh is not None, 'Mesh file required for _pcolor.'
        fmesh = self.mesh.load()
        xCell   = fmesh.variables['xCell'][:]
        yCell   = fmesh.variables['yCell'][:]
        xVertex = fmesh.variables['xVertex'][:]
        yVertex = fmesh.variables['yVertex'][:]
        dvEdge  = fmesh.variables['dvEdge'][:]
        dvEdge_small = 1.0
        dvEdge_max = dvEdge.max() + dvEdge_small
        x_period = fmesh.x_period
        y_period = fmesh.y_period
        if 'cell' in position:
            verticesOnCell  = fmesh.variables['verticesOnCell'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            # patches
            patches = []
            ncell = verticesOnCell.shape[0]
            for i in np.arange(ncell):
                idx_v = verticesOnCell[i,:nEdgesOnCell[i]]-1
                xp = xVertex[idx_v]
                yp = yVertex[idx_v]
                if any(np.abs(xp[0:-1]-xp[1:]) > dvEdge_max) or \
                   any(np.abs(yp[0:-1]-yp[1:]) > dvEdge_max):
                    xc = xCell[i]
                    yc = yCell[i]
                    for j in np.arange(nEdgesOnCell[i]):
                        if xp[j] - xc > dvEdge_max:
                            xp[j] -= x_period
                        elif xp[j] - xc < -dvEdge_max:
                            xp[j] += x_period
                        if yp[j] - yc > dvEdge_max:
                            yp[j] -= y_period
                        elif yp[j] - yc < -dvEdge_max:
                            yp[j] += y_period
                patches.append(Polygon(list(zip(xp,yp))))
        elif 'vertex' in position:
            cellsOnVertex = fmesh.variables['cellsOnVertex'][:]
            nEdgesOnCell    = fmesh.variables['nEdgesOnCell'][:]
            # patches
            patches = []
            idx_mask = []
            nvertex = cellsOnVertex.shape[0]
            for i in np.arange(nvertex):
                idx_c = cellsOnVertex[i,:]-1
                if any(idx_c == -1):
                    idx_mask.append(i)
                    continue
                xp = xCell[idx_c]
                yp = yCell[idx_c]
                if any(np.abs(xp[0:-1]-xp[1:]) > dvEdge_max) or \
                   any(np.abs(yp[0:-1]-yp[1:]) > dvEdge_max):
                    xc = xVertex[i]
                    yc = yVertex[i]
                    for j in np.arange(3):
                        if xp[j] - xc > dvEdge_max:
                            xp[j] = xp[j] - x_period
                        elif xp[j] - xc < -dvEdge_max:
                            xp[j] = xp[j] + x_period
                        if yp[j] - yc > dvEdge_max:
                            yp[j] = yp[j] - y_period
                        elif yp[j] - yc < -dvEdge_max:
                            yp[j] = yp[j] + y_period
                patches.append(Polygon(list(zip(xp,yp))))
            data = np.delete(data, idx_mask)
        else:
            raise ValueError('Unsupported position \'{}\''.format(position))
        # plot patch collection
        pc = PatchCollection(patches, **kwargs)
        pc.set_array(data)
        if len(patches) > 64:
            pc.set_linewidth(0.1)
        out = axis.add_collection(pc)
        return out

    def plot(self, axis=None, ptype='pcolor', levels=None,
             add_title=True, title=None, add_colorbar=True, cmap='viridis', **kwargs):
        """Plot horizontal map of a domain

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) other arguments
        :return: (basemap) figure

        """
        fig = self.plot_xy(axis=None, ptype='pcolor', levels=None, \
                  add_title=True, title=None, add_colorbar=True, cmap=cmap, **kwargs)
        return fig

    def plot_xy(self, zidx=0, axis=None, ptype='pcolor', levels=None,
             add_title=True, title=None, add_colorbar=True, cmap='viridis', **kwargs):
        """Plot horizontal map of a domain

        :zidx: (int) z index
        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) other arguments
        :return: (basemap) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # manually mapping levels to the colormap if levels is passed in,
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # data
        if self.ndim == 1:
            data = self.data
        elif self.ndim == 2:
            data = self.data[:,zidx]
        else:
            raise ValueError('Dimension mismatch.')
        # plot type
        if ptype == 'pcolor':
            fig = self._pcolor(data=data, axis=axis, position=self.position, norm=norm, \
                               cmap=plt.cm.get_cmap(cmap), alpha=1.0, **kwargs)
        elif ptype == 'contourf':
            fig = axis.tricontourf(self.x, self.y, data, levels=levels, extend='both', \
                                   norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
            # axis.tricontour(self.x, self.y, data, colors='k', levels=levels, extend='both', \
            #                 norm=norm, **kwargs)
        else:
            raise ValueError('Plot type {} not supported.'.format(ptype))
        # add title
        if add_title:
            if title is None:
                axis.set_title('{} ({})'.format(self.name, self.units))
            else:
                axis.set_title(title)
        # add colorbar
        if add_colorbar:
            cb = plt.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        # x- and y-limits and x- and y-labels
        xmax = self.x.max()
        xmin = self.x.min()
        ymax = self.y.max()
        ymin = self.y.min()
        axis.set_xlim([xmin, xmax])
        axis.set_ylim([ymin, ymax])
        axis.set_xlabel('x')
        axis.set_ylabel('y')
        return fig

    def plot_yx(self, zidx=0, axis=None, ptype='contourf', levels=None,
             add_title=True, title=None, add_colorbar=True, cmap='viridis', **kwargs):
        """Plot horizontal map of a domain

        :zidx: (int) z index
        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) other arguments
        :return: (basemap) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # manually mapping levels to the colormap if levels is passed in,
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # data
        if self.ndim == 1:
            data = self.data
        elif self.ndim == 2:
            data = self.data[:,zidx]
        else:
            raise ValueError('Dimension mismatch.')
        # plot type
        if ptype == 'pcolor':
            fig = self._pcolor(data=data, axis=axis, position=self.position, norm=norm, \
                               cmap=plt.cm.get_cmap(cmap), alpha=1.0, **kwargs)
        elif ptype == 'contourf':
            fig = axis.tricontourf(self.y, self.x, data, levels=levels, extend='both', \
                                   norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
            # axis.tricontour(self.y, self.x, np.transpose(data), colors='k', levels=levels, extend='both', \
            #                        norm=norm, **kwargs)
        else:
            raise ValueError('Plot type {} not supported.'.format(ptype))
        # add title
        if add_title:
            if title is None:
                axis.set_title('{} ({})'.format(self.name, self.units))
            else:
                axis.set_title(title)
        # add colorbar
        if add_colorbar:
            cb = plt.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        # x- and y-limits and x- and y-labels
        xmax = self.y.max()
        xmin = self.y.min()
        ymax = self.x.max()
        ymin = self.x.min()
        axis.set_xlim([xmin, xmax])
        axis.set_ylim([ymin, ymax])
        axis.set_xlabel('y')
        axis.set_ylabel('x')
        return fig

    def plot_xz(self, yfrac=0.5, **kwargs):
        """Plot xz-transect of a domain

        :yfrac: (float) normalized y coordinate [0, 1]
        :axis: (matplotlib.axes, optional) axis to plot figure on
        :**kwargs: (keyword arguments) arguments
        :return: (axis) axis of figure

        """
        xmax = self.x.max()
        xmin = self.x.min()
        ymax = self.y.max()
        ymin = self.y.min()
        cellarea = self.mesh.load().variables['areaCell'][:]
        d_cell = np.sqrt(cellarea.mean())
        npoints = int(np.ceil((xmax-xmin)/d_cell))
        print('Nearest neighbor interpolation to {} points.'.format(npoints))
        ymid = yfrac*(ymax-ymin)+ymin
        xx = np.linspace(xmin, xmax, npoints+1)
        yy = np.ones(xx.size)*ymid
        # select nearest neighbor
        pts = np.array(list(zip(xx, yy)))
        tree = spatial.KDTree(list(zip(self.x, self.y)))
        p = tree.query(pts)
        cidx = p[1]
        ax = self._plot_transect(xy=self.x[cidx], data=self.data[cidx,:], **kwargs)
        ax.set_xlabel('x')
        return ax

    def plot_yz(self, xfrac=0.5, **kwargs):
        """Plot yz-transect of a domain

        :xfrac: (float) normalized x coordinate [0, 1]
        :axis: (matplotlib.axes, optional) axis to plot figure on
        :**kwargs: (keyword arguments) arguments
        :return: (axis) axis of figure

        """
        xmax = self.x.max()
        xmin = self.x.min()
        ymax = self.y.max()
        ymin = self.y.min()
        cellarea = self.mesh.load().variables['areaCell'][:]
        d_cell = np.sqrt(cellarea.mean())
        npoints = int(np.ceil((ymax-ymin)/d_cell))
        print('Nearest neighbor interpolation to {} points.'.format(npoints))
        xmid = xfrac*(xmax-xmin)+xmin
        yy = np.linspace(ymin, ymax, npoints+1)
        xx = np.ones(yy.size)*xmid
        # select nearest neighbor
        pts = np.array(list(zip(xx, yy)))
        tree = spatial.KDTree(list(zip(self.x, self.y)))
        p = tree.query(pts)
        cidx = p[1]
        ax = self._plot_transect(xy=self.y[cidx], data=self.data[cidx,:], **kwargs)
        ax.set_xlabel('y')
        return ax

    def plot_xz_mean(self, **kwargs):
        """Plot xz-transect of y-mean of a domain

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :**kwargs: (keyword arguments) arguments
        :return: (axis) axis of figure

        """
        xmax = self.x.max()
        xmin = self.x.min()
        nc = self.data.shape[0]
        nz = self.data.shape[1]
        cellarea = self.mesh.load().variables['areaCell'][:]
        d_cell = np.sqrt(cellarea.mean())
        npoints = int(np.ceil((xmax-xmin)/d_cell))
        print('Average over y at {} bins in x.'.format(npoints))
        binbnd = np.linspace(xmin, xmax, npoints+1)
        binwidth = (xmax-xmin)/npoints
        dsum = np.zeros([npoints, nz])
        dwgt = np.zeros(npoints)
        mdata = np.zeros([npoints, nz])
        for i in np.arange(nc):
            idx = int((self.x[i]-binbnd[0]-1.e-6)/binwidth)
            dsum[idx,:] += self.data[i,:]*cellarea[i]
            dwgt[idx] += cellarea[i]
        for j in np.arange(nz):
            mdata[:,j] = dsum[:,j]/dwgt
        xx = 0.5*(binbnd[0:-1]+binbnd[1:])
        ax = self._plot_transect(xy=xx, data=mdata, **kwargs)
        ax.set_xlabel('x')
        return ax

    def plot_yz_mean(self, **kwargs):
        """Plot yz-transect of x-mean of a domain

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :**kwargs: (keyword arguments) arguments
        :return: (axis) axis of figure

        """
        ymax = self.y.max()
        ymin = self.y.min()
        nc = self.data.shape[0]
        nz = self.data.shape[1]
        cellarea = self.mesh.load().variables['areaCell'][:]
        d_cell = np.sqrt(cellarea.mean())
        npoints = int(np.ceil((ymax-ymin)/d_cell))
        print('Average over x at {} bins in y.'.format(npoints))
        binbnd = np.linspace(ymin, ymax, npoints+1)
        binwidth = (ymax-ymin)/npoints
        dsum = np.zeros([npoints, nz])
        dwgt = np.zeros(npoints)
        mdata = np.zeros([npoints, nz])
        for i in np.arange(nc):
            idx = int((self.y[i]-binbnd[0]-1.e-6)/binwidth)
            dsum[idx,:] += self.data[i,:]*cellarea[i]
            dwgt[idx] += cellarea[i]
        for j in np.arange(nz):
            mdata[:,j] = dsum[:,j]/dwgt
        yy = 0.5*(binbnd[0:-1]+binbnd[1:])
        ax = self._plot_transect(xy=yy, data=mdata, **kwargs)
        ax.set_xlabel('y')
        return ax

    def _plot_transect(self, xy=None, data=None, axis=None, ptype='contourf', levels=None, add_title=True, \
                      title=None, add_colorbar=True, cmap='viridis', **kwargs):
        """Plot transect of a domain

        :xy: (numpy array) horizontal dimension
        :data: (int) xy transect data
        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, scatter, contourf etc.
        :leveles: (list, optional) list of levels
        :add_title: (bool) do not add title if False
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) other arguments
        :return: (axis) axis of figure

        """
        # check dimension
        assert self.ndim == 2, '2D domain has no transect.'
        # check input
        assert xy is not None, 'horizontal dimension xy required.'
        assert data is not None, 'transect data required.'
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # manually mapping levels to the colormap if levels is passed in,
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # plot figure
        if ptype == 'contourf':
            fig = axis.contourf(xy, self.z, np.transpose(data), levels=levels, extend='both', \
                                norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
            # axis.contour(xy, self.z, np.transpose(data), colors='k', levels=levels, extend='both', \
            #                     norm=norm, **kwargs)
        elif ptype == 'pcolor':
            fig = axis.pcolor(xy, self.z, np.transpose(data), \
                              norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            raise ValueError('Plot type {} not supported.'.format(ptype))
        # add title
        if add_title:
            if title is None:
                axis.set_title('{} ({})'.format(self.name, self.units))
            else:
                axis.set_title(title)
        # add colorbar
        if add_colorbar:
            cb = plt.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-4, 4))
            cb.update_ticks()
        # add y-label
        axis.set_ylabel('Depth (m)')
        return axis

#--------------------------------
# MPASOProfile
#--------------------------------

class MPASOProfile(object):

    """MPASOProfile object"""

    def __init__(self, time=None, time_name='Time', time_units='s',
                       z=None, z_name='z', z_units='m',
                       data=None, data_name=None, data_units=None):
        """Initialize MPASOProfile

        :time: (1D numpy array/datetime object) time
        :time_name: (str, optional) name of time
        :time_units: (str, optional) units of time
        :z: (1D numpy array) vertical coordinate
        :z_name: (str) name of z
        :z_units: (str) units of z
        :data: (2D numpy array) data at each time and z
        :data_name: (str) name of variable
        :data_units: (str) units of variable

        """
        self.time = time
        self.time_name = time_name
        self.time_units = time_units
        self.z = z
        self.z_name = z_name
        self.z_units = z_units
        self.data = data
        self.data_name = data_name
        self.data_units = data_units
        # try:
        #     self.data_mean = np.mean(data, axis=0)
        # except TypeError:
        #     self.data_mean = None

    def save(self, path):
        """Save MPASOProfile object

        :path: (str) path of file to save
        :returns: none

        """
        np.savez(path, data=self.data, \
                 data_name=self.data_name, data_units=self.data_units, \
                 time=self.time, time_name=self.time_name, time_units=self.time_units, \
                 z=self.z, z_name=self.z_name, z_units=self.z_units)

    def load(self, path):
        """Load data to LESProfile object

        :path: (str) path of file to load
        :returns: (MPASOProfile object)

        """
        dat = np.load(path)
        self.__init__(data=dat['data'], \
                      data_name=str(dat['data_name']), data_units=str(dat['data_units']), \
                      time=dat['time'], time_name=str(dat['time_name']), time_units=str(dat['time_units']),\
                      z=dat['z'], z_name=str(dat['z_name']), z_units=str(dat['z_units']))
        return self

    def ddt(self):
        """Return the time derivative of MPASOProfile
        """
        nt = len(self.time)
        nz = len(self.z)
        dtime = np.zeros(nt-1)
        data = np.zeros([nt-1, nz], type(self.data))
        for i in np.arange(nt-1):
            dtime = (self.time[i+1] - self.time[i]).total_seconds()
            data[i,:] = (self.data[i+1,:]-self.data[i,:])/dtime
        # MPASOProfile
        out = MPASOProfile(time=self.time[1:], time_name=self.time_name, time_units=self.time_units,
                           z=self.z, z_name=self.z_name, z_units=self.z_units,
                           data=data, data_name='ddt '+self.data_name, data_units=self.data_units+'/s')
        return out

    def plot(self, axis=None, xlim=None, ylim=None,
                   xlabel=None, ylabel=None, title=None,
                   ptype='contourf', levels=None, cmap='viridis', **kwargs):
        """Plot the Hovmoller diagram (time - z)

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'Time' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :ptype: (str, optional) plot type, valid values: contourf (default), pcolor
        :leveles: (list, optional) list of levels
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.contourf() or
                                       matplotlib.pyplot.pcolor() depending on ptype
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        if levels is not None:
            # manually mapping levels to the colormap if levels is passed in,
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        else:
            norm = None
        # plot type
        if ptype == 'contourf':
            fig = axis.contourf(self.time, self.z, np.transpose(self.data), levels=levels, extend='both', \
                                norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        elif ptype == 'pcolor':
            fig = axis.pcolor(self.time, self.z, np.transpose(self.data), \
                              norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            raise ValueError('Plot type (ptype) should be \'contourf\' or \'pcolor\', got {}.'.format(ptype))
        # set title
        if title is not None:
            axis.set_title(title)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            if self.time_units is not None:
                axis.set_xlabel(self.time_name+' ('+self.time_units+')')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.z_name+' ('+self.z_units+')')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
            axis.set_ylim(ylim)
        # return figure
        return fig

    def plot_mean(self, axis=None, norm=1.0, znorm=1.0, xlim=None, ylim=None,
                        xlabel=None, ylabel=None, title=None, **kwargs):
        """Plot the mean profile

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :norm: (float) normalizing factor
        :znorm: (float) normalizing factor for vertical coordinate
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'self.name' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.plot()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot figure
        fig = axis.plot(self.data_mean*norm, self.z*znorm, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel(self.data_name+' ('+self.data_units+')')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.z_name+' ('+self.z_units+')')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
            axis.set_ylim(ylim)
        # return figure
        return fig


#--------------------------------
# MPASCICEData
#--------------------------------

class MPASCICEData(MPASOData):

    """MPASCICEData object"""
    pass

#--------------------------------
# MPASCICEMap
#--------------------------------

class MPASCICEMap(MPASOMap):

    """MPASCICEMap object"""
    pass

#--------------------------------
# Vertical transsect object
#--------------------------------
class VerticalTransect(object):

    """Vertical transect along the great circle defined by two endpoints"""

    def __init__(self, name=None, lon0=None, lat0=None, lon1=None, lat1=None, depth=None):
        """Initialize VerticalTransect

        :name: (str) trasect name
        :lon0: (float) longitude of endpoint 0
        :lat0: (float) latitude of endpoint 0
        :lon1: (float) longitude of endpoint 1
        :lat1: (float) latitude of endpoint 1
        :depth: (float) maximum depth of transect

        """
        self.name  = name
        if name == 'AR7W':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 304
            self.lat0  = 53.5
            self.lon1  = 312
            self.lat1  = 61
            self.depth = 4500.0
        elif name == 'Davis Strait':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 298.5
            self.lat0  = 66.5
            self.lon1  = 306
            self.lat1  = 67
            self.depth = 1500.0
        elif name == 'Hudson Strait':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 295.2
            self.lat0  = 60.4
            self.lon1  = 293.7
            self.lat1  = 61.9
            self.depth = 1000.0
        elif name == 'Nares Strait':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 284.2
            self.lat0  = 78.0
            self.lon1  = 287.5
            self.lat1  = 78.0
            self.depth = 1000.0
        elif name == 'Parry Channel' or name == 'Lancaster sound':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 281.2
            self.lat0  = 73.7
            self.lon1  = 279.7
            self.lat1  = 74.6
            self.depth = 1000.0
        elif name == 'Jones Sound':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 279.5
            self.lat0  = 75.6
            self.lon1  = 280
            self.lat1  = 76.2
            self.depth = 1000.0
        elif name == 'LabSea Center':
            print('Pre-defined transect \'{}\'.'.format(name))
            self.lon0  = 296
            self.lat0  = 63
            self.lon1  = 320
            self.lat1  = 50
            self.depth = 4500.0
        else:
            print('User defined transect \'{}\'.'.format(name))
            self.lon0  = lon0
            self.lat0  = lat0
            self.lon1  = lon1
            self.lat1  = lat1
            self.depth = depth

    def interpolate(self, npoints):
        """Interpolate along great circle."""
        lon_arr, lat_arr = gc_interpolate(self.lon0, self.lat0, self.lon1, self.lat1, npoints)
        return lon_arr, lat_arr

    def direction(self, npoints):
        """Direction of great circle (P0->P1) defined by the angle (in radius) counterclockwise from East """
        lon, lat = self.interpolate(npoints)
        lat0 = np.zeros(npoints)
        lon0 = np.zeros(npoints)
        lat1 = np.zeros(npoints)
        lon1 = np.zeros(npoints)
        lat0[1:-1] = lat[0:-2]
        lat1[1:-1] = lat[2:]
        lon0[1:-1] = lon[0:-2]
        lon1[1:-1] = lon[2:]
        dir_arr = gc_angle(lon0, lat0, lon1, lat1)
        dir_arr[0] = dir_arr[1]
        dir_arr[-1] = dir_arr[-2]
        return dir_arr

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
    elif region_name == 'Greenland':
        rg = region(lon_ll=296.0, lat_ll=46.0, lon_ur=356.0, lat_ur=85.0)
    elif region_name == 'test':
        rg = region(lon_ll=310.0, lat_ll=55.0, lon_ur=320.0, lat_ur=65.0)
    elif region_name == 'TropicalPacific':
        rg = region(lon_ll=130.0, lat_ll=-20.0, lon_ur=290.0, lat_ur=20.0)
    else:
        raise ValueError('Region {} not supported.'.format(region_name))
    return rg

def get_index_latlon(loni, lati, lon_arr, lat_arr, search_range=5.0):
    """Get the index of the location (loni, lati) in an array of
       locations (lon_arr, lat_arr)

    :loni: (float) Longitude of target location
    :lati: (float) Latitude of target location
    :lon_arr: (float) array of Longitude
    :lat_arr: (float) array of Latitude
    :search_range: (float) range of longitude and latitude for faster search

    """
    lon_mask = (lon_arr>=loni-search_range) & (lon_arr<=loni+search_range)
    lat_mask = (lat_arr>=lati-search_range) & (lat_arr<=lati+search_range)
    lonlat_mask = lon_mask & lat_mask
    lon_sub = lon_arr[lonlat_mask]
    lat_sub = lat_arr[lonlat_mask]
    pts = np.array([loni,lati])
    tree = spatial.KDTree(list(zip(lon_sub, lat_sub)))
    p = tree.query(pts)
    cidx = p[1]
    idx = np.argwhere(lon_arr==lon_sub[cidx])
    for i in idx[0][:]:
        if lat_arr[i] == lat_sub[cidx]:
            out = i
            break
    return out

def gc_radius():
    """Return the radius of Earth
    :returns: (float) radius of Earth in km

    """
    return 6371.0

def gc_angle(lon0, lat0, lon1, lat1):
    """Calculate the angle counterclockwise from east.
    :lon0: (float) longitude of point 1 in degrees
    :lat0: (float) latitude of point 1 in degrees
    :lon1: (float) longitude of point 2 in degrees
    :lat1: (float) latitude of point 2 in degrees
    :returns: (float) angle in degrees
    """
    dlon_r = np.radians(lon1-lon0)
    dlat_r = np.radians(lat1-lat0)
    angle = np.arctan2(dlat_r, dlon_r)
    return angle

def gc_angles(lon, lat):
    """A wrapper of gc_angle to compute the angle counterclockwise from east for an array of lon and lat
    :lon: (numpy array) array of longitudes
    :lat: (numpy array) array of latitudes

    """
    lat0 = np.zeros(lat.size)
    lon0 = np.zeros(lon.size)
    lat1 = np.zeros(lat.size)
    lon1 = np.zeros(lon.size)
    lat0[1:-1] = lat[0:-2]
    lat1[1:-1] = lat[2:]
    lon0[1:-1] = lon[0:-2]
    lon1[1:-1] = lon[2:]
    angles = gc_angle(lon0, lat0, lon1, lat1)
    angles[0] = angles[1]
    angles[-1] = angles[-2]
    return angles

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

def plot_basemap(region='Global', axis=None, projection=None):
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
    else:
        if projection is None:
            projection = 'cass'
        # regional map
        region_obj = region_latlon(region)
        lon_ll, lat_ll, lon_ur, lat_ur = region_obj.lon_ll, region_obj.lat_ll, region_obj.lon_ur, region_obj.lat_ur
        lon_c = 0.5*(lon_ll+lon_ur)
        lat_c = 0.5*(lat_ll+lat_ur)
        m = Basemap(projection=projection, llcrnrlon=lon_ll, llcrnrlat=lat_ll,
                urcrnrlon=lon_ur, urcrnrlat=lat_ur, resolution='l', lon_0=lon_c, lat_0=lat_c, ax=axis)
        # parallels and meridians
        mdlat = 10.0
        mdlon = 10.0
    # plot coastlines, draw label meridians and parallels.
    m.drawcoastlines(zorder=3)
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='gray',lake_color='lightgray', zorder=2)
    m.drawparallels(np.arange(-90.,91.,mdlat), labels=[1,0,0,1])
    m.drawmeridians(np.arange(-180.,181.,mdlon), labels=[1,0,0,1])
    return m

def plot_transect_normal(mpaso_data_x, mpaso_data_y, transect, name='Normal Component', **kwargs):
    # cross section of zonal and meridional component
    mpaso_vcsec_x = mpaso_data_x.get_transect(transect)
    mpaso_vcsec_y = mpaso_data_y.get_transect(transect)
    # show locations of data point along cross section
    lon_cs = mpaso_vcsec_x.lon
    lat_cs = mpaso_vcsec_x.lat
    # transect of normal component
    angles_cs = gc_angles(lon_cs, lat_cs)
    depth_cs = mpaso_vcsec_x.depth
    nd = depth_cs.size
    data_cs = np.zeros(mpaso_vcsec_x.data.shape)
    for i in np.arange(nd):
        data_cs[:,i] = -mpaso_vcsec_x.data[:,i]*np.sin(angles_cs)+mpaso_vcsec_y.data[:,i]*np.cos(angles_cs)
    mpaso_vcsec = MPASOVertCrossSection(data=data_cs, lon=lon_cs, lat=lat_cs, dist=mpaso_vcsec_x.dist,
                                         depth=mpaso_vcsec_x.depth, name=name,
                                         units=mpaso_vcsec_x.units)
    fig = mpaso_vcsec.plot(**kwargs)
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
