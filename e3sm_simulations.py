import os
import xml.etree.ElementTree as et

#--------------------------------
# E3SM Simulation
#--------------------------------
class E3SMSimulation(object):

    """E3SM simulation object"""

    def __init__(self, database=None, runname=None):
        """Initialize E3SMSimulation

        :database: (str) full file name of simulation database
        :runname: (str) short name of run

        """
        self.database = database
        self.runname = runname
        # read database
        tree = et.parse(database)
        root = tree.getroot()
        for child in root.findall('simulation'):
            name = child.attrib['name']
            if runname == name:
                for grandchild in child.getchildren():
                    setattr(self, grandchild.tag, grandchild.text)

    def get_path(self, comp='ocn', climo_yr0=None, climo_yr1=None, ts_yr0=None, ts_yr1=None, rest_yr=None):
        """Get paths the simulation

        :comp: (str) component
        :climo_yr0: (int) starting year of climatology data
        :climo_yr1: (int) ending year of climatology data
        :ts_yr0: (int) starting year of timeseries data
        :ts_yr1: (int) ending year of timeseries data
        :rest_yr: (int) year of restart data

        """
        # check machine
        hostname = os.uname()[1]
        assert self.machine in hostname, "Data hosted on {}, running on {}. Stop.".format(self.machine, hostname)
        # check arguments
        assert climo_yr0 is not None, "Input argument climo_yr0 is None. Stop."
        assert ts_yr0 is not None, "Input argument ts_yr0 is None. Stop."
        if climo_yr1 is None:
            climo_yr1 = climo_yr0+1
        if ts_yr1 is None:
            ts_yr1 = ts_yr0+1
        if rest_yr is None:
            rest_yr = ts_yr1+1
        # set paths
        data_root = self.dataroot+'/'+self.longname
        climo_root = self.climoroot+'/e3sm_climo/'+self.longname+'/{:04d}-{:04d}/'.format(climo_yr0, climo_yr1)+comp
        ts_root = self.climoroot+'/e3sm_ts/'+self.longname+'/{:04d}-{:04d}/'.format(ts_yr0, ts_yr1)+comp
        fig_root = os.environ['HOME']+'/work/e3sm_res_cmp/figures/'+self.runname+'/{:04d}-{:04d}'.format(climo_yr0, climo_yr1)
        if self.roottype == 'run':
            rst_root = data_root+'/run'
            mon_root = data_root+'/run'
        elif self.roottype == 'archive':
            rst_root = data_root+'/rest/{:04d}-01-01-00000'.format(rest_yr)
            mon_root = data_root+'/'+comp+'/hist'
        elif self.roottype == 'test':
            rst_root = data_root
            mon_root = data_root
        else:
            raise ValueError('Root type \'{}\' not supported'.format(self.roottype))
        path = {'rst_root': rst_root,
                'mon_root': mon_root,
                'climo_root': climo_root,
                'ts_root': ts_root,
                'fig_root': fig_root}
        return path


