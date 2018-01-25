# -*- coding: utf-8 -*-
"""Plotting Library for wrfout files.

author: Lukas Kugler

Major bug:  - PRIORITY LOW
kommt nicht zurecht, wenn nur ein wrfout file im ordner liegt!

"""
from __future__ import division
import os, sys
import numpy as np
from datetime import datetime

from netCDF4 import Dataset
import wrf
from wrf import (getvar, ALL_TIMES, interplevel, to_np, latlon_coords,
                 disable_xarray, enable_xarray, smooth2d, get_basemap,
                 CoordPair, vertcross)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('font',**{'family':'monospace'})
matplotlib.rcParams.update({'font.size': 9})

from mypickle import load_pickle, save_pickle
from aux import Timer, OOMFormatter, truncate_colormap

modulepath = os.path.dirname(__file__)

kts_per_mps = 1.94384   # kts = mps * kts_per_mps
c_p = 1005
R_d = 287.16

pix_width = 1080
pix_height = 810 #756

def _check_2D(data):
    if len(data.shape) == 3:
        N_images = data.shape[0]
    elif len(data.shape) == 2:
        N_images = 1
    else:
        raise NotImplementedError(['Data more than 2D! data.shape:', data.shape])
    return N_images

def PressReduction(QFE, ALT):
    """Pressure reduction according to ICAO.

    Args:
        QFE ... PSFC in hPa, shape (2 dimensional)
        ALT ... altitude (HGT) of surface in meters
    """
    if len(QFE.shape)>2:
        print 'input dimension too big; QFE.shape:', QFE.shape

    PMSL = np.empty(QFE.shape)
    gamma=  -0.0065
    P0=     1013.25
    g=      9.80665
    R_air = 287.16
    #PMSL[m,n]=  P_SFC[m,n]*np.exp(-g*HGT[m,n]/R_air/(T_SFC[m,n]-HGT[m,n]*g))
    PMSL= QFE*(1-gamma*ALT/(288.15*(QFE/P0)**(-R_air*gamma/g)))**(-g/R_air/gamma)
    return PMSL

class Load_Domain(object):
    """HIGH LEVEL THING

    import EasyWRF
    ThisWRF = EasyWRF.Load_Domain('./wrfout_dir/', 'DomainName')

    ThisWRF.plot('HGT')      # variables available within wrfout file
    ThisWRF.plot_derived('BT')  # calculated Brightness Temperature from OLR

    """
    def __init__(self, WRFout_data_directory, domainname="MyDomain",
                 coast=True,
                 river=False,
                 austrianborders=True,
                 bezirk=True
                 imgw=False):
        self.coast = coast
        self.river = river
        self.austrianborders = austrianborders
        self.bezirk = bezirk

        self.load_WRFout_data(WRFout_data_directory)
        self.domainname = domainname

        first_wrfout = self.data[0]
        self.init = first_wrfout.SIMULATION_START_DATE
        self.dx = first_wrfout.DX

        hgt = getvar(self.data, "HGT")
        self.hgt = hgt
        lats, lons = latlon_coords(hgt)

        try:
            fname = os.getcwd()+'/basemap_'+domainname
            self.bm = load_pickle(os.getcwd()+'/basemap_'+domainname)
            print 'precomputed basemap loaded.'

            int2 = lambda x: np.floor(x/100)*100
            if not (int2(np.nanmax(lats)) == int2(self.bm.latmax) and
                    int2(np.nanmin(lats)) == int2(self.bm.latmin)):
                print np.nanmax(lats), self.bm.latmax
                print np.nanmin(lats), self.bm.latmin
                raise RuntimeError('Saved basemap does not match up with wrfout files. May be resolved by using a new `domainname`')
        except IOError:
            print 'setting up basemap ... this may take up to two minutes'

            # basemap resolution
            if self.dx < 8000:
                bm_res = 'h'
            else:
                bm_res = 'i'

            with Timer('basemap loading time'):
                self.bm = get_basemap(hgt, resolution=bm_res)
            fname = os.getcwd()+'/basemap_'+domainname
            save_pickle(self.bm, fname=fname)
            print 'basemap saved to', fname

        self.x, self.y = self.bm(to_np(lons), to_np(lats))
        self.XY_gridpoints = hgt.size
        #print 'XY_gridpoints:', self.XY_gridpoints

        disable_xarray()
        dates = getvar(self.data, "times", timeidx=ALL_TIMES, method="cat")
        self.wrfdate = [datetime.strptime(str(each), "%Y-%m-%dT%H:%M:%S.000000000") for each in dates]
        enable_xarray()

    def load_WRFout_data(self, WRFout_data_directory):
        """Load wrfout files from given directory into list of Datasets."""
        wrfpath = WRFout_data_directory

        # last str element must be /
        if wrfpath[-1] is not '/':
            wrfpath+='/'

        wrffiles = []
        for each in os.listdir(wrfpath):
            if each.startswith('wrfout_d'):
                wrffiles.append(each)

        if len(wrffiles) < 2:
            print 'wrffiles:',wrffiles
            raise IOError('Not enough wrfout files.')

        wrflist = [wrfpath+each for each in wrffiles]
        wrflist.sort()

        datalist = []
        for i in range(len(wrflist)):
            datalist.append(Dataset(wrflist[i]))

        self.data = datalist
        print 'successfully loaded', len(datalist), 'wrfout files'

        # TODO: python indexing for 1 file?
        if len(datalist)==1:
            print 'WARN: for 1 wrfout file other python indexing neccessary -> will fail!'

    def _make_layout(self, figtitle, valid_datetime):
        """Prepare the figure by plotting titles, country boundaries etc."""
        self.figtitle = figtitle

        initdate = datetime.strptime(self.init, '%Y-%m-%d_%H:%M:%S')
        initdate = initdate.strftime('Init: %a, %d.%m.%Y %H:%M UTC')

        # Create the figure
        plt.ioff()
        fig, ax = plt.subplots(figsize=(pix_width/100, pix_height/100), dpi=100)

        if imgw:  # Plot UNIVIE Logo
            try:
                img = plt.imread(modulepath+'/../shapefiles/logo_imgw.png')
                fig.figimage(img, 0, 0)
            except Exception as e:
                print e, 'unable to plot imgw logo'

        # Plot grid lines dependent on resolution
        if self.dx < 12000:
            self.bm.drawparallels(np.arange(10, 60, 2.5), labels=[True, False, False, False])  # label left right top bottom
            self.bm.drawmeridians(np.arange(0, 50, 2.5), labels=[False, False, False, True])
        else:
            self.bm.drawparallels(np.arange(15, 70, 10), labels=[True,False,False,False])  # label left right top bottom
            self.bm.drawmeridians(np.arange(-40, 50, 10), labels=[False, False, False, True])

        if self.austrianborders:
            try:
                #self.bm.readshapefile(modulepath+'/../shapefiles/AUT_adm0', 'AUT_adm0', color='k', linewidth=.6)
                self.bm.readshapefile(modulepath+'/../shapefiles/AUT_adm1', 'AUT_adm1', color='k', linewidth=.4)
                if self.bezirk:
                    self.bm.readshapefile(modulepath+'/../shapefiles/AUT_adm2', 'AUT_adm2', color='k', linewidth=.1)
            except:
                print "unable to plot gi shape files"
                raise
        if self.river:
            self.bm.drawrivers(linewidth=0.6)
        if self.coast:
            self.bm.drawcoastlines(linewidth=0.6)
        self.bm.drawcountries(linewidth=0.6)

        # Adding 'initdate' and 'validdate' as textboxes
        validdate = valid_datetime.strftime('Valid: %a, %d.%m.%Y %H:%M UTC')
        wrf_version = 'WRF-V3.8'  # TODO: get from file
        resolution = (wrf_version+' @ ' + str(np.round(self.dx/1000,1)) + 'km')
        copyright = '$\copyright$ Lukas Kugler'

        ax.text(1, 1.04, initdate,
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=9)
        ax.text(1, 1.01, validdate,
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=9)
        ax.text(1.008, 0.85, resolution,
                verticalalignment='center', horizontalalignment='left',
                rotation=90,
                transform=ax.transAxes,
                color='black', fontsize=8)
        ax.text(0, 1.01, figtitle,
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes,
                fontsize=11)
        ax.text(1.008, 0.15, copyright,
                verticalalignment='center', horizontalalignment='left',
                rotation=90,
                transform=ax.transAxes,
                color='black', fontsize=8)

        fig.subplots_adjust(top=0.94, bottom=0.11)
        self.fig = fig
        self.ax = ax

        self.color_min = None
        self.color_max = None
        self.cbar_tickformat = None

    def _finalize_save(self, var_name, valid_datetime):
        """Saves the current figure."""
        init = datetime.strptime(self.init, '%Y-%m-%d_%H:%M:%S')
        savedir = './img/'+init.strftime('%Y%m%d%H')+'/'+self.domainname+'/'
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        timestamp = valid_datetime.strftime('%y%m%d_%H%M')
        filename = var_name + '_' + timestamp + '_' + self.domainname + '.png'
        fname = savedir + filename
        plt.savefig(fname, dpi=100)
        plt.close('all')
        print fname, 'saved'

    #####################################################################################
    def _load_var(self, var):
        return to_np(getvar(self.data, var, timeidx=ALL_TIMES, method="cat"))

    """Plotting Routines"""
    def _plot_heatmap(self, data_2D,
                        interpolation = 'nearest',
                        cmap='viridis',
                        color_min=None, color_max=None,
                        N_colors=14,
                        N_ticks = None,
                        cbar_tickformat=None,
                        cbar_label=None,
                        cbar_ticks=None):

        if N_ticks is None:
            N_ticks = N_colors+1

        if color_min is not None and color_max is not None and cbar_ticks is None:
            cbar_ticks = np.linspace(color_min, color_max, N_ticks, endpoint=True)

        """Plot the Heatmap"""
        norm = matplotlib.colors.BoundaryNorm(cbar_ticks, 256)
        self.heatmap = self.bm.imshow(data_2D, interpolation=interpolation,
                                      vmin=color_min, vmax=color_max,
                                      cmap=cmap, norm=norm, zorder=1)

        """Add the colorbar"""
        clb_psn = self.fig.add_axes([0.3,0.06,0.58,0.02])  # left, bottom, width, height
        if cbar_ticks is not None:
            clb = plt.colorbar(self.heatmap, orientation="horizontal",
                               cax=clb_psn,
                               ticks=cbar_ticks,
                               boundaries=cbar_ticks,
                               format=cbar_tickformat)
            clb.outline.set_linewidth(.4)
            plt.clim(cbar_ticks[0], cbar_ticks[-1])
            if cbar_label is not None:
                clb.set_label(cbar_label)

    ###################################################################################
    """LEGACY CODE"""

    # def _get_div10(self):
    #     if not hasattr(self, 'u10') or not hasattr(self, 'v10'):
    #         self._get_u10()
    #     if len(self.u1.shape)==3:
    #         dv = self.v1[:,1:,:]-self.v1[:,:-1,:]   # 2nd : North/south, 3rd: west-east!
    #         du = self.u1[:,:,1:]-self.u1[:,:,:-1]
    #         self.div1 = du/self.dx + dv/self.dx
    #
    #         self.divmax = np.nanmax(self.div1)
    #         self.divmin = np.nanmin(self.div1)
    #     else:
    #         raise NotImplementedError('len(self.u10.shape):',len(self.u10.shape))

    def _get_u10(self):
        print 'loading u1, v1 ...'
        self.u1 = to_np(getvar(self.data, "U", timeidx=ALL_TIMES, method="cat"))[:,0,:,:]                # 10m u-component of wind [ms^-1]
        self.v1 = to_np(getvar(self.data, "V", timeidx=ALL_TIMES, method="cat"))[:,0,:,:]
        self.u10 = to_np(getvar(self.data, "U10", timeidx=ALL_TIMES, method="cat"))              # 10m u-component of wind [ms^-1]
        self.v10 = to_np(getvar(self.data, "V10", timeidx=ALL_TIMES, method="cat"))

    ####################

    def _makeplot_sfcflux(self, i, data_2D, cmap='RdBu_r'):
        if not hasattr(self, 'u10') or not hasattr(self, 'v10'):
            self._get_u10()

        objective_resolution = 80 # wind barbs
        mean_xpoints = self.XY_gridpoints**(.5)
        nbarbs = 1 + int(mean_xpoints/objective_resolution)

        self.bm.barbs(self.x[::nbarbs,::nbarbs], self.y[::nbarbs,::nbarbs],
                 self.u10[i, ::nbarbs, ::nbarbs]*kts_per_mps,
                 self.v10[i, ::nbarbs, ::nbarbs]*kts_per_mps,
                 zorder=2, length=5, linewidth=0.5)

        # Make the contour outlines and filled contours for the smoothed sea level pressure.
        #bm.contour(self.x[1:,1:], self.y[1:,1:], div10, 7, colors="green", linewidths=0.5)

        color_min = -600
        color_max = 600
        N_colors = 12
        self.cbar_ticks = np.linspace(color_min, color_max, N_colors+1, endpoint=True)
        self.cbar_tickformat = "%d"  # OOMFormatter(-3, mathText=False)  #1e-3


        norm = matplotlib.colors.BoundaryNorm(self.cbar_ticks, 256)
        self.cf1 = self.bm.contourf(self.x, self.y, data_2D,
                                    self.cbar_ticks, vmin=color_min, vmax=color_max,
                                    cmap=cmap, norm=norm, zorder=1)
        self.color_min = color_min
        self.color_max = color_max

    def _plot_syn_boden(self, PMSL, U, V, T):
        color_min = -20
        color_max = 49
        N_colors =  23
        self.cbar_ticks = np.linspace(color_min, color_max, N_colors+1, endpoint=True)
        self.cbar_tickformat = "%d"

        norm = matplotlib.colors.BoundaryNorm(self.cbar_ticks, 256)
        self.cf1 = self.bm.contourf(self.x, self.y, T-273.15+290,
                                    self.cbar_ticks, vmin=color_min, vmax=color_max, alpha=0.7,
                                    cmap='gist_ncar', norm=norm, zorder=1)
        self.color_min = color_min
        self.color_max = color_max

        minP = np.nanmin(PMSL)
        maxP = np.nanmax(PMSL)

        linelist1 = np.arange(int(minP/5)*5, int(maxP/5)*5+1, 5)
        linelist2 = np.arange(int(minP/5)*5, int(maxP/5)*5+1, 1)

        con = self.bm.contour(self.x, self.y, PMSL, linelist1, colors="k", linewidths=2, zorder=3)
        self.bm.contour(self.x, self.y, PMSL, linelist2, colors="k", linewidths=1, zorder=2)
        plt.clabel(con, inline=1, fontsize=8, fmt='%d')

        objective_resolution = 24 # wind barbs
        mean_xpoints = self.XY_gridpoints**(.5)
        nbarbs = 1 + int(mean_xpoints/objective_resolution)

        self.bm.barbs(self.x[::nbarbs,::nbarbs], self.y[::nbarbs,::nbarbs],
                         U[::nbarbs, ::nbarbs]*kts_per_mps,
                         V[::nbarbs, ::nbarbs]*kts_per_mps,
                         zorder=2, length=6, linewidth=0.5)

    def _add_barbs(self, U, V):
        objective_resolution = 30 # wind barbs
        mean_xpoints = self.XY_gridpoints**(.5)
        nbarbs = 1 + int(mean_xpoints/objective_resolution)

        print U.shape, V.shape, self.x.shape, self.y.shape

        self.bm.barbs(self.x[::nbarbs,::nbarbs], self.y[::nbarbs,::nbarbs],
                        U[::nbarbs, ::nbarbs],
                        V[::nbarbs, ::nbarbs],
                        zorder=2, length=6, linewidth=0.5)

    def _plot_synoptik1(self, H500, PMSL, RELTOP):

        color_min = np.nanmin(RELTOP)
        color_max = np.nanmax(RELTOP)
        N_colors =  14
        self.cbar_ticks = np.linspace(color_min, color_max, N_colors+1, endpoint=True)
        self.cbar_tickformat = "%d"

        norm = matplotlib.colors.BoundaryNorm(self.cbar_ticks, 256)
        self.cf1 = self.bm.contourf(self.x, self.y, RELTOP,
                                    self.cbar_ticks, vmin=color_min, vmax=color_max, alpha=0.7,
                                    cmap='gist_ncar', norm=norm, zorder=1)

        linelist1 = np.arange(int(color_min/5)*50, int(color_max/5)*50+1, 500)
        con = self.bm.contour(self.x, self.y, H500, linelist1, colors="black", linewidths=1.2, zorder=2)

        minP = np.nanmin(PMSL)
        maxP = np.nanmax(PMSL)

        linelist1 = np.arange(int(minP/5)*5, int(maxP/5)*5+1, 5)
        con = self.bm.contour(self.x, self.y, PMSL, linelist1, colors="black", linewidths=0.8, zorder=3)
        plt.clabel(con, inline=1, fontsize=8, fmt='%d')

    def _plot_vorticity(self, Vort, minval, color_max):
        # round to order of magnitude
        self.color_min = minval
        self.color_max = color_max

        N_colors = 20

        self.cbar_ticks = np.linspace(self.color_min, self.color_max, N_colors+1, endpoint=True)
        self.cbar_tickformat = "%d"

        norm = matplotlib.colors.BoundaryNorm(self.cbar_ticks, 256)
        self.cf1 = self.bm.imshow(Vort, interpolation='nearest',
                                  vmin=self.color_min, vmax=self.color_max, alpha=0.7,
                                  cmap=plt.get_cmap("RdBu_r"), norm=norm, zorder=1)

    def _how2plot(self, wrfout_variable, data):
        """How to plot a wrfout variable"""
        # standard values
        cmap = 'viridis'
        color_min = np.nanmin(data)
        color_max = np.nanmax(data)
        print 'data min: '+str(color_min)+' max: '+str(color_max)
        N_colors = 10

        cbar_tickformat = "%d" #OOMFormatter(-3, mathText=False)  #1e-3
        cbar_label = None
        only_first_image = False

        if wrfout_variable == 'LU_INDEX':
            cmap = 'gist_ncar_r'
            N_colors = (color_max-color_min)
            cbar_tickformat = "%d"
            cbar_label = ''
            only_first_image = True

        if wrfout_variable == 'HGT':
            color_min = 0
            if color_max-color_min <= 800:
                dH = 50
                color_max = int(np.nanmax(data)/dH+1)*dH
                N_colors = (color_max-color_min)/dH  # 250 m
            elif color_max-color_min <= 1900:
                dH = 100
                color_max = int(np.nanmax(data)/dH+1)*dH
                N_colors = (color_max-color_min)/dH  # 250 m
            else:
                dH = 250
                color_max = int(np.nanmax(data)/dH+1)*dH
                N_colors = (color_max-color_min)/dH  # 250 m

            cbar_tickformat = "%d"
            cbar_label = 'm'
            only_first_image = True

            cmap = plt.get_cmap('terrain')
            cmap = truncate_colormap(cmap, 0, 1.225)

        if wrfout_variable == 'RAINC' or wrfout_variable == 'RAINNC':
            color_min = 0
            color_max = 100
            N_colors = 20
            cbar_tickformat = "%d"
            cbar_label = 'mm'

            cmap = 'gist_ncar_r'

        if wrfout_variable == 'OLR':
            cmap = 'gist_ncar_r'
            color_min = 100
            color_max = 320
            N_colors = 22

            cbar_tickformat = "%d" #OOMFormatter(-3, mathText=False)  #1e-3
            cbar_label = 'W/m$^2$'

        if wrfout_variable == 'SMOIS':
            cmap = 'RdBu'
            color_min = 0
            color_max = 1
            N_colors = 10

            cbar_tickformat = "%1.1f" #OOMFormatter(-3, mathText=False)  #1e-3
            cbar_label = 'm$^3$/m$^3$'

        if wrfout_variable == 'SH2O':
            cmap = 'RdBu'
            color_min = 0
            color_max = 1
            N_colors = 10

            cbar_tickformat = "%1.1f"
            cbar_label = 'm$^3$/m$^3$'

        return cmap, color_min, color_max, N_colors, cbar_tickformat, cbar_label, only_first_image

    ##################
    # EASY TO USE FUNCTIONs

    def plot(self, wrfout_variable):
        """ Simple to use routine for plotting builtin wrfout variable

        ### Example:
        wrfpath = './'
        ThisWRF = wrfout.Load_Domain(wrfpath, domainname='CEUR')

        ThisWRF.plot2Dvar('OLR')


        programming a routine shall be simple!
        the string of the wrfout-variable will be given as Input
        need to know:
            - type of plot (2D heatmap, ...)
            - minimum / maximum value of colorbar, otherwise will be guessed
            - number of colors
            - tickformat for colorbar
            - colorbar label

        """

        try:
            print 'loading '+str(wrfout_variable)+'...'
            data = (getvar(self.data, wrfout_variable, timeidx=ALL_TIMES, method="cat"))             # 10m u-component of wind [ms^-1]
        except:
            print 'wrfout_variable is', wrfout_variable
            print 'available are:', self.data[0].variables.keys()
            raise

        N_images = data.shape[0]
        if len(data.shape) == 4:
            data = data[:,0,:,:]
        elif len(data.shape) > 4 or len(data.shape) < 3:
            raise NotImplementedError(['Field is more than 2D! data.shape:', data.shape])

        """How to plot a variable in wrfout"""
        cmap, color_min, color_max, N_colors, cbar_tickformat, cbar_label, only_first_image = self._how2plot(wrfout_variable, data)

        if only_first_image:
            N_images = 1
            plotdata = data[0]

        """Plotting loop"""
        first_wrfout = self.data[0]
        figtitle = first_wrfout[wrfout_variable].description
        save_name = wrfout_variable

        for i in range(N_images):
            valid_datetime = self.wrfdate[i]
            self._make_layout(figtitle, valid_datetime)

            # kA
            if len(self.data) == 1:
                plotdata = data
            else:
                plotdata = data[i]

            self._plot_heatmap(plotdata,
                              cmap = cmap,
                              color_min = color_min,
                              color_max = color_max,
                              N_colors = N_colors,
                              cbar_tickformat = cbar_tickformat,
                              cbar_label = cbar_label)

            self._finalize_save(save_name, valid_datetime)


    def plot_derived(self, var='TE'):
        """ Simple to use routine for plotting derived variable
        from available list: 'TE' (Total Energy),

        ### Example:
        wrfpath = './'
        ThisWRF = wrfout.Load_Domain(wrfpath, domainname='CEUR')

        ThisWRF.plot2D_derived(var='TE')
        """
        avail = {'TE':'Total Energy',
                 'Synoptik': 'H500, MSLP, Relative Topography H500-H1000',
                 'SYN-SFC': 'MSLP + Wind & Theta @ lowest model level',
                 'SRH': '0-3 km Storm Relative Helicity',
                 'max_dbz': 'Colucolor_min Maximum Simulated Reflectivity',
                 'Column_SMOIS': 'Column Mean Soil Moisture +H2O',
                 'BT': 'Brightness Temperature'
                 }

        if var=='TE':
            try:
                print 'loading variables ...'
                U = _load('U')[:,:-1,:,:]      #  [:,:-1,:,:] ghört weg, is nur weil vertical grid dp unklar is
                V = _load('V')[:,:-1,:,:]
                P = _load('P')
                Ps = _load('PSFC')
                T = _load('T')[:,:-1,:,:]
                T00 = [float(self.data[i]['T00'][:]) for i in range(len(self.data))]
                for i in range(len(self.data)):
                    T[i,:,:,:] = T[i,:,:,:] + T00[i]
            except IOError as e:
                print(e, 'loading failed')

            # calculate derived variables
            # 1dim:time, 2dim:press, 3dim:northsouth, 4dim:westeast
            U = (U[:,:,:,1:] + U[:,:,:,:-1])/2   # arakawa C grid
            V = (V[:,:,1:,:] + V[:,:,:-1,:])/2   # arakawa C grid
            dp = P[:,:-1:,:,:] - P[:,1:,:,:]     # pressure diff btw sigma levels
            dA = self.dx**2                  # Area, dx*dy constant (not sure)

            Tr = 280    # reference Temperature; not 100% sure
            pr = 1e5   # reference pressure

            part1 = (U**2+V**2+c_p/Tr*T**2)*dp*dA
            part2 = ((np.log(Ps))**2)*dA
            TE = .5*np.sum(part1, axis=1) + .5*R_d*Tr*pr*part2   # for a colucolor_min of air

            dTE = TE[:,:,:] - TE[1,:,:]

            # Error check for plotting
            N_images = _check_2D(dTE)

            # plotting range
            minval = np.nanmin(dTE[1:,:,:])
            color_max = np.nanmax(dTE[1:,:,:])
            print 'minval, color_max', minval, color_max

            # plotting options
            figtitle = 'Total Energy increase since 1st image'
            save_name = 'TE'

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(dTE[i], N_colors=14)

                self._finalize_save(save_name, valid_datetime)
                plt.close('all')

        elif var=='SYN-SFC':
            try:
                print 'loading variables ...'
                U = _load('U')[:,0,:,:]
                V = _load('V')[:,0,:,:]
                T = _load('T')[:,0,:,:]
                PSFC = _load('PSFC')
            except IOError as e:
                print(e, 'loading failed')

            # umrechnen in PMSL
            print 'pressure reduction ...'
            # PMSL = PSFC/100
            # for i in range(PSFC.shape[0]):
            #     PMSL[i,:,:] = PressReduction(PMSL[i,:,:], self.hgt)
            #     PMSL[i,:,:] = smooth2d(PMSL[i,:,:], 4)   # glaetten

            PMSL = getvar(self.data, "slp", timeidx=ALL_TIMES)
            PMSL = smooth2d(PMSL, 6)

            # plotting options
            print 'plotting ...'
            figtitle = avail[var]
            save_name = 'SYN-SFC'

            N_images = _check_2D(PMSL)  # Error check for plotting

            if np.nanmax(PMSL)>1200 or np.nanmin(PMSL)<700:
                raise RuntimeError('PMSL not in hPa; max(p):'+str(np.nanmax(PMSL))+'min(p)'+str(np.nanmin(PMSL)))

            self.cbar_label = r'$^\circ$C'
            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_syn_boden(PMSL[i], U[i,:,:], V[i,:,:], T[i,:,:])

                self._finalize_save(save_name, valid_datetime)
                plt.close('all')

        elif var=='Synoptik':
            try:
                print 'loading variables ...'
                U = _load('U')[:,:-1,:,:]      #  [:,:-1,:,:] ghört weg, is nur weil vertical grid dp unklar is
                V = _load('V')[:,:-1,:,:]
                #P = _load('P')
                PSFC = _load('PSFC')
                T = _load('T')[:,:-1,:,:]
                T00 = [float(self.data[i]['T00'][:]) for i in range(len(self.data))]
                for i in range(len(self.data)):
                    T[i,:,:,:] = T[i,:,:,:] + T00[i]
            except IOError as e:
                print(e, 'loading failed')

            # umrechnen in PMSL
            print 'pressure reduction ...'
            PMSL = PSFC/100
            for i in range(PSFC.shape[0]):
                PMSL[i,:,:] = PressReduction(PMSL[i,:,:], self.hgt)
                PMSL[i,:,:] = smooth2d(PMSL[i,:,:], 5)   # glaetten

            # plotting options
            figtitle = avail[var]
            save_name = 'SYN'

            N_images = _check_2D(PMSL)  # Error check for plotting

            for i in range(N_images):

                P = getvar(self.data[i], "pressure")
                Z = getvar(self.data[i], "z", units="dm")
                H500 = interplevel(Z, P, 500)
                H850 = interplevel(Z, P, 850)
                H1000 = interplevel(Z, P, 1000)
                Reltop = H850-H1000

                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_synoptik1(to_np(H500), PMSL[i,:,:], to_np(Reltop))

                self._finalize_save(save_name, valid_datetime)
                plt.close('all')

        elif var=='SRH':
            try:
                print 'loading variables ...'
                U = _load('U')[:,0,:,:]
                V = _load('V')[:,0,:,:]
            except IOError as e:
                print(e, 'loading failed')

            SRH = wrf.helicity.get_srh(self.data, timeidx=wrf.ALL_TIMES)

            # plotting range
            minval = np.floor(np.nanmin(SRH)/100)*100
            color_max = np.ceil(np.nanmax(SRH)/100)*100

            minval = -250 #-max(abs(minval),abs(color_max))
            color_max = 250 #max(abs(minval),abs(color_max))
            print 'minval, color_max', minval, color_max

            # plotting options
            figtitle = avail[var]
            save_name = 'SRH'
            self.cbar_label = 'J/kg'

            N_images = _check_2D(SRH)  # Error check for plotting

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_vorticity(SRH[i], minval, color_max)

                self._finalize_save(save_name, valid_datetime)
                plt.close('all')

        elif var=='max_dbz':

            data2D = getvar(self.data, 'mdbz', timeidx=wrf.ALL_TIMES)
            #data2D = np.amax(data2D, axis=-3)
            # from wrf import uvmet
            # U, V = uvmet.get_uvmet(self.data, timeidx=wrf.ALL_TIMES)[:,0,:,:]
            # print U, V
            # plotting range
            data2D = np.array(data2D)
            data2D[data2D < 0.2] = np.nan

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            minval = np.log10(15)
            maxval = np.log10(140)
            cbar_ticks = np.logspace(minval, maxval, 12, endpoint=True)

            # plotting options
            figtitle = 'Max. Column Reflectivity'
            save_name = 'max_dbz'
            cbar_label = 'dBz'

            cmap = plt.get_cmap('gist_ncar')
            new_cmap = truncate_colormap(cmap, 0.2, 0.9)

            N_images = _check_2D(data2D)  # Error check for plotting

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(data2D[i], cmap=new_cmap,
                                interpolation = 'bilinear',
                                cbar_tickformat = "%0.1f",  #1e-3,
                                cbar_label = cbar_label,
                                cbar_ticks = cbar_ticks)
                #self._add_barbs(U[i], V[i])  # add wind lowest model level

                self._finalize_save(save_name, valid_datetime)

        elif var == 'Column_SMOIS':
            smois = _load('SMOIS')
            sh2o = _load('SH2O')
            # thickness_soil = self.data[0]['DZS'][:]  # thickness of soil layer
            data2D = np.mean(smois, axis=-3) + np.mean(sh2o, axis=-3)

            color_min = np.nanmin(data)
            color_max = np.nanmax(data)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = avail[var]
            save_name = 'Column_SMOIS'

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(data2D[i], cmap='RdBu',
                                 color_min=0, color_max=1,
                                 N_colors=10,
                                 cbar_tickformat="%1.1f",
                                 cbar_label="m$^3$/m$^3$")
                #self._add_barbs(U[i], V[i])  # add wind lowest model level

                self._finalize_save(save_name, valid_datetime)

        elif var == 'Column_QVAPOR':
            Q = _load('QVAPOR')
            cloud = _load('QCLOUD')
            rain = _load('QRAIN')+ _load('QICE')

            data2D = np.nansum(Q, axis=-3)
            cloud = np.nansum(cloud, axis=-3)
            rain = np.nansum(rain, axis=-3)

            cloud_contours = (1e-3*np.max(cloud), np.max(cloud))
            rain_contours = (1e-3*np.max(rain), np.max(rain))

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = 'Column Vapor'
            save_name = 'QVAPOR'

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                # self.bm.contourf(mask, interpolation='nearest',
                #                      cmap=plt.get_cmap('gray'), alpha=0.9, zorder=3)

                self.bm.contourf(self.x, self.y, cloud[i],
                                 cloud_contours, colors='none',
                                 hatches=['o'], zorder=5)

                self.bm.contourf(self.x, self.y, rain[i],
                                 rain_contours, colors='none',
                                 hatches=["\ "], zorder=6)

                self._plot_heatmap(data2D[i], cmap='gist_gray_r',
                                 interpolation = 'bilinear',
                                 color_min=0, color_max=0.26,
                                 N_colors=26, N_ticks=14,
                                 cbar_tickformat=OOMFormatter(-3, fformat="%d", mathText=False),  #1e-3,
                                 cbar_label="propto kg/m^2")
                #self._add_barbs(U[i], V[i])  # add wind lowest model level

                self._finalize_save(save_name, valid_datetime)

        elif var == 'Column_Q-Hydro':
            QC = _load('QCLOUD') + _load('QRAIN')+ _load('QICE')
            # thickness_soil = self.data[0]['DZS'][:]  # thickness of soil layer
            data2D = np.nansum(QC, axis=-3)

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = 'Column Cloud+Rain+Ice'
            save_name = 'Cloud+Rain+Ice'

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(data2D[i], cmap='gist_ncar_r',
                                 color_min=0, color_max=1e-2,
                                 N_colors=20,
                                 cbar_tickformat=OOMFormatter(-3, mathText=False),  #1e-3,
                                 cbar_label="propto kg/m^2")
                #self._add_barbs(U[i], V[i])  # add wind lowest model level

                self._finalize_save(save_name, valid_datetime)


        elif var == 'Updrafts':
            # p = getvar(self.data, "pressure", timeidx=ALL_TIMES, method="cat")
            # z = getvar(self.data, "z", timeidx=ALL_TIMES, units="dm", method="cat")
            U = getvar(self.data, "ua", timeidx=ALL_TIMES, units="kt", method="cat")[:,0,:,:]
            V = getvar(self.data, "va", timeidx=ALL_TIMES, units="kt", method="cat")[:,0,:,:]

            # U = to_np(interplevel(ua, p, 700))
            # V = to_np(interplevel(va, p, 700))
            #print U, V

            objective_resolution = 33 # wind barbs
            mean_xpoints = self.XY_gridpoints**(.5)
            nbarbs = 1 + int(mean_xpoints/objective_resolution)

            data = _load('W')
            data2D = np.nanmax(data, axis=-3)

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = 'Max. Updraft velocity'
            save_name = 'Updraft'

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                #self._add_barbs(U[i], V[i])  # add wind lowest model level

                U1 = U[i]
                V1 = V[i]

                self.bm.barbs(self.x[::nbarbs,::nbarbs], self.y[::nbarbs,::nbarbs],
                         U1[::nbarbs, ::nbarbs],
                         V1[::nbarbs, ::nbarbs],
                         zorder=2, length=5, linewidth=0.5)

                self._plot_heatmap(data2D[i],
                                   interpolation = 'nearest',
                                   cmap='gist_ncar_r',
                                   color_min=0, color_max=20,
                                   N_colors=20,
                                   cbar_tickformat="%d",  #1e-3,
                                   cbar_label="m/s")

                self._finalize_save(save_name, valid_datetime)
        elif var == 'BT':
            OLR = getvar(self.data, "OLR", timeidx=ALL_TIMES, method="cat")

            data2D = (OLR/(5.67*1e-8))**(1/4)-273.15

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = avail[var]
            save_name = 'BT'

            # cmap = plt.get_cmap('gist_ncar_r')
            # cmap = truncate_colormap(cmap, 0, 0.9)

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(data2D[i],
                                   interpolation = 'bilinear',
                                   cmap='gist_gray_r',
                                   color_min=-60, color_max=20,
                                   N_colors=20,
                                   cbar_tickformat="%d",  #1e-3,
                                   cbar_label="degC")

                self._finalize_save(save_name, valid_datetime)

        elif var == 'UDHEL':
            U = getvar(self.data, "ua", timeidx=ALL_TIMES, units="kt", method="cat")[:,0,:,:]
            V = getvar(self.data, "va", timeidx=ALL_TIMES, units="kt", method="cat")[:,0,:,:]

            objective_resolution = 33 # wind barbs
            mean_xpoints = self.XY_gridpoints**(.5)
            nbarbs = 1 + int(mean_xpoints/objective_resolution)
            data2D = wrf.helicity.get_uh(self.data, timeidx=wrf.ALL_TIMES)

            color_min = np.nanmin(data2D)
            color_max = np.nanmax(data2D)
            print 'data min: '+str(color_min)+' max: '+str(color_max)

            figtitle = 'Updraft Helicity'
            save_name = 'UDHEL'

            # cmap = plt.get_cmap('gist_ncar_r')
            # cmap = truncate_colormap(cmap, 0, 0.9)

            N_images = _check_2D(data2D)

            for i in range(N_images):
                valid_datetime = self.wrfdate[i]
                self._make_layout(figtitle, valid_datetime)

                self._plot_heatmap(data2D[i],
                                   interpolation = 'bilinear',
                                   cmap='viridis',
                                   color_min=0, color_max=250,
                                   N_colors=25,
                                   cbar_tickformat="%d",  #1e-3,
                                   cbar_label="m$^2$s$^{-2}$")

                self._finalize_save(save_name, valid_datetime)
        else:
            raise NotImplementedError(str(var)+
                                      'not available. Use one of these:'+
                                      str(avail))


    """Old functions"""
    def plot_SurfLay(self):
        hfx = to_np(getvar(self.data, "HFX", timeidx=ALL_TIMES, method="cat"))             # 10m u-component of wind [ms^-1]
        lhflux = to_np(getvar(self.data, "LH", timeidx=ALL_TIMES, method="cat"))
        grdflux = to_np(getvar(self.data, "GRDFLX", timeidx=ALL_TIMES, method="cat"))

        for i in range(len(self.data)):
            valid_datetime = self.wrfdate[i]

            figtitle = 'Upward Sensible Heat Flux at the surface'
            save_name = 'HFX'
            self._make_layout(figtitle, valid_datetime)
            self._makeplot_sfcflux(i, hfx[i], cmap='RdBu_r')
            self._finalize_save(save_name, valid_datetime)

            figtitle = 'Upward Latent Heat Flux at the surface'
            save_name = 'LHF'
            self._make_layout(figtitle, valid_datetime)
            self._makeplot_sfcflux(i, lhflux[i], cmap='RdBu')
            self._finalize_save(save_name, valid_datetime)

            figtitle = 'Downward Ground Heat Flux'
            save_name = 'GRDFLX'
            self._make_layout(figtitle, valid_datetime)
            self._makeplot_sfcflux(i, grdflux[i], cmap='RdBu_r')
            self._finalize_save(save_name, valid_datetime)


class Compare_Experiments(Load_Domain):
    """Compares fields wrfout.nc files

    EasyWRF.Compare_Experiments(('./case1-wrfouts/', './case2-wrfouts/'),
                                'testdomain',
                                variables = ['OLR', ],
                                )




    """
    def __init__(self,
                 controlrun = './control/',
                 others = ['./experiment1/', ],
                 domainname = "MyDomain",
                 variables = ['T2', ],
                 coast=True,
                 river=False,
                 austrianborders=True,
                 bezirk=True):

        # to clarify
        WRFout_controlrun_data_directory = controlrun
        WRFout_data_directories_list = others

        # load objects
        WRF_objects = [Load_Domain(eachdir, domainname,
                                coast=True,
                                river=False,
                                austrianborders=True,
                                bezirk=True) for eachdir in WRFout_data_directories_list]

        WRF_obj_cntrl = Load_Domain(WRFout_controlrun_data_directory, domainname,
                                coast=coast,
                                river=river,
                                austrianborders=austrianborders,
                                bezirk=bezirk)

        self.cntrl = WRF_obj_cntrl

        for dirname, WRF_object in zip(WRFout_data_directories_list, WRF_objects):
            print 'Comparing '+dirname+' to controlrun ...'
            case = str(dirname.split('/')[-2])

            for var in variables:
                print 'Variable: '+var
                figtitle = WRF_obj_cntrl.data[0][var].description+': '+case+'-Control'

                cntrl_data = WRF_obj_cntrl._load_var(var)
                data2 = WRF_object._load_var(var)
                dif = data2 - cntrl_data

                """Plot the control run of the experiment
                to compare with
                """
                color_range = (np.nanmin(cntrl_data), np.nanmax(cntrl_data))
                print 'data min: '+str(color_range[0])+' max: '+str(color_range[1])
                color_max = color_range[1]
                color_min = color_range[0]

                # schraffieren
                contour_boundaries = (np.nanmin(dif),
                                      np.percentile(dif,10),
                                      np.percentile(dif,20),
                                      np.percentile(dif,80),
                                      np.percentile(dif,90),
                                      np.nanmax(dif))

                hatches = ['///', '//', None, '+', '++']

                save_name = var+'_cntrl'
                cbar_label = WRF_obj_cntrl.data[0][var].units

                def mk_txt(a,b,c):
                    return ' '.join([a,str(np.round(b,1)),c])

                txt0 = '$\Delta >$'
                txt = mk_txt(txt0,contour_boundaries[4],cbar_label)
                hh4 = mpatches.Patch(facecolor='w',hatch=hatches[4],label=txt)
                txt = mk_txt(txt0,contour_boundaries[3],cbar_label)
                hh3 = mpatches.Patch(facecolor='w',hatch=hatches[3],label=txt)
                txt = mk_txt(txt0,contour_boundaries[2],cbar_label)
                hh2 = mpatches.Patch(facecolor='w',hatch=hatches[1],label=txt)
                txt = mk_txt(txt0,contour_boundaries[1],cbar_label)
                hh1 = mpatches.Patch(facecolor='w',hatch=hatches[0],label=txt)

                N_images = _check_2D(cntrl_data)

                for i in range(N_images):
                    valid_datetime = WRF_obj_cntrl.wrfdate[i]
                    WRF_obj_cntrl._make_layout(figtitle, valid_datetime)

                    WRF_obj_cntrl.bm.contourf(WRF_obj_cntrl.x, WRF_obj_cntrl.y, dif[i],
                                     contour_boundaries, colors='none',
                                     hatches=hatches, zorder=5)

                    WRF_obj_cntrl._plot_heatmap(cntrl_data[i],
                                       interpolation = 'nearest',
                                       cmap='viridis',
                                       color_min=color_min, color_max=color_max,
                                       N_colors=20,
                                       cbar_tickformat="%d",  #1e-3,
                                       cbar_label=cbar_label)

                    WRF_obj_cntrl.ax.legend(handles = [hh1,hh2,hh3,hh4], fontsize=9, loc=2)

                    self._save(case, save_name, valid_datetime)

                """ Difference Plots"""
                color_range = (np.nanmin(dif), np.nanmax(dif))
                print 'data min: '+str(color_range[0])+' max: '+str(color_range[1])
                color_max  = max(abs(color_range[1]),abs(color_range[0]))
                if color_max == 0:
                    print 'fields are the same - difference is zero!'
                color_min = -color_max

                save_name = var+'_dif'
                cbar_label = WRF_obj_cntrl.data[0][var].units

                N_images = _check_2D(dif)

                for i in range(N_images):
                    valid_datetime = WRF_obj_cntrl.wrfdate[i]
                    WRF_obj_cntrl._make_layout(figtitle, valid_datetime)

                    WRF_obj_cntrl._plot_heatmap(dif[i],
                                       interpolation = 'nearest',
                                       cmap='RdBu_r',
                                       color_min=color_min, color_max=color_max,
                                       N_colors=20,
                                       cbar_tickformat="%d",  #1e-3,
                                       cbar_label=cbar_label)

                    self._save(case, save_name, valid_datetime)

    def _save(self, case, var_name, dt_valid):
        init = datetime.strptime(self.cntrl.init, '%Y-%m-%d_%H:%M:%S')
        savedir = './img/'+init.strftime('%Y%m%d%H')+'/'+case+'/'
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        filename = var_name + '_' + dt_valid.strftime('%y%m%d_%H%M') + '.png'
        plt.savefig(savedir + filename, dpi=100)
        plt.close('all')
        print filename, 'saved'
