# -*- coding: utf-8 -*-

import time
import numpy as np
import matplotlib.ticker
import matplotlib.colors as mc

class Timer(object):
    """Provides simple timing of code for debugging.

    Example call:
        with Timer('Plotting Figure1'):
            plt.plot( some stuff )  # plotting
    """

    def __init__(self, name=' '):
        self.name = name

    def __enter__(self):
        self.elapsed = -time.time()

    def __exit__(self, type, value, traceback):
        print '['+str(self.name)+'] Elapsed: '+str(self.elapsed + time.time())+' secs'


class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    """Format colorbar by scientific notation 10^3 etc

    Example call:
        cbar_tickformat = OOMFormatter(-3, mathText=False)
    """

    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

def truncate_colormap(cmap, minval=0.0, color_max=1.0, n=100):
    """Create new colormap from a given one. From stackexchange.

    Trunkates colormaps somewhere btw older minimum color (0) and maximum color (1).
    If color_max>1, the last color extends after 1.

    Example:
        new_cmap = truncate_colormap(plt.get_cmap('terrain'),
                                     0.1, 0.75)
    """
    new_cmap = mc.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=color_max),
        cmap(np.linspace(minval, color_max, n)))
    return new_cmap
