'''Sub-module for plotting FITS images.'''
import os

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.pyplot import figure, show

import matplotlib.pyplot as plt


class BasicPlot():
    def __init__(self, fitsname):
        self.dpi = 300
        self.fitsimage = fitsname
        self.fitsdata = fits.getdata(fitsname).squeeze()
        self.data = self.fitsdata
        self.wcs = WCS(fits.getheader(fitsname)).celestial

    def plot2D(self):
        figsize = (self.fitsdata.shape[0] // self.dpi, self.fitsdata.shape[1] // self.dpi)
        fig = figure(figsize=figsize, dpi=self.dpi)
        ax = fig.add_subplot(111, projection=self.wcs)
        if self.data is not None:
            print('bla')
            im = ax.imshow(self.data, origin='lower', interpolation='none')
        else:
            ax.imshow(self.fitsdata, origin='lower', interpolation='none')
        ax.set(xlabel='Right ascension', ylabel='Declination')
        plt.colorbar(im)
        fig.savefig(self.fitsimage.replace('fits', 'png'), bbox_inches='tight', dpi=self.dpi)

    def savedata(self, outfile):
        fits.writeto(data=self.data, header=self.wcs.to_header(), filename=outfile, overwrite=True)