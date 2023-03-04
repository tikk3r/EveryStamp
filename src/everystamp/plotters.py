'''Sub-module for plotting FITS images.'''
import os

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.pyplot import figure, show

import matplotlib.pyplot as plt


class BasicPlot():
    ''' Creates a basic plot of a FITS file.'''
    def __init__(self, fitsname):
        ''' Initialise a basic plotting object for 2D FITS files.
        
        Args:
            fitsname : str
                Name of the FITS file that will be plotted.
        '''
        self.dpi = 300
        self.figsize = (12, 8)
        self.fitsimage = fitsname
        self.fitsdata = fits.getdata(fitsname).squeeze()
        if len(self.fitsdata.shape) > 2:
            raise NotImplementedError('Supplied FITS file appears to have dimensions besides RA and DEC. BasicPlot only supports 2D FITS files.')
        self.data = self.fitsdata
        self.wcs = WCS(fits.getheader(fitsname)).celestial

    def plot2D(self, plot_colourbar=False):
        ''' Save a 2D plot of the loaded FITS file.

        Args:
            plot_colourbar : bool
                Add a colour bar to the plot.
        '''
        figsize = [self.fitsdata.shape[0] // self.dpi, self.fitsdata.shape[1] // self.dpi]
        if figsize[0] < 12:
            figsize[0] = 12
        if figsize[1] < 8:
            figsize[1] = 8
        fig = figure(figsize=figsize, dpi=self.dpi)
        ax = fig.add_subplot(111, projection=self.wcs)
        if self.data is not None:
            im = ax.imshow(self.data, origin='lower', interpolation='none')
        else:
            ax.imshow(self.fitsdata, origin='lower', interpolation='none')
        ax.set(xlabel='Right ascension', ylabel='Declination')
        if plot_colourbar:
            plt.colorbar(im)
        fig.savefig(self.fitsimage.replace('fits', 'png'), bbox_inches='tight', dpi=self.dpi)

    def plot_noaxes(self):
        ''' Save a plot of the FITS image without any axes.'''
        figsize = [self.fitsdata.shape[0] // self.dpi, self.fitsdata.shape[1] // self.dpi]
        fig = figure(figsize=figsize)
        ax = fig.add_subplot(111)
        if self.data is not None:
            ax.imshow(self.data, origin='lower', interpolation='none')
        else:
            ax.imshow(self.fitsdata, origin='lower', interpolation='none')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        fig.savefig(self.fitsimage.replace('.fits', '.noaxes.png'), bbox_inches='tight', pad_inches=0, transparent=True, dpi=self.dpi)

    def savedata(self, outfile):
        ''' Save data of a BasicPlot object to a FITS file with the same WCS information.
        
        Any tonemapping applied to the original data will be carried over to the FITS file. Physical units will thus be lost.
        '''
        fits.writeto(data=self.data, header=self.wcs.to_header(), filename=outfile, overwrite=True)