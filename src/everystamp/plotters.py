"""Sub-module for plotting FITS images."""
import os

from aplpy import FITSFigure
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.image import imread
from matplotlib.pyplot import figure, show
from typing import Union

import matplotlib.pyplot as plt
import numpy


class BasicFITSPlot():
    """ Creates a basic plot of a FITS file."""
    def __init__(self, fitsname):
        """ Initialise a basic plotting object for 2D FITS files.
        
        Args:
            fitsname : str
                Name of the FITS file that will be plotted.
        """
        self.dpi = 300
        self.figsize = (12, 8)
        self.fitsimage = fitsname
        self.fitsdata = fits.getdata(fitsname).squeeze()
        if len(self.fitsdata.shape) > 2:
            raise NotImplementedError('Supplied FITS file appears to have dimensions besides RA and DEC. BasicPlot only supports 2D FITS files.')
        self.data = self.fitsdata
        self.wcs = WCS(fits.getheader(fitsname)).celestial
        self.figsize = [self.fitsdata.shape[0] // self.dpi, self.fitsdata.shape[1] // self.dpi]
        if self.figsize[0] < 12:
            self.figsize[0] = 12
        if self.figsize[1] < 8:
            self.figsize[1] = 8

    def plot2D(self, plot_colourbar=False, contour_image: numpy.ndarray = None, contour_levels: Union[int, list] = 7, cmap_min: float = None, cmap_max: float = None):
        """ Save a 2D plot of the loaded FITS file.

        Args:
            plot_colourbar : bool
                Add a colour bar to the plot.
            contour_image:
                Add contours based on this image.
            contour_levels:
                Number of contour levels to draw if an integer or contour levels if a list. Defaults to 5.
        """
        hdu = fits.PrimaryHDU(header=fits.getheader(self.fitsimage), data=self.data)
        f = FITSFigure(hdu, figsize=self.figsize)
        f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
        if contour_image:
                hdu_c = fits.open(contour_image)
                # f.show_contour(hdu_c, levels=contour_levels, colors='white', cmap='plasma')
                f.show_contour(hdu_c, levels=contour_levels, colors='C0')
        if plot_colourbar:
            plt.colorbar(im)
        f.savefig(self.fitsimage.replace('fits', 'png'), dpi=self.dpi)

    def plot_noaxes(self, cmap_min: float = None, cmap_max: float = None):
        """ Save a plot of the FITS image without any axes."""
        hdu = fits.PrimaryHDU(header=fits.getheader(self.fitsimage), data=self.data)
        f = FITSFigure(hdu, figsize=self.figsize)
        f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
        f.axis_labels.hide()
        f.tick_labels.hide()
        f.ticks.hide()
        f.savefig(self.fitsimage.replace('.fits', '.noaxes.png'), dpi=self.dpi)#bbox_inches='tight', pad_inches=0, transparent=True, dpi=self.dpi)
        return
        fig = figure(figsize=figsize)
        hdu = fits.PrimaryHDU(header=WCS(fits.getheader(self.fitsimage)).celestial.to_header(), data=self.data.squeeze())
        f = FITSFigure(hdu, figure=fig)
        f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
        plt.axis('off')
        fig.savefig(self.fitsimage.replace('fits', '.noaxes.png'), bbox_inches='tight', pad_inches=0, transparent=True, dpi=self.dpi)

    def savedata(self, outfile):
        """ Save data of a BasicPlot object to a FITS file with the same WCS information.
        
        Any tonemapping applied to the original data will be carried over to the FITS file. Physical units will thus be lost.
        """
        fits.writeto(data=self.data, header=self.wcs.to_header(), filename=outfile, overwrite=True)


class BasicImagePlot():
    """ Creates a basic plot of a FITS file."""
    def __init__(self, imname):
        """ Initialise a basic plotting object for 2D FITS files.
        
        Args:
            fitsname : str
                Name of the FITS file that will be plotted.
        """
        self.dpi = 300
        self.figsize = (12, 8)
        self.image = imname
        import cv2
        self.imdata = cv2.cvtColor(cv2.imread(imname), cv2.COLOR_BGR2RGB)
        self.data = self.imdata

    def plot2D(self):
        """ Save a plot of the FITS image without any axes."""
        figsize = [self.imdata.shape[0] // self.dpi, self.imdata.shape[1] // self.dpi]
        fig = figure(figsize=figsize)
        ax = fig.add_subplot(111)
        if self.data is not None:
            ax.imshow(self.data, origin='upper', interpolation='none')
            # if self.data.max() > 1:
            #     # Probably integer image.
            #     ax.imshow(self.data.astype('uint8'), origin='upper', interpolation='none')
            # elif self.data.max() <= 1:
            #     # Probably floating point image.
            #     ax.imshow(self.data.astype(float), origin='upper', interpolation='none')
        else:
            ax.imshow(self.imdata, origin='upper', interpolation='none')
            # if self.imdata.max() > 1:
            #     # Probably integer image.
            #     ax.imshow(self.imdata.astype('uint8'), origin='upper', interpolation='none')
            # elif self.imdata.max() <= 1:
            #     # Probably floating point image.
            #     ax.imshow(self.imdata.astype(float), origin='upper', interpolation='none')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        file_ext = '.' + self.image.split('.')[-1]
        fig.savefig(self.image.replace(file_ext, '.output.png'), bbox_inches='tight', pad_inches=0, transparent=True, dpi=self.dpi)