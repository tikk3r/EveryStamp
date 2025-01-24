"""Sub-module for plotting FITS images."""
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy
from aplpy import FITSFigure
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.pyplot import figure
import numpy as np
from typing import Union
import astropy.units as u


def find_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

    maskSup = 1e-7
    m = image_data[np.abs(image_data)>maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold*cut)[0]
        rms = np.std(m[ind])
        if np.abs((rms-rmsold) / rmsold) < diff: break
        rmsold = rms
    print(f'Noise : {str(round(rms * 1000, 4))} {u.mJy/u.beam}')
    return rms

class BasicFITSPlot:
    """Creates a basic plot of a FITS file."""

    def __init__(self, fitsname):
        """Initialise a basic plotting object for 2D FITS files.

        Args:
            fitsname : str
                Name of the FITS file that will be plotted.
        """
        self.dpi = 300
        self.figsize = (12, 8)
        self.fitsimage = fitsname
        self.fitsdata = fits.getdata(fitsname).squeeze()
        if len(self.fitsdata.shape) > 2:
            raise NotImplementedError(
                "Supplied FITS file appears to have dimensions besides RA and DEC. BasicPlot only supports 2D FITS files."
            )
        self.data = self.fitsdata
        self.wcs = WCS(fits.getheader(fitsname)).celestial
        self.figsize = [
            self.fitsdata.shape[0] // self.dpi,
            self.fitsdata.shape[1] // self.dpi,
        ]
        if self.figsize[0] < 12:
            self.figsize[0] = 12
        if self.figsize[1] < 8:
            self.figsize[1] = 8

    def plot2D(
        self,
        plot_colourbar=False,
        contour_image: Optional[numpy.ndarray] = None,
        contour_levels: Union[int, list] = 7,
        cmap_min: Optional[float] = None,
        cmap_max: Optional[float] = None,
        cmap: str = "grey",
    ):
        """Save a 2D plot of the loaded FITS file.

        Args:
            plot_colourbar : bool
                Add a colour bar to the plot.
            contour_image:
                Add contours based on this image.
            contour_levels:
                Number of contour levels to draw if an integer or contour levels if a list. Defaults to 5.
        """
        hdu = fits.PrimaryHDU(
            header=WCS(fits.getheader(self.fitsimage)).celestial.to_header(),
            data=self.data.squeeze(),
        )
        f = FITSFigure(hdu, figsize=self.figsize)
        if not cmap:
            f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
        else:
            print(f'Using colour map: {cmap}')
            f.show_colorscale(vmin=cmap_min, vmax=cmap_max, pmax=100, cmap=cmap)
        if contour_image:
            head = fits.getheader(contour_image)
            data = fits.getdata(contour_image)
            head = WCS(head).celestial.to_header()
            hdu_c = fits.PrimaryHDU(data=data, header=head)
            # f.show_contour(hdu_c, levels=contour_levels, colors='white', cmap='plasma')
            f.show_contour(hdu_c, levels=contour_levels, colors="C0")
        if plot_colourbar:
            f.add_colorbar()
        f.savefig(self.fitsimage.replace("fits", "png"), dpi=self.dpi)

    def plot_noaxes(self, cmap_min: float = None, cmap_max: float = None, cmap=None):
        """Save a plot of the FITS image without any axes."""
        figsize = [
            self.fitsdata.shape[0] // self.dpi,
            self.fitsdata.shape[1] // self.dpi,
        ]
        if figsize[0] < 12:
            figsize[0] = 12
        if figsize[1] < 8:
            figsize[1] = 8
        fig = figure(figsize=figsize)
        hdu = fits.PrimaryHDU(
            header=WCS(fits.getheader(self.fitsimage)).celestial.to_header(),
            data=self.data.squeeze(),
        )
        f = FITSFigure(hdu, figure=fig)
        if not cmap:
            f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
        else:
            f.show_colorscale(vmin=cmap_min, vmax=cmap_max, pmax=100, cmap=cmap)
        plt.axis("off")
        fig.savefig(
            self.fitsimage.replace("fits", ".noaxes.png"),
            bbox_inches="tight",
            pad_inches=0,
            transparent=True,
            dpi=self.dpi,
        )

    def savedata(self, outfile):
        """Save data of a BasicPlot object to a FITS file with the same WCS information.

        Any tonemapping applied to the original data will be carried over to the FITS file. Physical units will thus be lost.
        """
        fits.writeto(
            data=self.data,
            header=self.wcs.to_header(),
            filename=outfile,
            overwrite=True,
        )


class SRTPlot:
    """Create a line profile plot similar in style to the SRTPLOT task used by Hogbom 1974."""

    def __init__(self, fitsname):
        """Initialise a basic plotting object for 2D FITS files.

        Args:
            fitsname : str
                Name of the FITS file that will be plotted.
        """
        self.dpi = 300
        self.figsize = (12, 12)
        self.fitsimage = fitsname
        self.fitsdata = fits.getdata(fitsname).squeeze()
        if len(self.fitsdata.shape) > 2:
            raise NotImplementedError(
                "Supplied FITS file appears to have dimensions besides RA and DEC. BasicPlot only supports 2D FITS files."
            )
        self.data = self.fitsdata

    def plot2D(self, srt_lines: int = 25, srt_offset: float = 0.01, **kwargs):
        """Save an SRTPLOT style plot of the loaded FITS file.

        Args:
            lines : int
                Number of lines to plot.
            offset : float
                Offset between each line in data units.
        """
        f = plt.figure(figsize=self.figsize)
        ax = f.add_subplot(111)
        stride = self.data.shape[0] // srt_lines
        if stride == 0:
            stride += 1
        for i, line in enumerate(self.data[::stride, ::-1]):
            if i > 0:
                ax.fill_between(
                    list(range(len(line))),
                    (srt_lines - i) * srt_offset,
                    line + (srt_lines - i) * srt_offset,
                    color="w",
                    zorder=i + 2,
                )
            else:
                ax.fill_between(
                    list(range(len(line))),
                    0,
                    line + (srt_lines - i) * srt_offset,
                    color="w",
                    zorder=i + 2,
                )
            ax.plot(line + (srt_lines - i) * srt_offset, color="k", zorder=i + 2)
        ax.set_xlim(0, self.data.shape[1])
        plt.axis("off")
        f.savefig(
            self.fitsimage.replace(".fits", ".srtplot.png"),
            dpi=self.dpi,
            bbox_inches="tight",
            pad_inches=0,
        )


class BasicImagePlot:
    """ Creates a basic plot of a FITS file."""
    def __init__(self, imname, wcsimage=None):
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
        if wcsimage:
            self.wcs = WCS(fits.getheader(wcsimage)).celestial
        else:
            self.wcs = None

    def plot2D(self, plot_colourbar=False, contour_image: numpy.ndarray = None, contour_levels: Union[int, list] = 5, cmap: str = "gray", cmap_min: float = None, cmap_max: float = None):
        """ Save a plot of the FITS image without any axes."""
        figsize = [self.imdata.shape[0] // self.dpi, self.imdata.shape[1] // self.dpi]
        if figsize[0] < self.imdata.shape[0]:
            figsize[0] = 4
        if figsize[1] < self.imdata.shape[1]:
            figsize[1] = 4
        print("FIGURE SIZE: ", figsize)
        fig = figure(figsize=figsize, dpi=self.dpi)
        if self.wcs:
            ax = fig.add_subplot(111, projection=self.wcs)
            if self.data is not None:
                ax.imshow(self.data, interpolation='none', cmap=cmap)
            else:
                ax.imshow(self.imdata, interpolation='none', cmap=cmap)
        else:
            ax = fig.add_subplot(111)
            if self.data is not None:
                ax.imshow(self.data, origin='upper', interpolation='none', cmap=cmap)
            else:
                ax.imshow(self.imdata, origin='upper', interpolation='none', cmap=cmap)
        if contour_image:
            # Flip to get North up.
            cdata = np.flipud(fits.getdata(contour_image).squeeze())
            chead = fits.getheader(contour_image)
            wcs = WCS(chead).celestial
            # f.show_contour(hdu_c, levels=contour_levels, colors='white', cmap='plasma')
            crms = find_rms(cdata)
            clevels = np.arange(crms, np.percentile(cdata, 99.9), np.sqrt(2) * 70e-6)
            if type(contour_levels) is int:
                step = len(clevels) // contour_levels
                clevels = clevels[::step]
            ax.contour(cdata, levels=clevels, colors='w', transform=ax.get_transform(wcs), linewidths=2)
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)
        #plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.xlabel("Right ascension [J2000]", fontsize=16)
        plt.ylabel("Declination [J2000]", fontsize=16)
        #plt.margins(0, 0)
        #plt.gca().xaxis.set_major_locator(plt.NullLocator())
        #plt.gca().yaxis.set_major_locator(plt.NullLocator())
        file_ext = "." + self.image.split(".")[-1]
        fig.savefig(
            self.image.replace(file_ext, ".output.png"),
            bbox_inches="tight",
            transparent=True,
            dpi=self.dpi,
        )
