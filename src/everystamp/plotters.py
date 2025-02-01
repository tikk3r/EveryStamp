"""Sub-module for plotting FITS images."""

from typing import List, Optional, Union

import astropy.units as u
import colormaps
import matplotlib.pyplot as plt
import numpy
import aplpy
from aplpy import FITSFigure
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.pyplot import figure
import numpy as np
from aplpy import FITSFigure
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import (
    ImageNormalize,
    MinMaxInterval,
    PercentileInterval,
    LinearStretch,
    LogStretch,
    SqrtStretch,
)
from astropy.wcs import WCS
from blend_modes import (
    addition,
    hard_light,
    soft_light,
    lighten_only,
    darken_only,
    multiply,
    dodge,
    difference,
    subtract,
    grain_extract,
    grain_merge,
    divide,
    overlay,
    normal,
)
from matplotlib import colormaps as mplcm
from matplotlib.pyplot import figure
from PIL import Image


def find_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

    maskSup = 1e-7
    m = image_data[np.abs(image_data) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.0
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs((rms - rmsold) / rmsold) < diff:
            break
        rmsold = rms
    print(f"Noise : {str(round(rms * 1000, 4))} {u.mJy/u.beam}")
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
        self.load_rgb = False
        if len(self.fitsdata.shape) > 2:
            #raise NotImplementedError(
            #    "Supplied FITS file appears to have dimensions besides RA and DEC. BasicPlot only supports 2D FITS files."
            #)
            print("Generating RGB image")
            aplpy.make_rgb_image(self.fitsimage, "temp_rgb.png", embed_avm_tags=True, stretch_r="sqrt", stretch_g="sqrt", stretch_b="sqrt")
            self.load_rgb = True
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
        if self.load_rgb:
            f = FITSFigure("temp_rgb.png", figsize=self.figsize)
            f.show_rgb()
        else:
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
            cdata = fits.getdata(contour_image).squeeze()
            chead = fits.getheader(contour_image)
            crms = find_rms(cdata)
            contour_levels = np.arange(5*crms, np.percentile(cdata, 99.999), np.sqrt(2) * crms)
            hdu_c = fits.PrimaryHDU(data=cdata, header=chead)
            f.show_contour(hdu_c, levels=contour_levels, colors='white')
        if plot_colourbar:
            f.add_colorbar()
        f.savefig(self.fitsimage.replace("fits", "png"), dpi=self.dpi)

    def plot_noaxes(self, cmap_min: float = None, cmap_max: float = None, cmap=None, contour_image=None):
        """Save a plot of the FITS image without any axes."""
        figsize = [
            self.fitsdata.shape[0] // self.dpi,
            self.fitsdata.shape[1] // self.dpi,
        ]
        if figsize[0] < 12:
            figsize[0] = 12
        if figsize[1] < 8:
            figsize[1] = 8
        #fig = figure(figsize=figsize)
        if self.load_rgb:
            f = FITSFigure("temp_rgb.png", figsize=self.figsize)
            f.show_rgb()
        else:
            hdu = fits.PrimaryHDU(
                header=WCS(fits.getheader(self.fitsimage)).celestial.to_header(),
                data=self.data.squeeze(),
            )
            f = FITSFigure(hdu, figsize=self.figsize)
            if not cmap:
                f.show_grayscale(vmin=cmap_min, vmax=cmap_max, pmax=100)
            else:
                f.show_colorscale(vmin=cmap_min, vmax=cmap_max, pmax=100, cmap=cmap)
        if contour_image:
            cdata = fits.getdata(contour_image).squeeze()
            chead = fits.getheader(contour_image)
            crms = find_rms(cdata)
            contour_levels = np.arange(5*crms, np.percentile(cdata, 99.999), np.sqrt(2) * crms)
            hdu_c = fits.PrimaryHDU(data=cdata, header=chead)
            f.show_contour(hdu_c, levels=contour_levels, colors='white')
        plt.axis("off")
        f.savefig(
            self.fitsimage.replace("fits", ".noaxes.png"),
            transparent=True,
            dpi=self.dpi,
        )

    def savedata(self, outfile):
        """Save data of a BasicPlot object to a FITS file with the same WCS information

        Any tonemapping applied to the original data will be carried over to the FITS file. Physical units will thus be lost.
        """
        fits.writeto(
            data=self.data,
            header=self.wcs.to_header(),
            filename=outfile,
            overwrite=True,
        )


class BlendPlot:
    """Creates a composite image using blending modes.

    Attributes:
        backgroundg
        foregroundg
        blend_cmapsg
        centreg
        radiusg
        rmscutg
        blend_modesg
        blend_cmapsg
        blend_opacities

    Methods:
        blend
        load_preset
        prepare_images
        set_blends
    """

    def __init__(
        self,
        background: str,
        foreground: List[str],
        cmaps: List[str],
        centre: SkyCoord,
        radius: float,
        rmscut: float,
    ) -> None:
        """Initialises the BlendPlot.

        Args:
            background (str):
            foreground (str):
            cmaps (list[str]):
            centre (SkyCoord):
            radius (float):
            rmscut (float):
        """
        Image.MAX_IMAGE_PIXELS = None
        self.background = background
        self.foreground = foreground
        self.blend_cmaps = cmaps
        self.centre = centre
        self.radius = radius
        self.rmscut = rmscut

    def prepare_images(self) -> None:
        """Prepare images for blending.

        This involves plotting all the specified images in the same WCS projection, with the specified rms noise cut.
        Plots are saved as PNGs in the current directory.
        """
        print("Preparing background image.")
        fig = FITSFigure(self.background, figsize=(8, 8), dpi=150)
        fig.show_rgb(interpolation="none")
        print(
            f"Recentring on {self.centre.ra.value}, {self.centre.dec.value} {self.radius}"
        )
        fig.recenter(self.centre.ra.value, self.centre.dec.value, self.radius)
        fig.axis_labels.hide()
        fig.tick_labels.hide()
        fig.ticks.hide()
        fig.savefig("temp_background.png")

        for i, fg in enumerate(self.foreground):
            print("Preparing foreground image.")
            h = fits.getheader(fg)
            wcs = WCS(h).celestial
            d = fits.getdata(fg).squeeze()
            rms = find_rms(d)
            d[d < self.rmscut * rms] = np.nan
            norm = ImageNormalize(
                d, interval=PercentileInterval(99.99), stretch=SqrtStretch()
            )

            figf = FITSFigure(self.background, figsize=(8, 8), dpi=150)
            if self.blend_cmaps[i] in list(mplcm):
                cm = self.blend_cmaps[i]
            else:
                cm = eval("colormaps." + self.blend_cmaps[i])
            figf.ax.imshow(
                d,
                transform=figf.ax.get_transform(wcs),
                cmap=cm,
                norm=norm,
                interpolation="none",
            )
            figf.recenter(self.centre.ra.value, self.centre.dec.value, self.radius)
            figf.axis_labels.hide()
            figf.tick_labels.hide()
            figf.ticks.hide()
            figf.savefig(f"temp_foreground_{i:02d}.png", transparent=True)

            fig.ax.imshow(
                d, transform=fig.ax.get_transform(wcs), cmap="afmhot", norm=norm
            )
            fig.savefig(f"temp_reference{i:02d}.png", dpi=150)

        del fig, figf

    def blend(self) -> None:
        """Blend the background and foreground images together using the specified modes.

        Returns:
            None
        """
        img_blend = np.array(Image.open("temp_background.png")).astype(float)
        for i, (bms, bas) in enumerate(zip(self.blend_modes, self.blend_opacities)):
            img_fg = np.array(Image.open(f"temp_foreground_{i:02d}.png")).astype(float)
            if len(bas) == 1:
                opacs = bas * len(bms.split(","))
            elif (len(bas) > 1) and (len(bas) != len(bms.split(","))):
                raise ValueError(
                    "Blend layers and opacities must match in length if nested opacities are given."
                )
            else:
                opacs = bas
            for bm, ba in zip(bms.split(","), opacs):
                print(f"Blending {bm} with {ba} opacity.")
                match bm:
                    case "add":
                        img_blend = addition(img_blend, img_fg, ba)
                    case "softlight":
                        img_blend = soft_light(img_blend, img_fg, ba)
                    case "hardlight":
                        img_blend = hard_light(img_blend, img_fg, ba)
                    case "lighten_only":
                        img_blend = lighten_only(img_blend, img_fg, ba)
                    case "darken_only":
                        img_blend = darken_only(img_blend, img_fg, ba)
                    case "multiply":
                        img_blend = multiply(img_blend, img_fg, ba)
                    case "dodge":
                        img_blend = dodge(img_blend, img_fg, ba)
                    case "difference":
                        img_blend = difference(img_blend, img_fg, ba)
                    case "subtract":
                        img_blend = subtract(img_blend, img_fg, ba)
                    case "grain_extract":
                        img_blend = grain_extract(img_blend, img_fg, ba)
                    case "grain_merge":
                        img_blend = grain_merge(img_blend, img_fg, ba)
                    case "divide":
                        img_blend = divide(img_blend, img_fg, ba)
                    case "overlay":
                        img_blend = overlay(img_blend, img_fg, ba)
                    case "normal":
                        img_blend = normal(img_blend, img_fg, ba)
        Image.fromarray(np.uint8(img_blend)).save(
            self.foreground[0].replace(".fits", "_blend.png")
        )

    def load_preset(self, preset: str) -> None:
        """Loads a preset of blending modes, colour maps and opacities.

        Args:
            preset (str): name of the preset to load.

        Returns:
            None

        Raises:
            ValueError: if an unknown preset is specified.

        """
        match preset:
            case "opt+x-ray+lofar":
                self.blend_modes = ["add,softlight", "add,add,overlay"]
                self.blend_cmaps = ["c_7_16", "solar"]
                self.blend_opacities = [0.35, 0.6]
                self.rmscut = 5.0
            case "opt+lofar_hot":
                self.blend_modes = ["add,add"]
                self.blend_cmaps = ["afmhot"]
                self.blend_opacities = [0.6]
                self.rmscut = 5.0
            case "opt+lofar_solar":
                self.blend_modes = ["add,softlight"]
                self.blend_cmaps = ["solar"]
                self.blend_opacities = [0.6]
                self.rmscut = 5.0
            case _:
                raise ValueError("Unknown preset requested.")

    def set_blends(
        self,
        blend_modes: List[str],
        blend_cmaps: List[str],
        blend_opacities: List[float | List[float]],
    ) -> None:
        """Set the blending modes, colour maps and opacities of the current plot.

        Args:
            blend_modes (list[str]): blending modes to apply to each image.
            blend_cmaps (list[str]): colour maps to use for each image.
            blend_opacities (list[float]): opacities with which to blend each layer.

        Returns:
            None
        """
        self.blend_modes = blend_modes
        self.blend_cmaps = blend_cmaps
        self.blend_opacities = blend_opacities


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
    """Creates a basic plot of an image (not FITS) file."""

    def __init__(self, imname, wcsimage=None):
        """Initialise a basic plotting object for 2D FITS files.

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
        self.wcsimage = wcsimage
        if wcsimage:
            print(f"Loading WCS from {wcsimage}")
            self.wcs = WCS(fits.getheader(wcsimage)).celestial
        else:
            self.wcs = None

    def plot2D(
        self,
        plot_colourbar=False,
        contour_image: numpy.ndarray = None,
        contour_levels: Union[int, list] = 5,
        cmap: str = "gray",
        cmap_min: float = None,
        cmap_max: float = None,
    ):
        """Save a plot of the FITS image without any axes."""
        figsize = [self.imdata.shape[0] // self.dpi, self.imdata.shape[1] // self.dpi]
        if figsize[0] < self.imdata.shape[0]:
            figsize[0] = 8
        if figsize[1] < self.imdata.shape[1]:
            figsize[1] = 8
        print("FIGURE SIZE: ", figsize)
        f = FITSFigure(self.wcsimage, figsize=self.figsize, dimensions=[0, 1], slices=[0])
        f.show_rgb(self.image)
        if contour_image:
            ## Flip to get North up.
            ##cdata = np.flipud(fits.getdata(contour_image).squeeze())
            cdata = fits.getdata(contour_image).squeeze()
            chead = fits.getheader(contour_image)
            crms = find_rms(cdata)
            clevels = np.arange(5*crms, np.percentile(cdata, 99.999), np.sqrt(2) * crms)
            f.show_contour(contour_image, levels=clevels, colors="w")
        f.savefig(self.image + "_plot.png", dpi=self.dpi)
        return
        fig = figure(figsize=figsize, dpi=self.dpi)
        if self.wcs:
            ax = fig.add_subplot(111, projection=self.wcs)
            #if self.data is not None:
            #    ax.imshow(np.flipud(self.data), interpolation='none', cmap=cmap)
            #else:
            #    ax.imshow(np.flipud(self.imdata), interpolation='none', cmap=cmap)
        else:
            raise ValueError("BOOM")
            ax = fig.add_subplot(111)
            if self.data is not None:
                ax.imshow(
                    np.flipud(self.data),
                    origin="lower",
                    interpolation="none",
                    cmap=cmap,
                )
            else:
                ax.imshow(
                    np.flipud(self.imdata),
                    origin="lower",
                    interpolation="none",
                    cmap=cmap,
                )
        if contour_image:
            # Flip to get North up.
            # cdata = np.flipud(fits.getdata(contour_image).squeeze())
            cdata = fits.getdata(contour_image).squeeze()
            chead = fits.getheader(contour_image)
            wcs = WCS(chead).celestial
            crms = find_rms(cdata)
            crms = 30e-6
            clevels = np.arange(5*crms, np.percentile(cdata, 99.999), np.sqrt(2) * crms)
            print("Contour levels:", clevels)
            if type(contour_levels) is int:
                step = len(clevels) // contour_levels
                if step > 0:
                    clevels = clevels[::step]
                else:
                    clevels = 5
            ax.contour(cdata, levels=clevels, colors='k', transform=ax.get_transform(wcs), linewidths=2, zorder=999)
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)
        #plt.gca().set_axis_off()
        #plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.xlabel("Right ascension [J2000]", fontsize=16)
        plt.ylabel("Declination [J2000]", fontsize=16)
        # plt.margins(0, 0)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        file_ext = "." + self.image.split(".")[-1]
        fig.savefig(
            self.image.replace(file_ext, ".output.png"),
            bbox_inches="tight",
            transparent=True,
            dpi=self.dpi,
        )
