"""Sub-module for trimming FITS images."""

from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.wcs import WCS

def make_cutout_2D(image, pos, size, outfile):
    head = fits.getheader(image)
    data = fits.getdata(image).squeeze()
    if len(data.shape) > 2:
        raise ValueError('Image data has dimensions other than RA and DEC, which is not yet supported.')
    wcs = WCS(head).celestial

    cutout = Cutout2D(data, pos, size, wcs=wcs)
    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(outfile)

