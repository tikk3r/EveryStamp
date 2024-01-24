"""Sub-module for trimming FITS images."""

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import pyregion
from numpy import array
from typing import Union
import sys


def make_cutout_2D(image: str = None, pos: Union[list, tuple] = None,
                   size: Union[list, tuple] = None, outfile: str = None):
    """
    Make 2D cutout with astropy
    ---------------------------
    :param image: image name
    :param pos: central position
    :param size: image size
    :param outfile: output fits file
    """
    head = fits.getheader(image)
    data = fits.getdata(image).squeeze()
    if len(data.shape) > 2:
        raise ValueError(
            "Image data has dimensions other than RA and DEC, which is not yet supported."
        )
    wcs = WCS(head).celestial

    cutout = Cutout2D(data, pos, size, wcs=wcs)
    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(outfile)


def make_cutout_region(image: str = None, region: str = None, outfile: str = None):
    """
    Make 2D cutout with pyregion
    ---------------------------
    :param image: image name
    :param region: region file
    :param outfile: output fits file
    """

    hdu = fits.open(image)
    head = hdu[0].header
    data = hdu[0].data

    while data.ndim > 2:
        data = data[0]

    r = pyregion.open(region).as_imagecoord(header=head)
    mask = r.get_mask(hdu=hdu[0], shape=(head["NAXIS1"], head["NAXIS2"]))

    assert len(r) == 1, 'Multiple regions in 1 file given, only one allowed'

    shape = array(r[0].coord_list)

    # circle
    if len(shape) == 3:
        cutout = Cutout2D(data=data * mask,
                          position=(shape[0], shape[1]),
                          size=(shape[2], shape[2]),
                          wcs=WCS(head, naxis=2),
                          mode='partial')
    # square
    elif len(shape) > 3:
        print(shape)

        cutout = Cutout2D(data=data * mask,
                          position=(shape[0], shape[1]),
                          size=(shape[3], shape[2]),
                          wcs=WCS(head, naxis=2),
                          mode='partial')
    else:
        sys.exit("ERROR: Should not arrive here")

    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(outfile)
