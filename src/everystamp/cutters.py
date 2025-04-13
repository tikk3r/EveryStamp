"""Sub-module for trimming FITS images."""

from typing import Optional, Union

import numpy as np
import pyregion

# import pyregion
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.units import Quantity
from astropy.wcs import WCS
from numpy import array


def make_cutout_2D(
    image: str,
    pos: SkyCoord,
    size: Quantity,
    outfile: str,
):
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


def make_cutout_2D_fast(
    image: str,
    pos: SkyCoord,
    size: Quantity,
    outfile: str,
):
    with fits.open(image) as hdul:
        hdu = hdul[0]
        wcs = WCS(hdu.header).celestial

        pix_pos = wcs.world_to_pixel(pos)

        width_pix, height_pix = wcs.world_to_pixel(
            SkyCoord(pos.ra + size, pos.dec + size)
        )
        width_pix -= pix_pos[0]
        height_pix -= pix_pos[1]

        # Calculate the bounding box in pixel coordinates
        x_min = int(np.floor(pix_pos[0] - width_pix / 2))
        x_max = int(np.ceil(pix_pos[0] + width_pix / 2))
        y_min = int(np.floor(pix_pos[1] - height_pix / 2))
        y_max = int(np.ceil(pix_pos[1] + height_pix / 2))

        cutout_data = hdu.data[..., y_min:y_max, x_max:x_min]

        hdu.header["NAXIS1"] = x_max - x_min
        hdu.header["NAXIS2"] = y_max - y_min
        hdu.header["CRPIX1"] -= x_max
        hdu.header["CRPIX2"] -= y_min

        hdu_out = fits.PrimaryHDU(data=cutout_data, header=hdu.header)
        hdu_out.writeto(outfile)

    return hdu_out


def make_cutout_region(image: str, region: str, outfile: str):
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

    assert len(r) == 1, "Multiple regions in 1 file given, only one allowed"

    shape = array(r[0].coord_list)

    # circle
    if len(shape) == 3:
        cutout = Cutout2D(
            data=data * mask,
            position=(shape[0], shape[1]),
            size=(shape[2], shape[2]),
            wcs=WCS(head, naxis=2),
            mode="partial",
        )
    # square
    elif len(shape) > 3:
        print(shape)

        cutout = Cutout2D(
            data=data * mask,
            position=(shape[0], shape[1]),
            size=(shape[3], shape[2]),
            wcs=WCS(head, naxis=2),
            mode="partial",
        )
    else:
        raise RuntimeError("Non-square or circular region found.")

    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(outfile)
