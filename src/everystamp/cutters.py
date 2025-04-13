"""Sub-module for trimming FITS images."""

import numpy as np
import pyregion

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import Header
from astropy.nddata import Cutout2D
from astropy.units import Quantity
from astropy.wcs import WCS
from numpy import array

def copy_header_keywords(in_header: Header, out_header: Header):
    for key in in_header.keys():
        if key not in out_header.keys():
            out_header[key] = in_header[key]

def make_cutout_2D(
    image: str,
    pos: SkyCoord,
    size: Quantity,
    outfile: str,
):
    """Make a square cutout of given size using Astropy's Cutout2D function.

    Args:
        image: name of the FITS image to make a cutout of.
        pos: position on which to centre the cutout.
        size: width of the cutout.
        outfile: name of the trimmed output image.
    """
    head = fits.getheader(image)
    data = fits.getdata(image).squeeze()
    if len(data.shape) > 2:
        raise ValueError(
            "Image data has dimensions other than RA and DEC, which is not yet supported."
        )
    wcs = WCS(head).celestial

    cutout = Cutout2D(data, pos, size, wcs=wcs)
    head_new = cutout.wcs.to_header()
    copy_header_keywords(head, head_new)
    hdu = fits.PrimaryHDU(data=cutout.data, header=head_new)
    hdu.writeto(outfile)


def make_cutout_2D_fast(
    image: str,
    pos: SkyCoord,
    size: Quantity,
    outfile: str,
):
    """Make a square cutout of given size by slicing instead of Astropy's Cutout2D function.

    Args:
        image: name of the FITS image to make a cutout of.
        pos: position on which to centre the cutout.
        size: width of the cutout.
        outfile: name of the trimmed output image.
    """
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


def make_cutout_region(image: str, region: str, outfile: str):
    """Make a cutout from a FITS image using a DS9 region file.

    Args:
        image: name of the FITS image to make a cutout of.
        region: the region to which the image will be trimmed, in DS9 format.
        outfile: name of the trimmed output image.

    Raises:
        NotImplementedError: for cubes with 3 or more dimensions.
        RuntimeError: if the region cannot be parsed for whatever reason.
    """
    hdu = fits.open(image)
    head = hdu[0].header
    data = hdu[0].data.squeeze()

    if data.ndim > 2:
        raise NotImplementedError("Region cutouts for >2D images not supported.")

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

    head_new = cutout.wcs.to_header()
    copy_header_keywords(head, head_new)
    hdu = fits.PrimaryHDU(data=cutout.data, header=head_new)
    hdu.writeto(outfile)
