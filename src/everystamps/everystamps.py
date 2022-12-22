#!/usr/bin/env python
''' Python library aiming to provide a wrapper around various astronomical surveys that offer cutouts.'''
__version__ = 'v0.0.0'
__author__ = 'Frits Sweijen'
__license__ = 'GPLv3'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Legacystamps {:s} by {:s}'.format(__version__, __author__))
    parser.add_argument('--ra', type=float, required=True, help='Right ascension of cutout centre in degrees.')
    parser.add_argument('--dec', type=float, required=True, help='Declination of cutout centre in degrees.')
    parser.add_argument('--bands', type=str, required=True, help='Bands to download. Allowed values are g, r and z. Multiple bands can be specified as a single string. In the case of a JPEG image a colour image will be generated. In the case of a FITS image a FITS cube will be downloaded.')
    parser.add_argument('--mode', type=str, required=False, default='jpeg', help='Image type to retrieve. Can be "jpeg", "fits" or "both" to retrieve either a JPEG image, FITS file or both. Default value is jpeg.')
    parser.add_argument('--size', type=float, required=False, default=0.01, help='Cutout size in degrees.')
    parser.add_argument('--layer', type=str, required=False, default='ls-dr9', help='Layer to make a cutout from. Default value is ls-dr9. Examples are ls-dr9, sdss or unwise-neo4. See Legacy documentation for all possibilies.')
    parser.add_argument('--autoscale', required=False, default=False, dest='autoscale', action='store_true', help='Automatically change the pixel size if the resulting image would exceed the server maximum of 3000x3000 pixels.')
    parser.add_argument('--download-dir', type=str, required=False, default='', dest='ddir', help='Directory to store downloaded files. If not given will download to $PWD.')
    args = parser.parse_args()
    download(args.ra, args.dec, args.bands, mode=args.mode, size=args.size, layer=args.layer, autoscale=args.autoscale, ddir=args.ddir)