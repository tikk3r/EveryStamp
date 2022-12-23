#!/usr/bin/env python
''' Python library aiming to provide a wrapper around various astronomical surveys that offer cutouts.'''
__version__ = 'v0.0.0'
__author__ = 'Frits Sweijen'
__license__ = 'GPLv3'

import logging
logging.basicConfig(format='[%(name)s] %(asctime)s - %(levelname)s: %(message)s', level=logging.INFO)
logger = logging.getLogger('EveryStamp')


def main():
    '''Main entry point if called as a standalone executable.
    '''
    import argparse
    parser = argparse.ArgumentParser(description='EveryStamps {:s} by {:s}'.format(__version__, __author__), add_help=True, usage=argparse.SUPPRESS)
    parser._action_groups.pop()

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('--survey', type=str, required=True, choices=['legacy', 'pan-starrs', 'vlass'], help='Survey from which to download the cutout.')
    required_args.add_argument('--ra', type=float, required=True, help='Right ascension of cutout centre in degrees.')
    required_args.add_argument('--dec', type=float, required=True, help='Declination of cutout centre in degrees.')
    required_args.add_argument('--size', type=float, required=False, default=0.01, help='Cutout size in degrees.')

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('--download_dir', type=str, required=False, default='', dest='ddir', help='Directory to store downloaded files. If not given will download to $PWD.')
    optional_args.add_argument('--mode', type=str, required=False, default='jpeg', choices=['jpeg', 'fits', 'both'], help='Image type to retrieve. Can be "jpeg", "fits" or "both" to retrieve either a JPEG image, FITS file or both. Default value is jpeg.')

    legacy_args = parser.add_argument_group('[DESI Legacy Imaging Surveys]')
    legacy_args.add_argument('--legacy_bands', type=str, required=False, help='Bands to download. Allowed values are g, r and z. Multiple bands can be specified as a single string. In the case of a JPEG image a colour image will be generated. In the case of a FITS image a FITS cube will be downloaded. Default: grz')
    legacy_args.add_argument('--legacy_layer', type=str, required=False, default='ls-dr9', help='Layer to make a cutout from. Default value is ls-dr9. Examples are ls-dr9, sdss or unwise-neo4. See Legacy documentation for all possibilies.')
    legacy_args.add_argument('--legacy_autoscale', required=False, default=False, action='store_true', help='Automatically change the pixel size if the resulting image would exceed the server maximum of 3000x3000 pixels.')

    ps_args = parser.add_argument_group('[Pan-STARRS]')
    ps_args.add_argument('--ps_bands', type=str, required=False, default='gri', help='Bands to download. Allowed values are g, r and i. Multiple bands can be specified as a single string. Default: gri')

    vlass_args = parser.add_argument_group('[VLASS]')
    vlass_args.add_argument('--vlass_ms', type=str, required=False, default='', help='Measurement Set to take the cutout position from.')
    vlass_args.add_argument('--vlass_consider_QA_rejected', type=bool, required=False, default=False, help='Also consider tiles that failed the Quality Assurance checks.')

    args = parser.parse_args()
    logger.info('Survey is %s', args.survey)
    if args.survey == 'legacy':
        from everystamp.downloaders import LegacyDownloader
        ld = LegacyDownloader()
        ld.download(ra=args.ra, dec=args.dec, bands=args.legacy_bands, mode=args.mode, size=args.size, layer=args.legacy_layer, autoscale=args.legacy_autoscale, ddir=args.ddir)
    elif args.survey == 'panstarrs':
        from everystamp.downloaders import PanSTARRSDownloader
        pd = PanSTARRSDownloader()
        pd.download(ra=args.ra, dec=args.dec, bands=args.ps_bands, mode=args.mode, size=args.size, ddir=args.ddir)
    elif args.survey == 'vlass':
        from everystamp.downloaders import VLASSDownloader
        vd = VLASSDownloader()
        vd.download(ra=args.ra, dec=args.dec, crop=True, consider_QA_rejected=args.vlass_consider_QA_rejected, ddir=args.ddir)


if __name__ == '__main__':
    main()
