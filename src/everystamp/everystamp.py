#!/usr/bin/env python
''' Python library aiming to provide a wrapper around various astronomical surveys that offer cutouts.'''
__version__ = '1.0.0'
__author__ = 'Frits Sweijen'
__license__ = 'GPLv3'

import argparse
import logging
logging.basicConfig(format='[%(name)s] %(asctime)s - %(levelname)s: %(message)s', level=logging.INFO)
logger = logging.getLogger('EveryStamp')

from astroquery.skyview import SkyView
from collections.abc import Iterable
import requests

def flatten(xs): 
    for x in xs: 
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)): 
            yield from flatten(x) 
        else: 
            yield x 


def add_args_download(parser):
    custom_surveys = ['legacy', 'pan-starrs', 'vlass', 'lolss', 'lotss', 'tgss']
    try:
        skyview_surveys = list(flatten(list(SkyView.survey_dict.values())))
    except requests.exceptions.ConnectionError:
        logger.warning('Failed to get SkyView surveys. SkyView cutouts will not be available.')
        skyview_surveys = []
    allowed_surveys = custom_surveys + skyview_surveys
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('--survey', type=str, required=True, choices=allowed_surveys, help='Survey from which to download the cutout.')
    required_args.add_argument('--ra', type=float, required=True, help='Right ascension of cutout centre in degrees.')
    required_args.add_argument('--dec', type=float, required=True, help='Declination of cutout centre in degrees.')
    required_args.add_argument('--mode', type=str, required=True, default='jpeg', choices=['jpeg', 'fits', 'both'], help='Image type to retrieve. Can be "jpeg", "fits" or "both" to retrieve either a JPEG image, FITS file or both. Default value is jpeg.')

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('--download_dir', type=str, required=False, default='', dest='ddir', help='Directory to store downloaded files. If not given will download to $PWD.')
    optional_args.add_argument('--size', type=float, required=False, default=0.01, help='Cutout size in degrees.')

    legacy_args = parser.add_argument_group('[DESI Legacy Imaging Surveys]')
    legacy_args.add_argument('--legacy_bands', type=str, required=False, help='Bands to download. Allowed values are g, r and z. Multiple bands can be specified as a single string. In the case of a JPEG image a colour image will be generated. In the case of a FITS image a FITS cube will be downloaded. Default: grz')
    legacy_args.add_argument('--legacy_layer', type=str, required=False, default='ls-dr9', help='Layer to make a cutout from. Default value is ls-dr9. Examples are ls-dr9, sdss or unwise-neo4. See Legacy documentation for all possibilies.')
    legacy_args.add_argument('--legacy_autoscale', required=False, default=False, action='store_true', help='Automatically change the pixel size if the resulting image would exceed the server maximum of 3000x3000 pixels.')

    ps_args = parser.add_argument_group('[Pan-STARRS]')
    ps_args.add_argument('--ps_bands', type=str, required=False, default='gri', help='Bands to download. Allowed values are g, r and i. Multiple bands can be specified as a single string. Default: gri')

    vlass_args = parser.add_argument_group('[VLASS]')
    vlass_args.add_argument('--vlass_ms', type=str, required=False, default='', help='Measurement Set to take the cutout position from.')
    vlass_args.add_argument('--vlass_consider_QA_rejected', type=bool, required=False, default=False, help='Also consider tiles that failed the Quality Assurance checks.')

    lolss_args = parser.add_argument_group('[LoLSS]')
    lolss_args.add_argument('--lolss_release', type=str, required=False, default='pdr', choices=['pdr'], help='Data release to download from.')

    lotss_args = parser.add_argument_group('[LoTSS]')
    lotss_args.add_argument('--lotss_release', type=str, required=False, default='dr1', choices=['pdr', 'dr1', 'dr2'], help='Data release to download from.')


def add_args_plot(parser):
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('--image', type=str, required=False, help='FITS image to plot.')

    required_args = parser.add_argument_group('Optional arguments')
    required_args.add_argument('--gamma', type=float, default=1.0, required=False, help='Gamma compress (<1) or expand (>1) an image.')
    required_args.add_argument('--CLAHE', action='store_true', default=False, required=False, help='Apply contrast-limited adaptive histogram equalisation.')
    required_args.add_argument('--CLAHE-gridsize', default=5, type=int, required=False, help='Grid size to use for CLAHE.')
    required_args.add_argument('--CLAHE-cliplim', default=1.0, type=float, required=False, help='Clip limit to use for CLAHE.')


def process_args_download(args):
    logger.info('Survey is %s', args.survey)
    if args.survey == 'legacy':
        from everystamp.downloaders import LegacyDownloader
        ld = LegacyDownloader()
        ld.download(ra=args.ra, dec=args.dec, bands=args.legacy_bands, mode=args.mode, size=args.size, layer=args.legacy_layer, autoscale=args.legacy_autoscale, ddir=args.ddir)
    elif args.survey == 'pan-starrs':
        from everystamp.downloaders import PanSTARRSDownloader
        pd = PanSTARRSDownloader()
        pd.download(ra=args.ra, dec=args.dec, bands=args.ps_bands, mode=args.mode, size=args.size, ddir=args.ddir)
    elif args.survey == 'vlass':
        if args.mode == 'both' or args.mode == 'jpeg':
            raise ValueError('VLASS download does not support JPEG (yet).')
        from everystamp.downloaders import VLASSDownloader
        vd = VLASSDownloader()
        vd.download(ra=args.ra, dec=args.dec, crop=True, consider_QA_rejected=args.vlass_consider_QA_rejected, ddir=args.ddir)
    elif args.survey == 'lolss':
        if args.mode == 'both' or args.mode == 'jpeg':
            raise ValueError('LoLLS download does not support JPEG (yet).')
        from everystamp.downloaders import VODownloader
        vd = VODownloader(url='https://vo.astron.nl/lolss/q/cutout/siap.xml', name='LoLSS')
        vd.download(ra=args.ra, dec=args.dec, size=args.size, ddir=args.ddir)
    elif args.survey == 'lotss':
        if (args.mode == 'both' or args.mode == 'jpeg') and (args.lotss_release != 'dr2'):
            raise ValueError('LoTSS {:s} download does not support JPEG (yet).'.format(args.lotss_release.upper()))
        from everystamp.downloaders import VODownloader
        if args.lotss_release == 'pdr':
            vd = VODownloader(url='https://vo.astron.nl/lofartier1/q_img/cutout/siap.xml', name='LoTSS-PDR')
        elif args.lotss_release == 'dr1':
            vd = VODownloader(url='https://vo.astron.nl/hetdex/lotss-dr1-img/cutout/siap.xml', name='LoTSS-DR1')
            vd.download(ra=args.ra, dec=args.dec, size=args.size, ddir=args.ddir)
        elif args.lotss_release == 'dr2':
            from everystamp.downloaders import HiPSDownloader
            vd = HiPSDownloader(hips='astron.nl/P/lotss_dr2_high', name='LoTSS-DR2')
            vd.download(ra=args.ra, dec=args.dec, size=args.size, ddir=args.ddir, mode=args.mode.replace('e', ''), pixsize=1.5)
    elif args.survey == 'tgss':
        if args.mode == 'both' or args.mode == 'jpeg':
            raise ValueError('TGSS download does not support JPEG (yet).')
        from everystamp.downloaders import VODownloader
        vd = VODownloader(url='https://vo.astron.nl/tgssadr/q_fits/cutout/siap.xml', name='TGSS')
        vd.download(ra=args.ra, dec=args.dec, size=args.size, ddir=args.ddir)
    else:
        if args.mode == 'both' or args.mode == 'jpeg':
            raise ValueError('SkyView download does not support JPEG (yet).')
        from everystamp.downloaders import SkyViewDownloader
        sd = SkyViewDownloader(args.survey)
        sd.download(ra=args.ra, dec=args.dec, size=args.size, ddir=args.ddir)

def process_args_plot(args):
    logger.info('Plotting image %s', args.image)
    from everystamp.plotters import BasicPlot
    from everystamp.tonemapping import gamma, make_nonnegative
    import numpy as np
    bp = BasicPlot(args.image)
    if args.CLAHE:
        import cv2
        bp.data = make_nonnegative(bp.fitsdata)
        bp.data /= np.nanmax(bp.data)
        bp.data *= 2**16
        bp.data = bp.data.astype(np.uint16)
        clahe = cv2.createCLAHE(clipLimit=args.CLAHE_cliplim, tileGridSize=(args.CLAHE_gridsize, args.CLAHE_gridsize))
        bp.data = clahe.apply(bp.data)
    if args.gamma:
        bp.data = gamma(bp.data, args.gamma)
    bp.plot2D()


def main():
    '''Main entry point if called as a standalone executable.
    '''
    parser = argparse.ArgumentParser(description='EveryStamp {:s} by {:s}'.format(__version__, __author__))
    parser._action_groups.pop()

    subparsers = parser.add_subparsers(dest='cmd', description='Description of sub commands.')
    subparser_dl = subparsers.add_parser('download', usage='everystamp download --survey SURVEY --ra RA --dec DEC --mode MODE --size SIZE', description='Download a cutout from a user-specified survey. See everystamp download -h for all available parameters.', help='Download a cutout from a specified survey.')
    add_args_download(subparser_dl)

    subparser_plot = subparsers.add_parser('plot', description='Plot a given FITS image. See everystamp plot -h for more information.', help='Plot a user-supplied FITS image.')
    add_args_plot(subparser_plot)
    
    args = parser.parse_args()

    if args.cmd == 'download':
        process_args_download(args)
    if args.cmd == 'plot':
        process_args_plot(args)

if __name__ == '__main__':
    main()