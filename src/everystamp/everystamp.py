#!/usr/bin/env python
''' Python library aiming to provide a wrapper around various astronomical surveys that offer cutouts.'''
__version__ = '1.3.0'
__author__ = 'Frits Sweijen'
__license__ = 'GPLv3'
from typing import Generator

import argparse
import logging
logging.basicConfig(format='[%(name)s] %(asctime)s - %(levelname)s: %(message)s', level=logging.INFO)
logger = logging.getLogger('EveryStamp')
import os
import shutil
import subprocess

from astroquery.skyview import SkyView #type: ignore
from collections.abc import Iterable
from everystamp.tonemapping import lhdr
import astropy.visualization
import cv2 #type: ignore
import requests

# Check if LuminanceHDR is installed.
HAS_LHDR = lhdr.has_luminance_hdr() #True if not subprocess.getstatusoutput('luminance-hdr-cli') else False
#HAS_LHDR = shutil.which('luminance-hdr-cli')


def flatten(xs: Iterable) -> Generator:
    ''' Generator to flatten a list of nested lists.

    Args:
        xs : list
            The list to flatten.

    Yields:
        x : iterable
            An iterable that will generate the flattened list.
    '''
    for x in xs: 
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)): 
            yield from flatten(x) 
        else: 
            yield x 


def _add_args_download(parser):
    ''' Add arguments to the argparse instance for downloading cutouts.
    
    Args:
        parser : ArgumentParser
            ArgumentParser instance to which to add entries.
    '''
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
    optional_args.add_argument('--download_dir', type=str, required=False, default=os.getcwd(), dest='ddir', help='Directory to store downloaded files. If not given will download to $PWD.')
    optional_args.add_argument('--size', type=float, required=False, default=0.01, help='Cutout size in degrees.')
    optional_args.add_argument('--skyview_pixsize', type=float, required=False, default=1.0, help='Pixel size in arcsec for SkyView cutouts.')

    legacy_args = parser.add_argument_group('[DESI Legacy Imaging Surveys]')
    legacy_args.add_argument('--legacy_bands', default='grz', type=str, required=False, help='Bands to download. Allowed values are g, r and z. Multiple bands can be specified as a single string. In the case of a JPEG image a colour image will be generated. In the case of a FITS image a FITS cube will be downloaded. Default: grz')
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


def _add_args_plot(parser):
    ''' Add arguments to the argparse instance for plotting images.
    
    Args:
        parser : ArgumentParser
            ArgumentParser instance to which to add entries.
    '''
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('--image', type=str, required=False, help='FITS image to plot.')

    required_args = parser.add_argument_group('Optional arguments')
    required_args.add_argument('--gamma', type=float, default=1.0, required=False, help='Gamma compress (<1) or expand (>1) an image after tone mapping.')
    required_args.add_argument('--CLAHE', action='store_true', default=False, required=False, help='Apply contrast-limited adaptive histogram equalisation.')
    required_args.add_argument('--CLAHE-gridsize', default=5, type=int, required=False, help='Grid size to use for CLAHE.')
    required_args.add_argument('--CLAHE-cliplim', default=1.0, type=float, required=False, help='Clip limit to use for CLAHE.')
    required_args.add_argument('--stretch', default=None, type=str, required=False, choices=['log', 'sqrt', 'squared', 'asinh', 'sinh'], help='Stretch an image with a certian function.')

    required_args.add_argument('--contour_image', default=None, required=False, help='Plot the given image as contours over the main image.')

    if HAS_LHDR:
        required_args.add_argument('--hdr-tonemap', default=None, type=str, choices=['ashikhmin', 'drago', 'duran', 'fattal', 'ferradans', 'ferwerda', 'kimkautz', 'lischinski', 'mantiuk06', 'mantiuk08', 'pattanaik', 'reinhard02', 'reinhard05', 'vanhateren'], required=False, help='HDR tonemapping to apply')

        hdr_ashikhmin_args = parser.add_argument_group('HDR Tone mapping -- Ashikhmin et al. 2002 arguments')
        hdr_ashikhmin_args.add_argument('--ashikhmin-eq2', default=True, type=bool, required=False, help='Equation 2?')
        hdr_ashikhmin_args.add_argument('--ashikhmin-simple', default=True, type=bool, required=False, help='Simple?')
        hdr_ashikhmin_args.add_argument('--ashikhmin-local_threshold', default=None, type=float, required=False, help='Local threshold.')

        hdr_drago_args = parser.add_argument_group('HDR Tone mapping -- Drago et al. 2003 arguments')
        hdr_drago_args.add_argument('--drago-bias', default=0.85, type=float, required=False, help='Bias parameter controlling the exponent base.')

        hdr_fattal_args = parser.add_argument_group('HDR Tone mapping -- Fattal et al. 2002 arguments')
        hdr_fattal_args.add_argument('--fattal-alpha', default=None, type=float, required=False, help='Controls which gradient magnitude is preserved.')
        hdr_fattal_args.add_argument('--fattal-beta', default=None, type=float, required=False, help='Controls local detail enhancement.')
        hdr_fattal_args.add_argument('--fattal-colour_saturation', default=None, type=float, required=False, help='Controls colour saturation.')
        hdr_fattal_args.add_argument('--fattal-noise', default=None, type=float, required=False, help='Controls when local detail enhancement is reduced.')

        hdr_ferradans_args = parser.add_argument_group('HDR Tone mapping -- Ferradans et al. 2011 arguments')
        hdr_ferradans_args.add_argument('--ferradans-rho', default=-2, type=float, required=False, help='Controls overall lightness.')
        hdr_ferradans_args.add_argument('--ferradans-inverse_alpha', default=5, type=float, required=False, help='Controls local detail enhancement.')

        hdr_ferwerda_args = parser.add_argument_group('HDR Tone mapping -- Ferwerda et al. 1996 arguments')
        hdr_ferwerda_args.add_argument('--ferwerda-multiplier', default=None, type=float, required=False, help='Ferwerda multiplier factor.')
        hdr_ferwerda_args.add_argument('--ferwerda-luminance_adaptation', default=None, type=float, required=False, help='Ferwerda adaptive luminance factor.')

        hdr_hateren_args = parser.add_argument_group('HDR Tone mapping -- van Hateren 2006 arguments')
        hdr_hateren_args.add_argument('--vanhateren-pupil_area', default=None, type=float, required=False, help='Pupil area.')

        hdr_kimkautz_args = parser.add_argument_group('HDR Tone mapping -- Kim and Kautz 2008 arguments')
        hdr_kimkautz_args.add_argument('--kimkautz-c1', default=None, type=float, required=False, help='Kim and Kautz c1 factor.')
        hdr_kimkautz_args.add_argument('--kimkautz-c2', default=None, type=float, required=False, help='Kim and Kautz c2 factor.')

        hdr_lischinski_args = parser.add_argument_group('HDR Tone mapping -- Lischinski 2006 arguments')
        hdr_lischinski_args.add_argument('--lischinski-alpha', default=None, type=float, required=False, help='Alpha.')

        hdr_mantiuk06_args = parser.add_argument_group('HDR Tone mapping -- Mantiuk et al. 2006 arguments')
        hdr_mantiuk06_args.add_argument('--mantiuk06-contrast', default=None, type=float, required=False, help='Contrast factor.')
        hdr_mantiuk06_args.add_argument('--mantiuk06-saturation', default=None, type=float, required=False, help='Saturation factor.')
        hdr_mantiuk06_args.add_argument('--mantiuk06-detail', default=None, type=float, required=False, help='Detail factor.')
        hdr_mantiuk06_args.add_argument('--mantiuk06-contrast_equalisation', default=True, type=bool, required=False, help='Equalise contrast?')

        hdr_mantiuk08_args = parser.add_argument_group('HDR Tone mapping -- Mantiuk et al. 2008 arguments')
        hdr_mantiuk08_args.add_argument('--mantiuk08-contrast_enhancement', default=None, type=float, required=False, help='Contrast enhancement factor.')
        hdr_mantiuk08_args.add_argument('--mantiuk08-colour_saturation', default=None, type=float, required=False, help='Colour saturation factor.')
        hdr_mantiuk08_args.add_argument('--mantiuk08-luminance_level', default=0.0, type=float, required=False, help='Luminance level.')
        hdr_mantiuk08_args.add_argument('--mantiuk08-set_luminance', default=True, type=bool, required=False, help='Set luminance level?')

        hdr_duran_args = parser.add_argument_group('HDR Tone mapping -- Durand and Dorsey et al. 2002 arguments')
        hdr_duran_args.add_argument('--duran-sigma_spatial', default=None, type=float, required=False, help='Spatial kernel size.')
        hdr_duran_args.add_argument('--duran-sigma_range', default=None, type=float, required=False, help='Range kernel size.')
        hdr_duran_args.add_argument('--duran-base_contrast', default=None, type=float, required=False, help='Base contrast.')

        hdr_reinhard02_args = parser.add_argument_group('HDR Tone mapping -- Reinhard et al. 2002 arguments')
        hdr_reinhard02_args.add_argument('--reinhard02-key', default=None, type=float, required=False, help='Key.')
        hdr_reinhard02_args.add_argument('--reinhard02-phi', default=None, type=float, required=False, help='Phi.')
        hdr_reinhard02_args.add_argument('--reinhard02-use_scales', default=None, type=float, required=False, help='Use scales?')
        hdr_reinhard02_args.add_argument('--reinhard02-range', default=None, type=float, required=False, help='Range.')
        hdr_reinhard02_args.add_argument('--reinhard02-low', default=None, type=float, required=False, help='Low.')
        hdr_reinhard02_args.add_argument('--reinhard02-high', default=None, type=float, required=False, help='High.')

        hdr_reinhard05_args = parser.add_argument_group('HDR Tone mapping -- Reinhard et al. 2005 arguments')
        hdr_reinhard05_args.add_argument('--reinhard05-brightness', default=None, type=float, required=False, help='Brightness.')
        hdr_reinhard05_args.add_argument('--reinhard05-chroma', default=None, type=float, required=False, help='Chroma.')
        hdr_reinhard05_args.add_argument('--reinhard05-lightness', default=None, type=float, required=False, help='Lightness.')

        hdr_pattanaik_args = parser.add_argument_group('HDR Tone mapping -- Pattanaik arguments')
        hdr_pattanaik_args.add_argument('--pattanaik-multiplier', default=None, type=float, required=False, help='Multiplier.')
        hdr_pattanaik_args.add_argument('--pattanaik-local_tonemap', default=True, type=bool, required=False, help='Multiplier.')
        hdr_pattanaik_args.add_argument('--pattanaik-auto_lum', default=True, type=bool, required=False, help='Automatic luminance?')
        hdr_pattanaik_args.add_argument('--pattanaik-cone_level', default=None, type=float, required=False, help='Cone level.')
        hdr_pattanaik_args.add_argument('--pattanaik-rod_level', default=None, type=float, required=False, help='Rod level.')
    else:
        logger.warning('Cannot find luminance-hdr-cli. HDR tone mapping functionality will not be available unless LuminanceHDR is (correctly) installed.')


def _process_args_download(args):
    ''' Process arguments to the argparse instance for downloading cutouts.
    
    Args:
        parser : ArgumentParser
            ArgumentParser instance to which to add entries.
    '''
    if args.ddir and (not os.path.exists(args.ddir)):
        logger.info('Download directory does not exist, creating it')
        os.mkdir(args.ddir)
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
        sd.download(ra=args.ra, dec=args.dec, size=args.size, pixsize=args.skyview_pixsize, ddir=args.ddir)


def _process_args_plot(args):
    ''' Process arguments to the argparse instance for plotting images.
    
    Args:
        parser : ArgumentParser
            ArgumentParser instance to which to add entries.
    '''
    logger.info('Plotting image %s', args.image)
    from everystamp.plotters import BasicFITSPlot, BasicImagePlot
    from everystamp.tonemapping import gamma, make_nonnegative
    import numpy as np
    from everystamp.tonemapping.lhdr import ashikhmin, drago, duran, fattal, ferradans, ferwerda, kimkautz, lischinski, mantiuk06, mantiuk08, reinhard02, reinhard05, pattanaik, vanhateren
    if args.image.lower().endswith('fits'):
        bp = BasicFITSPlot(args.image)
    else:
        # Probably an image format.
        bp = BasicImagePlot(args.image)
    if HAS_LHDR and args.hdr_tonemap:
        if args.hdr_tonemap == 'ashikhmin':
            logger.info('Tonemapping image with ashikhmin')
            bp.data = ashikhmin(bp.data, eq2=args.ashikhmin_eq2, simple=args.ashikhmin_simple, local_threshold=args.ashikhmin_local_threshold)
        if args.hdr_tonemap == 'fattal':
            logger.info('Tonemapping image with fattal')
            bp.data = fattal(bp.data, alpha=args.fattal_alpha, beta=args.fattal_beta, colour_saturation=args.fattal_colour_saturation, noise=args.fattal_noise)
        if args.hdr_tonemap == 'drago':
            logger.info('Tonemapping image with drago')
            bp.data = drago(bp.data, bias=args.drago_bias)
        if args.hdr_tonemap == 'ferradans':
            logger.info('Tonemapping image with ferradans')
            bp.data = ferradans(bp.data, rho=args.ferradans_rho, inv_alpha=args.ferradans_inverse_alpha)
        if args.hdr_tonemap == 'ferwerda':
            logger.info('Tonemapping image with ferwerda')
            bp.data = ferwerda(bp.data, multiplier=args.ferwerda_multiplier, luminance_adaptation=args.ferwerda_luminance_adaptation)
        if args.hdr_tonemap == 'kimkautz':
            logger.info('Tonemapping image with kimkautz')
            bp.data = kimkautz(bp.data, c1=args.kimkautz_c1, c2=args.kimkautz_c2)
        if args.hdr_tonemap == 'lischinski':
            logger.info('Tonemapping with lischinski')
            bp.data = lischinski(bp.data, alpha=args.lischinski_alpha)
        if args.hdr_tonemap == 'mantiuk06':
            logger.info('Tonemapping image with mantiuk06')
            bp.data = mantiuk06(bp.data, contrast=args.mantiuk06_contrast, saturation=args.mantiuk06_saturation, detail=args.mantiuk06_detail, contrast_equalisation=args.mantiuk06_contrast_equalisation)
        if args.hdr_tonemap == 'mantiuk08':
            logger.info('Tonemapping image with mantiuk08')
            bp.data = mantiuk08(bp.data, contrast_enhancement=args.mantiuk08_colour_saturation, colour_saturation=args.mantiuk08_colour_saturation, luminance_level=args.mantiuk08_luminance_level, set_luminance=args.mantiuk08_set_luminance)
        if args.hdr_tonemap == 'duran':
            logger.info('Tonemapping image with duran')
            bp.data = duran(bp.data, sigma_spatial=args.duran_sigma_spatial, sigma_range=args.duran_sigma_range, base_contrast=args.duran_base_contrast)
        if args.hdr_tonemap == 'pattanaik':
            logger.info('Tonemapping image with pattanaik')
            bp.data = pattanaik(bp.data, multiplier=args.pattanaik_multiplier, local_tonemap=args.pattanaik_local_tonemap, auto_lum=args.pattanaik_auto_lum, cone_level=args.pattanaik_cone_level, rod_level=args.pattanaik_rod_level)
        if args.hdr_tonemap == 'reinhard02':
            logger.info('Tonemapping image with reinhard02')
            bp.data = reinhard02(bp.data, key=args.reinhard02_key, phi=args.reinhard02_phi, use_scales=args.reinhard02_use_scales, range=args.reinhard02_range, low=args.reinhard02_low, high=args.reinhard02_high)
        if args.hdr_tonemap == 'reinhard05':
            logger.info('Tonemapping image with reinhard05')
            bp.data = reinhard05(bp.data, brightness=args.reinhard05_brightness, chroma=args.reinhard05_chroma, lightness=args.reinhard05_lightness)
        if args.hdr_tonemap == 'vanhateren':
            logger.info('Tonemapping with vanhateren')
            bp.data = vanhateren(bp.data, pupil_area=args.vanhateren_pupil_area)
    
    if args.CLAHE:
        if args.image.lower().endswith('fits'):
            bp.data = make_nonnegative(bp.fitsdata)
            bp.data /= np.nanmax(bp.data)
            bp.data *= 2**16
            bp.data = bp.data.astype(np.uint16)
            clahe = cv2.createCLAHE(clipLimit=args.CLAHE_cliplim, tileGridSize=(args.CLAHE_gridsize, args.CLAHE_gridsize))
            bp.data = clahe.apply(bp.data)
        else:
            clahe = cv2.createCLAHE(clipLimit=args.CLAHE_cliplim, tileGridSize=(args.CLAHE_gridsize, args.CLAHE_gridsize))
            Lab = cv2.cvtColor(bp.data, cv2.COLOR_RGB2Lab)
            L, a, b = cv2.split(Lab)
            L_clahe = clahe.apply(L)
            Lab_clahe = cv2.merge((L_clahe, a, b))
            cv2.cvtColor(Lab_clahe, cv2.COLOR_Lab2RGB)

    if args.gamma != 1:
        logger.info('Applying gamma stretch of {:f}'.format(args.gamma))
        if args.image.lower().endswith('fits'):
            bp.data = gamma(bp.data / np.nanmax(bp.data), args.gamma)
        else:
            Lab = cv2.cvtColor(bp.data.astype(np.uint8), cv2.COLOR_RGB2Lab)
            L, a, b = cv2.split(Lab)
            lookUpTable = np.empty((1,256), np.uint8)
            for i in range(256):
                lookUpTable[0,i] = np.clip(pow(i / 255.0, args.gamma) * 255.0, 0, 255)
            L_gamma = cv2.LUT(L, lookUpTable)
            Lab_gamma = cv2.merge((L_gamma, a, b))
            bp.data = cv2.cvtColor(Lab_gamma, cv2.COLOR_Lab2RGB)
    
    if args.stretch:
        # 'log10', 'sqrt', 'squared', 'asinh', 'sinh'
        if args.stretch == 'log':
            stretch = astropy.visualization.LogStretch()
        elif args.stretch == 'sqrt':
            stretch = astropy.visualization.SqrtStretch()
        elif args.stretch == 'squared':
            stretch = astropy.visualization.SquaredStretch()
        elif args.stretch == 'asinh':
            stretch = astropy.visualization.AsinhStretch()
        elif args.stretch == 'sinh':
            stretch = astropy.visualization.SinhStretch()
        bp.data = stretch(bp.data / np.nanmax(bp.data), min)
    
    if args.contour_image:
        bp.plot2D(contour_image = args.contour_image)
    else:
        bp.plot2D()
    if args.image.lower().endswith('fits'):
        bp.savedata(args.image.replace('.fits', '.tonemapped.fits'))
        bp.plot_noaxes()


def main():
    '''Main entry point if called as a standalone executable.
    '''
    parser = argparse.ArgumentParser(description='EveryStamp {:s} by {:s}'.format(__version__, __author__))
    parser._action_groups.pop()

    subparsers = parser.add_subparsers(dest='cmd', description='Description of sub commands.')
    subparser_dl = subparsers.add_parser('download', usage='everystamp download --survey SURVEY --ra RA --dec DEC --mode MODE --size SIZE', description='Download a cutout from a user-specified survey. See everystamp download -h for all available parameters.', help='Download a cutout from a specified survey.')
    _add_args_download(subparser_dl)

    subparser_plot = subparsers.add_parser('plot', description='Plot a given FITS image. See everystamp plot -h for more information.', help='Plot a user-supplied FITS image.')
    _add_args_plot(subparser_plot)

    args = parser.parse_args()

    if args.cmd == 'download':
        _process_args_download(args)
    if args.cmd == 'plot':
        _process_args_plot(args)


if __name__ == '__main__':
    main()
