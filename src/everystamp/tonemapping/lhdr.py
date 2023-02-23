'''Sub-module containing wrappers to LuminosityHDR's HDR tonemapping programs.'''
from typing import Union

import os
import subprocess

from astropy.io import fits
from skimage import io
import numpy

BASECOMMAND = 'luminance-hdr-cli'


def run_command(cmd: Union[list, str]) -> None:
    ''' Run a command through subprocess.
    
    Args:
        cmd : list
            List of arguments building up the command or a string of the full command.
    
    Returns:
        None
    '''
    print('Running command:', ' '.join(cmd))
    try:
        if type(cmd) is list:
            subprocess.run(cmd)
        elif type(cmd) is str:
            subprocess.run(cmd.split(' '))
        else:
            raise ValueError('Malformed command: ', cmd)
    except subprocess.CalledProcessError as e:
        print("Tone mapping failed with error code", e.returncode)
        print(e.output)


def _load_tonemapped_tmpdata(name: str, as_gray: bool = True) -> numpy.ndarray:
    ''' Load the tonemapped output of LuminanceHDR ad return it as an array.
    
    Args:
        name : str
            File path to the data to load.
        as_gray : bool
            Treat the image as gray scale or not. Default: True.
            
    Returns:
        data : numpy.ndarray
            NumPy array containing the image.
    '''
    return io.imread(name, as_gray=as_gray)


def _store_tmpfile(data: numpy.ndarray, name: str) -> str:
    ''' Store the data to tonemap to a temporary FITS file to pass to LuminanceHDR.
    
    Args:
        data : numpy.ndarray
            NumPy array containing the image.

    Returns:
        path : str
            Absolute path to the temporary file.
    '''
    fits.writeto(data=data, filename=name, overwrite=True)
    #img.save(name, compression=None, x_resolution=data.shape[0], y_resolution=data.shape[1])
    return os.path.abspath(name)


def drago(data: numpy.ndarray, bias: float = 0.85) -> numpy.ndarray:
    ''' Tonemap the image using gradient domain compression as described in Fattal et al. 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        bias : float
            Bias the logarithmic base towards lower or higher values. Default: 0.85.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_drago.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo drago '
    if bias is not None:
        cmd += f'--tmoDrgBias {bias} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    # os.remove(tmpname)
    # os.remove(tmpname_out)
    return data_tm


def duran(data: numpy.ndarray, sigma_spatial: float = None, sigma_range: float = None, base_contrast: float = None):
    ''' Tonemap the image using gradient domain compression as described in Durand and Dorsey et al. 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        sigma_spatial : float
            Spatial kernal size. Default: None.
        sigma_range : float
            Range kernel size. Default: None.
        base_contrast : float
            Base contrast. Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_durand.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo durand '
    if sigma_spatial is not None:
        cmd += f'--tmoDurSigmaS {sigma_spatial} '
    if sigma_range is not None:
        cmd += f'--tmoDurSigmaR {sigma_range} '
    if base_contrast is not None:
        cmd += f'--tmoDurBase {base_contrast} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    # os.remove(tmpname)
    # os.remove(tmpname_out)
    return data_tm


def fattal(data: numpy.ndarray, alpha: float = None, beta: float = None, colour_saturation: float = None, noise: float = None) -> numpy.ndarray:
    ''' Tonemap the image using gradient domain compression as described in Fattal et al. 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        alpha : float
            Controls which gradient magnitude is preserved when beta is 1. Default is None.
        beta : float
            Controls local detail enhancement.
        colour_saturation : float
            Controls colour saturation.
        noise : float
            Controls the threshold for what is seen as noise, where detail enhancement is reduced.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_fattal.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo fattal '
    if alpha is not None:
        cmd += f'--tmoFatAlpha {alpha} '
    if beta is not None:
        cmd += f'--tmoFatBeta {beta} '
    if colour_saturation is not None:
        cmd += f'--tmoFatAlpha {colour_saturation} '
    if noise is not None:
        cmd += f'--tmoFatAlpha {noise} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def ferradans(data: numpy.ndarray, rho: float = -2, inv_alpha: float = 5) -> numpy.ndarray:
    ''' Tonemap the image using the method described in Ferradans et al. 2011.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        rho : float
            Controls overall lightness. Larger values yield a brighter image. Default is -2.
        beta : float
            Controls detail enhancement. Larger values yield more detail enhancement. Default is 5.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_ferradans.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo ferradans '
    if rho is not None:
        cmd += f'--tmoFerRho {rho} '
    if inv_alpha is not None:
        cmd += f'--tmoFerInvAlpha {inv_alpha} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def ferwerda(data, multiplier: float = None, luminance_adaptation: float = None) -> numpy.ndarray:
    ''' Tonemap the image using the method described in Ferwerda et al. 1996.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        multiplier : float
            Default is None.
        luminance_adaptation : float
            Default is None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_ferwerda.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo ferradans '
    if multiplier is not None:
        cmd += f'--tmoFerwerdaMul {multiplier} '
    if luminance_adaptation is not None:
        cmd += f'--tmoFerwerdaAdaptLum {luminance_adaptation} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def kimkautz(data, c1: float = None, c2: float = None) -> numpy.ndarray:
    ''' Tonemap the image using the human vision based method described in Kim and Kaus 2008.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        c1 : float
            Default is None.
        c2 : float
            Default is None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_kimkautz.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo kimkautz '
    if c1 is not None:
        cmd += f'--tmoKimKautzC1 {c1} '
    if c2 is not None:
        cmd += f'--tmoKimKautzC2 {c2} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def mantiuk06(data, contrast: float = None, saturation: float = None, detail: float = None, contrast_equalisation: bool = False) -> numpy.ndarray:
    ''' Tonemap the image using the human vision based method described in Mantiuk 2006.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        contrast : float
            Default is None.
        saturation : float
            Default is None.
        detail : float
            Default is None.
        contrast_equalisation : bool
            Default is False.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_mantiuk06.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo mantiuk06 '
    if contrast is not None:
        cmd += f'--tmoM06Contrast {contrast} '
    if saturation is not None:
        cmd += f'--tmoM06Saturation {saturation} '
    if detail is not None:
        cmd += f'--tmoM06Detail {detail} '
    cmd += f'--tmoM06ContrastEqual {contrast_equalisation} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def mantiuk08(data, contrast_enhancement: float = None, colour_saturation: float = None, luminance_level: float = None, set_luminance: bool = False) -> numpy.ndarray:
    ''' Tonemap the image using the human vision based method described in Mantiuk 2008.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        contrast_enhancement : float
            Default is None.
        colour_saturation : float
            Default is None.
        luminance_level : float
            Used when set_luminance is True. Default is None.
        set_luminance : bool
            Default is False

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_mantiuk08.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo mantiuk08 '
    if colour_saturation is not None:
        cmd += f'--tmoM08ColorSaturation {colour_saturation} '
    if contrast_enhancement is not None:
        cmd += f'--tmoM08ContrastEnh {contrast_enhancement} '
    if luminance_level is not None:
        cmd += f'--tmoM08LuminanceLvl {luminance_level} '
    cmd += f'--tmoM08SetLuminance {set_luminance} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def reinhard02(data, key: float = None, phi: float = None, use_scales: bool = True, range: float = None, low: float = None, high: float = None) -> numpy.ndarray:
    ''' Tonemap the image using the human vision based method described in Mantiuk 2008.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        key : float
            Default: None.
        phi : float
            Default: None.
        use_scales : bool
            Default: True.
        range : float
            Default: None.
        low : float
            Default: None.
        high : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_reinhard02.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo reinhard02 '
    if key is not None:
        cmd += f'--tmoR02Key {key} '
    if phi is not None:
        cmd += f'--tmoR02Phi {phi} '
    if use_scales is not None:
        cmd += f'--tmoR02Scales {use_scales} '
    if range is not None:
        cmd += f'--tmoR02Num {range} '
    if low is not None:
        cmd += f'--tmoR02Low {low} '
    if high is not None:
        cmd += f'--tmoR02High {high} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def reinhard05(data, brightness: float = None, chroma: float = None, lightness: float = None) -> numpy.ndarray:
    ''' Tonemap the image using the human vision based method described in Reinhard 2005.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        brightness : float
            Default: None.
        chroma : float
            Default: None.
        lightness : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    '''
    tmpname = _store_tmpfile(data, 'tmp_reinhard05.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = BASECOMMAND + ' -e 0 --tmo reinhard05 '
    if brightness is not None:
        cmd += f'--tmoR05Brightness {brightness} '
    if chroma is not None:
        cmd += f'--tmoR05Chroma {chroma} '
    if lightness is not None:
        cmd += f'--tmoR05Lightness {lightness} '
    cmd += f'-o {tmpname_out} {tmpname}'
    run_command(cmd.split(' '))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm

