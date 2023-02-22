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

