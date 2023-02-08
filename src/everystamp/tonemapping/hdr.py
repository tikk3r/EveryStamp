'''Sub-module containing wrappers to LuminosityHDR's HDR tonemapping programs.'''
import os
import subprocess

from skimage import io
from astropy.io import fits
import numpy

BASECOMMAND = 'luminance-hdr-cli'


def run_command(cmd):
    print('Running command:', ' '.join(cmd))
    try:
        subprocess.run(cmd)
    except subprocess.CalledProcessError as e:
        print("Tone mapping failed with error code", e.returncode)
        print(e.output)


def _load_tonemapped_tmpdata(name: str) -> numpy.ndarray:
    return io.imread(name, as_gray=True)


def _store_tmpfile(data: numpy.ndarray, name: str) -> str:
    fits.writeto(data=data, filename=name, overwrite=True)
    #img.save(name, compression=None, x_resolution=data.shape[0], y_resolution=data.shape[1])
    return os.path.abspath(name)


def fattal(data: numpy.ndarray, alpha: float = None, beta: float = None, colour: float = None, noise: float = None):
    tmpname = _store_tmpfile(data, 'tmp_fattal.fits')
    tmpname_out = tmpname.replace('.fits', '.tonemapped.tiff')
    cmd = [BASECOMMAND, '--tmo fattal']
    if alpha is not None:
        cmd.append(f'--tmoFatAlpha {alpha}')
    if beta is not None:
        cmd.append(f'--tmoFatBeta {beta}')
    if colour is not None:
        cmd.append(f'--tmoFatAlpha {alpha}')
    if noise is not None:
        cmd.append(f'--tmoFatAlpha {alpha}')
    cmd.append('-o ' + tmpname_out)
    cmd.append('-e 0')
    cmd.append(tmpname)
    run_command(cmd)
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    return data_tm