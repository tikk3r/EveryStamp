"""Sub-module containing wrappers to LuminosityHDR's HDR tonemapping programs."""
import os
import subprocess
from typing import Optional, Union

import cv2
import numpy
from astropy.io import fits
from PIL import Image
from skimage import io

BASECOMMAND = "luminance-hdr-cli"
if "EVERYSTAMP_LUMINANCE_HDR" in os.environ:
    BASECOMMAND = os.environ["EVERYSTAMP_LUMINANCE_HDR"]


def has_luminance_hdr():
    x = subprocess.run(
        BASECOMMAND,
        shell=True,
        env=dict(os.environ),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return not x.returncode


def run_command(cmd: Union[list, str]) -> None:
    """Run a command through subprocess.

    Args:
        cmd : list
            List of arguments building up the command or a string of the full command.

    Returns:
        None
    """
    print("Running command:", " ".join(cmd))
    try:
        if isinstance(cmd, list):
            subprocess.run(cmd)
        elif isinstance(cmd, str):
            subprocess.run(cmd.split(" "))
        else:
            raise ValueError("Malformed command: ", cmd)
    except subprocess.CalledProcessError as e:
        print("Tone mapping failed with error code", e.returncode)
        print(e.output)


def _load_tonemapped_tmpdata(name: str, as_gray: bool = True) -> numpy.ndarray:
    """Load the tonemapped output of LuminanceHDR ad return it as an array.

    Args:
        name : str
            File path to the data to load.
        as_gray : bool
            Treat the image as gray scale or not. Default: True.

    Returns:
        data : numpy.ndarray
            NumPy array containing the image.
    """
    # return imageio.imread(name)
    return io.imread(name, as_gray=as_gray)
    # return io.imread(name, as_gray=as_gray)


def _store_tmpfile(data: numpy.ndarray, name: str, header=None) -> str:
    """Store the data to tonemap to a temporary FITS file to pass to LuminanceHDR.

    Args:
        data : numpy.ndarray
            NumPy array containing the image.

    Returns:
        path : str
            Absolute path to the temporary file.
    """
    if name.lower().endswith("fits"):
        print(data.shape)
        fits.writeto(data=data.squeeze(), filename=name, header=header, overwrite=True)
    elif name.lower().endswith("exr"):
        cv2.imwrite(name, data.astype("float32"))
    else:
        im = Image.fromarray(data)
        im.convert(mode="RGB").save(name)
    return os.path.abspath(name)


def ashikhmin(
    data,
    eq2: bool = True,
    simple: Optional[float] = None,
    local_threshold: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Ashikmin 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        eq2 : bool
            Default: True.
        simple : bool
            Default: True.
        local_threshold : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    """
    tmpname = _store_tmpfile(data, "tmp_ashikmin.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo ashikmin "
    cmd += f"--tmoAshEq2 {eq2} "
    cmd += f"--tmoAshSimple {simple} "
    if local_threshold is not None:
        cmd += f"--tmoAshLocal {local_threshold} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def drago(data: numpy.ndarray, bias: float = 0.85) -> numpy.ndarray:
    """Tonemap the image using gradient domain compression as described in Fattal et al. 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        bias : float
            Bias the logarithmic base towards lower or higher values. Default: 0.85.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    """
    tmpname = _store_tmpfile(data, "tmp_drago.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo drago "
    if bias is not None:
        cmd += f"--tmoDrgBias {bias} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    # os.remove(tmpname)
    # os.remove(tmpname_out)
    return data_tm


def duran(
    data: numpy.ndarray,
    sigma_spatial: Optional[float] = None,
    sigma_range: Optional[float] = None,
    base_contrast: Optional[float] = None,
):
    """Tonemap the image using gradient domain compression as described in Durand and Dorsey et al. 2002.

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
    """
    tmpname = _store_tmpfile(data, "tmp_durand.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo durand "
    if sigma_spatial is not None:
        cmd += f"--tmoDurSigmaS {sigma_spatial} "
    if sigma_range is not None:
        cmd += f"--tmoDurSigmaR {sigma_range} "
    if base_contrast is not None:
        cmd += f"--tmoDurBase {base_contrast} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    # os.remove(tmpname)
    # os.remove(tmpname_out)
    return data_tm


def fattal(
    data: numpy.ndarray,
    alpha: Optional[float] = None,
    beta: Optional[float] = None,
    colour_saturation: Optional[float] = None,
    noise: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using gradient domain compression as described in Fattal et al. 2002.

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
    """
    tmpname = _store_tmpfile(data, "tmp_fattal.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo fattal "
    if alpha is not None:
        cmd += f"--tmoFatAlpha {alpha} "
    if beta is not None:
        cmd += f"--tmoFatBeta {beta} "
    if colour_saturation is not None:
        cmd += f"--tmoFatAlpha {colour_saturation} "
    if noise is not None:
        cmd += f"--tmoFatAlpha {noise} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def ferradans(
    data: numpy.ndarray, rho: float = -2, inv_alpha: float = 5
) -> numpy.ndarray:
    """Tonemap the image using the method described in Ferradans et al. 2011.

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
    """
    tmpname = _store_tmpfile(data, "tmp_ferradans.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo ferradans "
    if rho is not None:
        cmd += f"--tmoFerRho {rho} "
    if inv_alpha is not None:
        cmd += f"--tmoFerInvAlpha {inv_alpha} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def ferwerda(
    data,
    multiplier: Optional[float] = None,
    luminance_adaptation: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using the method described in Ferwerda et al. 1996.

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
    """
    tmpname = _store_tmpfile(data, "tmp_ferwerda.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo ferradans "
    if multiplier is not None:
        cmd += f"--tmoFerwerdaMul {multiplier} "
    if luminance_adaptation is not None:
        cmd += f"--tmoFerwerdaAdaptLum {luminance_adaptation} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def kimkautz(
    data, c1: Optional[float] = None, c2: Optional[float] = None
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Kim and Kaus 2008.

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
    """
    tmpname = _store_tmpfile(data, "tmp_kimkautz.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo kimkautz "
    if c1 is not None:
        cmd += f"--tmoKimKautzC1 {c1} "
    if c2 is not None:
        cmd += f"--tmoKimKautzC2 {c2} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def lischinski(data, alpha: Optional[float] = None):
    """Tonemap the image using the method described in Lischinski et al. 2006.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        alpha : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    """
    tmpname = _store_tmpfile(data, "tmp_lischinski.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo lischinski "
    if alpha is not None:
        cmd += f"--tmoLischinskiAlpha {alpha} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def mantiuk06(
    data,
    contrast: Optional[float] = None,
    saturation: Optional[float] = None,
    detail: Optional[float] = None,
    contrast_equalisation: bool = False,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Mantiuk 2006.

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
    """
    tmpname = _store_tmpfile(data, "tmp_mantiuk06.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo mantiuk06 "
    if contrast is not None:
        cmd += f"--tmoM06Contrast {contrast} "
    if saturation is not None:
        cmd += f"--tmoM06Saturation {saturation} "
    if detail is not None:
        cmd += f"--tmoM06Detail {detail} "
    cmd += f"--tmoM06ContrastEqual {contrast_equalisation} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def mantiuk08(
    data,
    contrast_enhancement: Optional[float] = None,
    colour_saturation: Optional[float] = None,
    luminance_level: Optional[float] = None,
    set_luminance: bool = False,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Mantiuk 2008.

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
    """
    tmpname = _store_tmpfile(data, "tmp_mantiuk08.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo mantiuk08 "
    if colour_saturation is not None:
        cmd += f"--tmoM08ColorSaturation {colour_saturation} "
    if contrast_enhancement is not None:
        cmd += f"--tmoM08ContrastEnh {contrast_enhancement} "
    if luminance_level is not None:
        cmd += f"--tmoM08LuminanceLvl {luminance_level} "
    cmd += f"--tmoM08SetLuminance {set_luminance} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    # os.remove(tmpname)
    # os.remove(tmpname_out)
    return data_tm


def reinhard02(
    data,
    key: Optional[float] = None,
    phi: Optional[float] = None,
    use_scales: bool = True,
    range: Optional[float] = None,
    low: Optional[float] = None,
    high: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Mantiuk 2008.

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
    """
    tmpname = _store_tmpfile(data, "tmp_reinhard02.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo reinhard02 "
    if key is not None:
        cmd += f"--tmoR02Key {key} "
    if phi is not None:
        cmd += f"--tmoR02Phi {phi} "
    if use_scales is not None:
        cmd += f"--tmoR02Scales {use_scales} "
    if range is not None:
        cmd += f"--tmoR02Num {range} "
    if low is not None:
        cmd += f"--tmoR02Low {low} "
    if high is not None:
        cmd += f"--tmoR02High {high} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def reinhard05(
    data,
    brightness: Optional[float] = None,
    chroma: Optional[float] = None,
    lightness: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Reinhard 2005.

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
    """
    tmpname = _store_tmpfile(data, "tmp_reinhard05.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo reinhard05 "
    if brightness is not None:
        cmd += f"--tmoR05Brightness {brightness} "
    if chroma is not None:
        cmd += f"--tmoR05Chroma {chroma} "
    if lightness is not None:
        cmd += f"--tmoR05Lightness {lightness} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def pattanaik(
    data,
    multiplier: Optional[float] = None,
    local_tonemap: bool = True,
    auto_lum: bool = True,
    cone_level: Optional[float] = None,
    rod_level: Optional[float] = None,
) -> numpy.ndarray:
    """Tonemap the image using the human vision based method described in Pattanaik et al. 2000 and Pattanaik et al. 2002.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        multiplier : float
            Default: None.
        local_tonemap : bool
            Default: True.
        auto_lum : bool
            Default: True.
        cone_level : float
            Default: None.
        rod_level : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    """
    tmpname = _store_tmpfile(data, "tmp_pattanaik.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo pattanaik "
    if multiplier is not None:
        cmd += f"--tmoPatMultiplier {multiplier} "
    if cone_level is not None:
        cmd += f"--tmoPatCone {cone_level} "
    if rod_level is not None:
        cmd += f"--tmoPatCone {rod_level} "
    cmd += f"--tmoPatLocal {local_tonemap} "
    cmd += f"--tmoPatAutoLum {auto_lum} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm


def vanhateren(data, pupil_area: Optional[float] = None):
    """Tonemap the image using the method described in van Hateren 2006.

    Parameters set to None will take their default values as set in LuminanceHDR.

    Args:
        data : numpy.ndarray
            Input data to tonemap.
        pupil_area : float
            Default: None.

    Returns:
        data_tm : numpy.ndarray
            Tonemapped data.
    """
    tmpname = _store_tmpfile(data, "tmp_vanhateren.fits")
    tmpname_out = tmpname.replace(".fits", ".tonemapped.tiff")
    cmd = BASECOMMAND + " -e 0 --tmo vanhateren "
    if pupil_area is not None:
        cmd += f"--tmoVanHaterenPupilArea {pupil_area} "
    cmd += f"-o {tmpname_out} {tmpname}"
    run_command(cmd.split(" "))
    data_tm = _load_tonemapped_tmpdata(tmpname_out)
    os.remove(tmpname)
    os.remove(tmpname_out)
    return data_tm
