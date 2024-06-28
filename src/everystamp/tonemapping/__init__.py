'''Sub-module for tone mapping images.'''
import numpy as np


def make_nonnegative(data):
    ''' Makes an array non-negative by adding the minimum entry to all entries.

    Args:
        data : numpy.ndarray
            Array to make non-negative.

    Returns:
        data_scaled : numpy.ndarray
            Non-negative version of data.
    '''
    minimum = np.nanmin(data)
    if minimum < 0:
        data_scaled = data - np.nanmin(data)
        return data_scaled
    else:
        return data


def normalise(x, vmin=None, vmax=None, clip=True):
    """ Normalise an array of values according to y = (x - vmin / (vmax - vmin).
    
    Args:
        x : numpy.ndarray
            Array of values to normalise.
        vmin : float
            Lower bound for normalisation.
        vmax : float
            Upper bound for normalisation.
        clip : bool
            Clip values to the range [0:1]?
    """
    if not vmin:
        vmin = np.nanmin(x)
    if not vmax:
        vmax = np.nanmax(x)
    y = (x - vmin) / (vmax - vmin)
    if clip:
        y = np.clip(y, 0, 1)
        return y
    else:
        return y


def gamma(data, gamma):
    ''' Apply a gamma correction to the image, raising it to the power 1/gamma.

    Args:
        data : numpy.ndarray
            Data which to gamma correct.

    Returns:
        gamma : float
            The gamma correction factor.
    '''
    return data ** (1 / gamma)
