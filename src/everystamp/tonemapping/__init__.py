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