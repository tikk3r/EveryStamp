'''Sub-module for tone mapping images.'''
import numpy as np


def make_nonnegative(data):
    minimum = np.nanmin(data)
    if minimum < 0:
        data_scaled = data - np.nanmin(data)
        return data_scaled
    else:
        return data


def gamma(data, gamma):
    return data ** gamma