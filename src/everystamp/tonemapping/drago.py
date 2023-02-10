'''Sub-module to perform HDR tonemapping according to Drago  et al. 2003

https://resources.mpi-inf.mpg.de/tmo/logmap/logmap.pdf.
'''
import cv2
import numpy as np
import scipy


def bias(t, b):
    biast = t ** (np.log(b) / np.log(0.5))
    return biast


def map_drago(data, b):
    # World luminance.
    Lw = data
    # Maximum world luminance.
    Lwmax = np.nanmax(data)
    # Maximum luminance in the output image.
    Ldmax = 1000#0.1 * Lwmax

    Ld = (Ldmax * 0.01) / np.log10(Lwmax + 1) * np.log(Lw + 1) / np.log(2.0 + bias(Lw / Lwmax, b) * 8)
    # Ld = (np.log(Lw + 1)) / np.log(2.0 + bias(Lw / Lwmax, b) * 8)
    # Ld = (Lw / Lwmax) * (np.log(Lw + 1)) / np.log(2.0 + bias(Lw / Lwmax, b) * 8)
    return Ld