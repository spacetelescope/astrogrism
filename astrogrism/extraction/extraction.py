from astropy.io import fits
import asdf
import numpy as np
import pathlib

from astrogrism.HST.dispersion_models import DISPXY_Model, DISPXY_Extension
from astrogrism.HST.transform_models import (AstrogrismForwardGrismDispersion,
                                             AstrogrismBackwardGrismDispersion)

__all__ = ['extract_2d_spectrum',]

pkg_dir = pathlib.Path(__file__).parent.parent.absolute()

def extract_2d_spectrum(data, ll_x, ll_y, ur_x, ur_y, ll_l = 1.1,
                      ur_l = 1.7, order = 1, grism=None):
    """
    Function to do a simple box cutout around a 2D spectrum.

    The input lower and upper x and y bounds are pixel coordinates from the
    direct image. The upper and lower wavelength bounds are in microns.

    TODO: Change this to take grism as input instead of lower and upper wavelengths
    and specwcs reference file, with default wavelength range based on grism.
    Currently assumes the G141 reference file and wavelength ranged that I've been
    working with.
    """
    asdf.get_config().add_extension(DISPXY_Extension())
    config_dir = pkg_dir / 'config' / 'HST'
    if grism == "G800L":
        instrument = "ACS"
    else:
        instrument = "WFC3"
    specwcs_ref = config_dir / "{}_{}_specwcs.asdf".format(instrument, grism)

    specwcs = asdf.open(specwcs_ref).tree
    displ = specwcs['displ']
    dispx = specwcs['dispx']
    dispy = specwcs['dispy']
    invdispl = specwcs['invdispl']
    invdispx = specwcs['invdispx']
    invdispy = specwcs['invdispy']
    orders = specwcs['order']

    det2det = AstrogrismForwardGrismDispersion(orders,
                                               lmodels=displ,
                                               xmodels=invdispx,
                                               ymodels=dispy)
    det2det.inverse = AstrogrismBackwardGrismDispersion(orders,
                                                        lmodels=invdispl,
                                                        xmodels=dispx,
                                                        ymodels=dispy)

    ll = det2det.inverse.evaluate(ll_x, ll_y, ll_l, order)
    ur = det2det.inverse.evaluate(ur_x, ur_y, ur_l, order)

    print(ll, ur)

    return data[int(ll[1]):int(ur[1])+1, int(ll[0]):int(ur[0])+1]

