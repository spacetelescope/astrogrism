from astrogrism import GrismObs
import grismconf

import pathlib
import tempfile
from zipfile import ZipFile
from astropy.utils.data import download_file
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()
G102_IMAGE_FILE = str(test_dir / 'data'/ 'IRG102_icwz15e7q_flt.fits') # noqa


def test_wfc3_g102_astropywcs():
    """
    Tests the Astrogrism Detector > World Transform against
    the built in Astropy WCS model
    """
    # Define Pixel Grid
    xx = np.arange(0, 1014, 20)
    yy = np.arange(0, 1014, 20)
    # Download file and initialize wcs
    fn = G102_IMAGE_FILE
    grism_image_hdulist = fits.open(fn)
    astropy_wcs = WCS(grism_image_hdulist['SCI'].header)
    # Convert Pixel Grid to Sky Grid
    astropy_coords = astropy_wcs.pixel_to_world(xx, yy)
    # Extract sky coordinates
    astropywcs_ra = astropy_coords.ra.value
    astropywcs_dec = astropy_coords.dec.value

    # Init GrismObs
    fn = G102_IMAGE_FILE
    grismobs = GrismObs(fn)
    # Retrieve Transform
    image2world = grismobs.geometric_transforms.get_transform('detector',
                                                              'world')
    # Calculate sky coordinates on grid
    astrogrism_ra, astrogrism_dec = image2world(xx, yy, 0, 0)[:2]

    # Compare results
    np.testing.assert_allclose(astrogrism_ra, astropywcs_ra, atol=5e-02)
    np.testing.assert_allclose(astrogrism_dec, astropywcs_dec, atol=5e-02)


def _image2grism(x_center, y_center, wavelengths, grism_file=None):
    if not grism_file:
        raise NotImplementedError("Grism FLT File is required for now")
    fn = G102_IMAGE_FILE
    grismobs = GrismObs(fn)
    image2grism = grismobs.geometric_transforms.get_transform('detector',
                                                              'grism_detector')
    x, y = list(), list()
    for val in wavelengths:
        dispersion = image2grism.evaluate(x_center,
                                          y_center,
                                          wavelength=val,
                                          order=1)
        x.append(dispersion[0])
        y.append(dispersion[1])
    return x, y
