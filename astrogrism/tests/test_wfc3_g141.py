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
G141_IMAGE_FILE = 'https://github.com/npirzkal/aXe_WFC3_Cookbook/raw/main/cookbook_data/G141/ib6o23rsq_flt.fits' # noqa


def test_wfc3_g141_astropywcs():
    """
    Tests the Astrogrism Detector > World Transform against
    the built in Astropy WCS model
    """
    # Define Pixel Grid
    xx = np.arange(0, 1014, 20)
    yy = np.arange(0, 1014, 20)
    # Download file and initialize wcs
    fn = download_file(G141_IMAGE_FILE, cache=True)
    grism_image_hdulist = fits.open(fn)
    astropy_wcs = WCS(grism_image_hdulist['SCI'].header)
    # Convert Pixel Grid to Sky Grid
    astropy_coords = astropy_wcs.pixel_to_world(xx, yy)
    # Extract sky coordinates
    astropywcs_ra = astropy_coords.ra.value
    astropywcs_dec = astropy_coords.dec.value

    # Init GrismObs
    fn = download_file(G141_IMAGE_FILE, cache=True)
    grismobs = GrismObs(fn)
    # Retrieve Transform
    image2world = grismobs.geometric_transforms.get_transform('detector',
                                                              'world')
    # Calculate sky coordinates on grid
    astrogrism_ra, astrogrism_dec = image2world(xx, yy, 0, 0)[:2]

    # Compare results
    np.testing.assert_allclose(astrogrism_ra, astropywcs_ra, atol=5e-02)
    np.testing.assert_allclose(astrogrism_dec, astropywcs_dec, atol=5e-02)


def test_wfc3_g141_grizli():
    """
    Tests the Astrogrism Detector > Grism Transform against
    grizli's grism transform
    """
    # TODO: Calculate Grizli values dynamically rather than rely on stored vals

    # Calculate Astrogrism dispersed Xs and Ys for given wavelengths
    grizli_wavelengths_file = pathlib.Path(test_dir,
                                           'grizli_wfc3_g141_solutions',
                                           'grizli_wavelength_solutions.npy')
    grizli_wavelengths = np.load(grizli_wavelengths_file)
    astrogrism_x, astrogrism_y = _image2grism(500,
                                              500,
                                              grizli_wavelengths,
                                              G141_IMAGE_FILE)

    # Compare results
    np.testing.assert_allclose(astrogrism_x, 500 + np.arange(500), atol=5e-02)
    grizli_y_file = pathlib.Path(test_dir,
                                 'grizli_wfc3_g141_solutions',
                                 'grizli_y_solutions.npy')
    grizli_y = np.load(grizli_y_file)
    np.testing.assert_allclose(astrogrism_y, grizli_y, atol=5e-02)


def _image2grism(x_center, y_center, wavelengths, grism_file=None):
    if not grism_file:
        raise NotImplementedError("Grism FLT File is required for now")
    fn = download_file(grism_file, cache=True)
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
