from astrogrism import GrismObs
import grismconf

import pathlib
from astropy.utils.data import download_file
import numpy as np

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()
G141_IMAGE_FILE = 'https://github.com/npirzkal/aXe_WFC3_Cookbook/raw/main/cookbook_data/G141/ib6o23rsq_flt.fits' # noqa


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
