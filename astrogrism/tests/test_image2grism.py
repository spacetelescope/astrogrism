import pathlib
import tempfile
from zipfile import ZipFile

from astrogrism import GrismObs
import astropy.units as u
from astropy.utils.data import download_file
from astropy.tests.helper import assert_quantity_allclose
import grismconf
import numpy as np
import pytest

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()
testdata = [
    ('G141', str(test_dir / 'data' / 'IRG141_ib6o23rsq_flt.fits')),
    ('G102', str(test_dir / 'data' / 'IRG102_icwz15e7q_flt.fits'))
]


@pytest.mark.parametrize("grism,grism_image", testdata)
def test_wfc3_grismconf(grism, grism_image):
    """
    Tests the Astrogrism Detector > Grism Transform against
    grismconf's grism transform
    """
    # Download grismconf repo
    fn = download_file(
        'https://github.com/npirzkal/GRISM_WFC3/archive/refs/heads/master.zip',
        cache='update')

    with pathlib.Path(tempfile.TemporaryDirectory().name) as grismconfdir:
        with ZipFile(fn, 'r') as grismconfarchive:
            grismconfarchive.extractall(grismconfdir)

        # Initialize Grismconf and initial parameters
        grismconf_lookup = {
            'G141': pathlib.Path('GRISM_WFC3-master/IR/G141.conf'),
            'G102': pathlib.Path('GRISM_WFC3-master/IR/G102.conf')
        }

        C = grismconf.Config(grismconfdir / grismconf_lookup[grism])

    # Define the center coordinates
    x_center = 500
    y_center = 500

    # Test translation on 20 steps, arbitrarily
    t = np.arange(0, 1, 0.05)

    # Calculate X and Y offsets
    xoffsets = C.DISPX('+1', x_center, y_center, t)
    yoffsets = C.DISPY('+1', x_center, y_center, t)
    # Grab corresponding wavelengths
    grismconf_wavelengths = C.DISPL('+1', x_center, y_center, t)*u.AA
    # Calculate dispersed Xs and Ys
    grismconf_x = xoffsets+x_center
    grismconf_y = yoffsets+y_center

    # Calculate Astrogrism dispersed Xs and Ys for given wavelengths
    astrogrism_x, astrogrism_y = _image2grism(x_center,
                                              y_center,
                                              grismconf_wavelengths,
                                              grism_image)

    # Compare results
    np.testing.assert_allclose(astrogrism_x, grismconf_x, atol=5e-02)
    np.testing.assert_allclose(astrogrism_y, grismconf_y, atol=5e-02)

    # Test Roundtripping
    roundtrip_wavelengths = _grism2image(x_center,
                                         y_center,
                                         astrogrism_x,
                                         astrogrism_y,
                                         grism_image)
    assert_quantity_allclose(roundtrip_wavelengths, grismconf_wavelengths)


def _image2grism(x_center, y_center, wavelengths, grism_file=None):
    if not grism_file:
        raise NotImplementedError("Grism FLT File is required for now")

    grismobs = GrismObs(grism_file)
    image2grism = grismobs.geometric_transforms.get_transform('detector',
                                                              'grism_detector')
    x, y = list(), list()
    for val in wavelengths:
        print(f"Val: {val}, x_cen: {x_center}, y_cen: {y_center}")
        dispersion = image2grism.evaluate(x_center, y_center, val, 1)
        x.append(dispersion[0])
        y.append(dispersion[1])

    return x, y


def _grism2image(x_center, y_center, x_dispersion, y_dispersion, grism_file=None):
    if not grism_file:
        raise NotImplementedError("Grism FLT File is required for now")

    grismobs = GrismObs(grism_file)
    grism2image = grismobs.geometric_transforms.get_transform('grism_detector',
                                                              'detector')
    if len(x_dispersion) != len(y_dispersion):
        raise AttributeError("Non-equal quantities of coordinates provided")

    wavelengths = list()
    for x, y in zip(x_dispersion, y_dispersion):
        _, _, wavelength, _ = grism2image.evaluate(x, y, x0=x_center, y0=y_center, order=1)
        wavelengths.append(wavelength)

    return wavelengths
