from astrogrism.simulate import generate_synthetic_spectrum

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

import pytest

# Taken from wavelengthrange file (0th Order)
# We should eventually modify the test to actually read the wavelengthrange reference file
TEST_GRISMS = [
    ('G141', 1.0402*u.um, 1.6998*u.um),
    ('G102', 0.7401*u.um, 1.2297*u.um)
]


@pytest.mark.parametrize("grism,grism_min,grism_max", TEST_GRISMS)
def test_synthetic_spectrum_grism_bounds(grism, grism_min, grism_max):
    spectrum = generate_synthetic_spectrum(grism)
    assert_quantity_allclose(spectrum.wavelength.min().quantity, grism_min, rtol=1e-1)
    assert_quantity_allclose(spectrum.wavelength.max().quantity, grism_max, rtol=1e-1)


def test_single_chip_detector_warning():
    with pytest.warns(RuntimeWarning, match="grism does not have multiple detectors."):
        generate_synthetic_spectrum('G141', detector=1, verbose=True)

    with pytest.warns(RuntimeWarning, match="grism does not have multiple detectors"):
        generate_synthetic_spectrum('G102', detector=2, verbose=True)


def test_invalid_detector():
    with pytest.raises(ValueError, match="Invalid detector"):
        generate_synthetic_spectrum("G141", 3)


def test_invalid_grism():
    with pytest.raises(ValueError, match="Unrecognized grism"):
        generate_synthetic_spectrum("X999")
