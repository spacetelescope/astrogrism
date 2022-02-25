from astrogrism.simulate import generate_simulation_spectrum

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

import pytest

# Taken from wavelengthrange file
# We should eventual modify the test to actually read the wavelengthrange reference file
TEST_GRISMS = [
    ('G141', 1.0402*u.um, 1.7845*u.um),
    ('G102', 0.7401*u.um, 0.9878*u.um)
]


@pytest.mark.parametrize("grism,grism_min,grism_max", TEST_GRISMS)
def test_synthetic_spectrum_grism_bounds(grism, grism_min, grism_max):
    spectrum = generate_simulation_spectrum(grism)
    assert_quantity_allclose(spectrum.wavelength.min(), grism_min)
    assert_quantity_allclose(spectrum.wavelength.max(), grism_max)


def test_single_chip_detector_warning():
    with pytest.warns(RuntimeWarning):
        generate_simulation_spectrum('G141')

    with pytest.warns(RuntimeWarning):
        generate_simulation_spectrum('G102')


def test_invalid_detector():
    generate_simulation_spectrum("G141", 3)
