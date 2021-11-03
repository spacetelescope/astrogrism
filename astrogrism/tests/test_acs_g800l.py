from astrogrism import GrismObs

import pathlib
import numpy as np
from astropy.tests.helper import assert_quantity_allclose
import astropy.units as u

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()


def test_acs_g800l_roundtrip():
    """
    Tests Astrogrism transforms round tripping as expected between grism
    detector, direct detector and world frames.
    """
    test_file = pathlib.Path(test_dir, "data", "acs_test_file.fits")
    grism_obs = GrismObs(str(test_file))

    assert grism_obs.geometric_transforms["CCD1"].available_frames == ['grism_detector',
                                                                       'detector', 'world']

    # Check that detector -> grism is working before checking full transform
    d2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("detector", "grism_detector")
    d2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("detector", "grism_detector")

    d2g1_expected = (1084.130391789169, 2044.4123782198496, 1024.0, 2048.0, 1.0)
    d2g2_expected = (1080.3572828816784, 2045.548766180373, 1024.0, 2048.0, 1.0)

    np.testing.assert_allclose(d2g1(1024.0, 2048.0, 0.7, 1.0), d2g1_expected, atol=5e-2)
    np.testing.assert_allclose(d2g2(1024.0, 2048.0, 0.7, 1.0), d2g2_expected, atol=5e-2)

    # Now test transforming all the way end to end
    world_ref = (264.1018298204877, -32.90801703429679, 0.7*u.Unit("micron"), 1.0)

    g2w1 = grism_obs.geometric_transforms["CCD1"].get_transform("grism_detector", "world")
    w2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("world", "grism_detector")

    g2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("grism_detector", "world")
    w2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "grism_detector")

    w2g1_expected = (2103.325418, 1020.241261, 2047.007260, 1022.979833, 1.0)
    w2g2_expected = (2099.471492, 1020.803474, 2047.001747, 1023.005213, 1.0)

    np.testing.assert_allclose(w2g1(*world_ref), w2g1_expected, atol=5e-5)
    np.testing.assert_allclose(w2g2(*world_ref), w2g2_expected, atol=5e-5)

    # Use rtol here because the wavelength doesn't round trip perfectly
    g2w1_res = g2w1(*w2g1_expected)
    g2w2_res = g2w2(*w2g2_expected)

    [assert_quantity_allclose(g2w1_res[i], world_ref[i], rtol=0.005) for i in
     range(len(world_ref))]
    [assert_quantity_allclose(g2w2_res[i], world_ref[i], rtol=0.005) for i in
     range(len(world_ref))]
