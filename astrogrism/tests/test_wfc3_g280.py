from astrogrism import GrismObs

import pathlib
import numpy as np

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()
grism_image_file = 'https://github.com/npirzkal/aXe_WFC3_Cookbook/raw/main/cookbook_data/G141/ib6o23rsq_flt.fits' # noqa


def test_wfc3_g280_roundtrip():
    """
    Tests Astrogrism transforms round tripping as expected between grism
    detector, direct detector and world frames.
    """
    test_file = pathlib.Path(test_dir,"test_data", "icwz15e7q_flt.fits")
    grism_obs = GrismObs(str(test_file))

    assert grism_obs.geometric_transforms["CCD1"].available_frames == ['grism_detector',
                                                                       'detector', 'world']

    # Check that detector -> grism is working before checking full transform
    d2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("detector", "grism_detector")
    d2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("detector", "grism_detector")

    d2g1_expected = (1247.1169881075878, 1680.2088037732872, 1500.0, 1500.0, 1.0)
    d2g2_expected = (750.5270115330625, 1671.8505273099233, 1000.0, 1500.0, 1.0)

    np.testing.assert_allclose(d2g1(1500.0, 1500.0, 5000, 1.0), d2g1_expected, atol=5e-2)
    np.testing.assert_allclose(d2g2(1000.0, 1500.0, 5000, 1.0), d2g2_expected, atol=5e-2)

    # Now test transforming all the way end to end
    world_ref = (206.4318333333, 26.41859444444, 4000, 1.0)

    g2w1 = grism_obs.geometric_transforms["CCD1"].get_transform("grism_detector", "world")
    w2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("world", "grism_detector")

    g2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("grism_detector", "world")
    w2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "grism_detector")

    w2g1_expected = (1869.7387435684468, 1245.8356517642842, 2047.2667314355476,
                     1070.8144590508862, 1.0)
    w2g2_expected = (1873.0879758597869, 1237.8689879427377, 2047.2665683303317,
                     1070.814770017469, 1.0)

    np.testing.assert_allclose(w2g1(*world_ref), w2g1_expected, atol=5e-5)
    np.testing.assert_allclose(w2g2(*world_ref), w2g2_expected, atol=5e-5)

    # Use rtol here because the wavelength doesn't round trip perfectly
    np.testing.assert_allclose(g2w1(*w2g1_expected), world_ref, rtol=0.005)
    np.testing.assert_allclose(g2w2(*w2g2_expected), world_ref, rtol=0.005)
