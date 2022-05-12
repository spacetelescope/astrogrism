from astrogrism import GrismObs

import pathlib
import numpy as np
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()


def test_wfc3_g280_roundtrip():
    """
    Tests Astrogrism transforms round tripping as expected between grism
    frame, direct frame and world frames.
    """
    test_file = pathlib.Path(test_dir, "data", "uvis_test_file.fits")
    grism_obs = GrismObs(str(test_file))

    assert grism_obs.geometric_transforms["CCD1"].available_frames == ['grism_frame',
                                                                       'direct_frame', 'world']
    assert grism_obs.geometric_transforms["CCD2"].available_frames == ['grism_frame',
                                                                       'direct_frame', 'world']

    # Check that direct_frame -> grism is working before checking full transform
    d2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("direct_frame", "grism_frame")
    d2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("direct_frame", "grism_frame")

    d2g1_expected = (1247.1169881075878, 1680.2088037732872, 1500.0, 1500.0, 1.0)
    d2g2_expected = (750.5270115330625, 1671.8505273099233, 1000.0, 1500.0, 1.0)

    np.testing.assert_allclose(d2g1(1500.0, 1500.0, 5000, 1.0), d2g1_expected, atol=5e-2)
    np.testing.assert_allclose(d2g2(1000.0, 1500.0, 5000, 1.0), d2g2_expected, atol=5e-2)

    # Now test transforming all the way end to end
    world_ref1 = (206.45684221694, 26.41827449813, 4000*u.AA, 1.0)
    world_ref2 = (206.4312669986, 26.41859812035, 4000*u.AA, 1.0)

    g2w1 = grism_obs.geometric_transforms["CCD1"].get_transform("grism_frame", "world")
    w2g1 = grism_obs.geometric_transforms["CCD1"].get_transform("world", "grism_frame")

    g2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("grism_frame", "world")
    w2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "grism_frame")

    w2g1_expected = (1869.5620631381307, 1199.8189959488238, 2046.9999576026378,
                     1025.000003223617, 1.0)
    w2g2_expected = (1872.897871807071, 1191.880590457291, 2047.000039365022,
                     1024.9999995682701, 1.0)

    np.testing.assert_allclose(w2g1(*world_ref1), w2g1_expected, atol=5e-5)
    np.testing.assert_allclose(w2g2(*world_ref2), w2g2_expected, atol=5e-5)

    # Use rtol here because the wavelength doesn't round trip perfectly
    g2w1_res = g2w1(*w2g1_expected)
    g2w2_res = g2w2(*w2g2_expected)
    [assert_quantity_allclose(g2w1_res[i], world_ref1[i], rtol=0.005) for i in
     range(len(world_ref1))]
    [assert_quantity_allclose(g2w2_res[i], world_ref2[i], rtol=0.005) for i in
     range(len(world_ref2))]

    # Test world <-> direct_frame transforms
    w2d_expected = (2047, 1025, 4000*u.AA, 1.0)

    d2w1 = grism_obs.geometric_transforms["CCD1"].get_transform("direct_frame", "world")
    w2d1 = grism_obs.geometric_transforms["CCD1"].get_transform("world", "direct_frame")

    d2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("direct_frame", "world")
    w2d2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "direct_frame")

    [assert_quantity_allclose(w2d1(*world_ref1)[i], w2d_expected[i], rtol=5e-5) for
     i in range(len(world_ref1))]
    [assert_quantity_allclose(w2d2(*world_ref2)[i], w2d_expected[i], rtol=5e-5) for
     i in range(len(world_ref2))]

    d2w1_res = d2w1(*w2d_expected)
    d2w2_res = d2w2(*w2d_expected)
    [assert_quantity_allclose(d2w1_res[i], world_ref1[i], rtol=0.005) for i in
     range(len(world_ref1))]
    [assert_quantity_allclose(d2w2_res[i], world_ref2[i], rtol=0.005) for i in
     range(len(world_ref2))]
