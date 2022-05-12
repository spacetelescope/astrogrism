from astrogrism import GrismObs

import pathlib
import numpy as np
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

# TODO: Switch to importlib
test_dir = pathlib.Path(__file__).parent.absolute()


def test_uvis_subarray():
    """
    Tests Astrogrism transforms round tripping as expected between grism
    frame, direct frame and world frames.
    """
    test_file = pathlib.Path(test_dir, "data", "uvis_subarray_test_file.fits")
    grism_obs = GrismObs(str(test_file))
    grism_hdr = grism_obs.grism_image[1].header

    assert grism_obs.geometric_transforms["CCD2"].available_frames == ['grism_frame',
                                                                       'direct_frame', 'world']

    # Check that direct_frame -> grism is working before checking full transform
    d2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("direct_frame", "grism_frame")

    d2g2_expected = (418.3477, 1673.770, 600.0, 1500.0, 1.0)

    np.testing.assert_allclose(d2g2(600.0, 1500.0, 4000, 1.0), d2g2_expected, atol=5e-2)

    # Now test transforming all the way end to end
    world_ref2 = (grism_hdr["CRVAL1"], grism_hdr["CRVAL2"], 4000*u.AA, 1.0)

    g2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("grism_frame", "world")
    w2g2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "grism_frame")

    for key, value in grism_hdr.items():
        print(f"{key}: {value}")
    w2g2_expected = (1.874515e+03, -3.091660e+02, grism_hdr["CRPIX1"]-1, grism_hdr["CRPIX2"]-1, 1.0)

    np.testing.assert_allclose(w2g2(*world_ref2), w2g2_expected, rtol=5e-6)

    # Use rtol here because the wavelength doesn't round trip perfectly
    g2w2_res = g2w2(*w2g2_expected)
    [assert_quantity_allclose(g2w2_res[i], world_ref2[i], rtol=0.005) for i in
     range(len(world_ref2))]

    # Test world <-> direct_frame transforms
    w2d_expected = (grism_hdr["CRPIX1"]-1, grism_hdr["CRPIX2"]-1, 4000*u.AA, 1.0)

    d2w2 = grism_obs.geometric_transforms["CCD2"].get_transform("direct_frame", "world")
    w2d2 = grism_obs.geometric_transforms["CCD2"].get_transform("world", "direct_frame")

    [assert_quantity_allclose(w2d2(*world_ref2)[i], w2d_expected[i], rtol=5e-5) for
     i in range(len(world_ref2))]

    d2w2_res = d2w2(*w2d_expected)
    [assert_quantity_allclose(d2w2_res[i], world_ref2[i], rtol=0.005) for i in
     range(len(world_ref2))]


def test_ir_subarray():
    """
    Tests Astrogrism transforms round tripping as expected between grism
    frame, direct frame and world frames.
    """
    test_file = pathlib.Path(test_dir, "data", "wfc3ir_subarray_test_file.fits")
    grism_obs = GrismObs(str(test_file))
    grism_hdr = grism_obs.grism_image[1].header

    assert grism_obs.geometric_transforms.available_frames == ['grism_frame',
                                                               'direct_frame', 'world']

    # Check that direct_frame -> grism is working before checking full transform
    d2g = grism_obs.geometric_transforms.get_transform("direct_frame", "grism_frame")

    d2g_expected = (94.127118654, 137.60409411, 137.0, 137.0, 1.0)

    np.testing.assert_allclose(d2g(137, 137, 7000*u.AA, 1), d2g_expected, atol=5e-2)

    # Now test transforming all the way end to end
    world_ref = (grism_hdr["CRVAL1"], grism_hdr["CRVAL2"], 7000*u.AA, 1.0)

    g2w = grism_obs.geometric_transforms.get_transform("grism_frame", "world")
    w2g = grism_obs.geometric_transforms.get_transform("world", "grism_frame")

    w2g_expected = (94.1271186539101, 137.60409410896412, grism_hdr["CRPIX1"]-1,
                    grism_hdr["CRPIX2"]-1, 1.0)

    np.testing.assert_allclose(w2g(*world_ref), w2g_expected, rtol=5e-6)

    # Use rtol here because the wavelength doesn't round trip perfectly
    g2w_res = g2w(*w2g_expected)
    [assert_quantity_allclose(g2w_res[i], world_ref[i], rtol=0.005) for i in
     range(len(world_ref))]

    # Test world <-> direct_frame transforms
    w2d_expected = (grism_hdr["CRPIX1"]-1, grism_hdr["CRPIX2"]-1, 7000*u.AA, 1.0)

    d2w = grism_obs.geometric_transforms.get_transform("direct_frame", "world")
    w2d = grism_obs.geometric_transforms.get_transform("world", "direct_frame")

    [assert_quantity_allclose(w2d(*world_ref)[i], w2d_expected[i], rtol=5e-5) for
     i in range(len(world_ref))]

    d2w_res = d2w(*w2d_expected)
    [assert_quantity_allclose(d2w_res[i], world_ref[i], rtol=0.005) for i in
     range(len(world_ref))]
