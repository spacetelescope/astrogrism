## Base class for Astrogrism
from importlib.resources import path as resource_path

import asdf
from astropy.io import fits
from astropy.modeling import models
from astropy import units as u
from gwcs import wcs as gwcs
from gwcs import coordinate_frames as cf
import numpy as np
import pathlib

from .HST.transform_models import WFC3IRForwardGrismDispersion, WFC3IRBackwardGrismDispersion
from .HST.dispersion_models import DISPXY_Extension

# Will remove this hardcoded path once I have this packaged up for use with importlib
#pkg_dir = pathlib.Path('/Users/rosteen/projects/astrogrism_sandbox')
pkg_dir = pathlib.Path(__file__).parent.absolute()

class GrismObs():
    """
    Base class for astrogrism package. Stores all necessary information about
    a single grism observation.
    """

    def __init__(self, grism_image, direct_image=None, telescope=None, instrument=None,
                 detector=None, filter=None):

        for image in ('grism_image', 'direct_image'):
            image_obj = locals().get(image)
            if isinstance(image_obj, str):
                setattr(self, image, fits.open(pathlib.Path(image_obj)))
            elif isinstance(image_obj, (fits.HDUList, type(None))):
                setattr(self, image, image_obj)
            else:
                raise TypeError(f"Unrecognized type: {image} must be filepath or FITS HDUList ")

        # For now, the grism_image is a required argument. Eventually we want to
        # break that, but for now, raise an error
        if not grism_image:
            raise TypeError("grism_image is a required argument")

        # Parse grism image file header for meta info
        self.grism_header = self.grism_image["PRIMARY"].header

        # Attempt to retrieve any information missing from the header (e.g. SIP)
        # Should probably make these properties instead.
        if telescope is None:
            self.telescope = self.grism_header.get("TELESCOP")

        if instrument is None:
            self.instrument = self.grism_header.get("INSTRUME")

        if filter is None:
            self.filter = self.grism_header.get("FILTER")


        # Build GWCS geometric transform pipeline
        self._build_geometric_transforms()

    def _build_geometric_transforms(self):

        """
        Build transform pipeline under the hood so the user doesn't need
        to worry about it.

        TODO:
        - Try to get SIP coefficients from grism observation header before
        resorting to premade file
        """

        # Register custom asdf extension
        asdf.get_config().add_extension(DISPXY_Extension())

        # Get paths to premade configuration files
        #config_dir = "{}/config/{}/".format(pkg_dir, self.telescope)
        config_dir = pkg_dir / self.telescope / 'config'

        if self.telescope == "HST":
            if self.filter in ("G102", "G141"):
                instrument = self.instrument + "_IR"
            elif self.filter == "G280":
                instrument = self.instrument + "_UV"
        else:
            instrument = self.instrument
        sip_file = config_dir / "{}_distortion.fits".format(instrument)
        spec_wcs_file = config_dir / "{}_{}_specwcs.asdf".format(
                                                       self.instrument,
                                                       self.filter)

        # Build the grism_detector <-> detector transforms
        specwcs = asdf.open(str(spec_wcs_file)).tree
        displ = specwcs['displ']
        dispx = specwcs['dispx']
        dispy = specwcs['dispy']
        invdispl = specwcs['invdispl']
        invdispx = specwcs['invdispx']
        invdispy = specwcs['invdispy']
        orders = specwcs['order']

        gdetector = cf.Frame2D(name='grism_detector',
                       axes_order=(0, 1),
                       unit=(u.pix, u.pix))
        det2det = WFC3IRForwardGrismDispersion(orders,
                                               lmodels=displ,
                                               xmodels=invdispx,
                                               ymodels=dispy)
        det2det.inverse = WFC3IRBackwardGrismDispersion(orders,
                                                        lmodels=invdispl,
                                                        xmodels=dispx,
                                                        ymodels=dispy)

        grism_pipeline = [(gdetector, det2det)]

        # Now add the detector -> world transform
        sip_hdus = fits.open(str(sip_file))

        acoef = dict(sip_hdus[1].header['A_*'])
        a_order = acoef.pop('A_ORDER')
        bcoef = dict(sip_hdus[1].header['B_*'])
        b_order = bcoef.pop('B_ORDER')

        # Get the inverse SIP polynomial coefficients from file
        apcoef = dict(sip_hdus[1].header['AP_*'])
        bpcoef = dict(sip_hdus[1].header['BP_*'])

        try:
            ap_order = apcoef.pop('AP_ORDER')
            bp_order = bpcoef.pop('BP_ORDER')
        except ValueError:
            raise

        crpix = [sip_hdus[1].header['CRPIX1'], sip_hdus[1].header['CRPIX2']]

        crval = [self.grism_image[1].header['CRVAL1'],
                 self.grism_image[1].header['CRVAL2']]

        cdmat = np.array([[sip_hdus[1].header['CD1_1'], sip_hdus[1].header['CD1_2']],
                  [sip_hdus[1].header['CD2_1'], sip_hdus[1].header['CD2_2']]])

        # Gather Forward SIP Model Polynomials
        a_polycoef = {}
        for key in acoef:
            a_polycoef['c' + key.split('A_')[1]] = acoef[key]

        b_polycoef = {}
        for key in bcoef:
            b_polycoef['c' + key.split('B_')[1]] = bcoef[key]

        # Gather Inverse SIP Model Polynomials
        ap_polycoef = {}
        for key in apcoef:
            ap_polycoef['c' + key.split('AP_')[1]] = apcoef[key]

        bp_polycoef = {}
        for key in bpcoef:
             bp_polycoef['c' + key.split('BP_')[1]] = bpcoef[key]

        # Construct models
        a_poly = models.Polynomial2D(a_order, **a_polycoef)
        b_poly = models.Polynomial2D(b_order, **b_polycoef)
        ap_poly = models.Polynomial2D(ap_order, **ap_polycoef)
        bp_poly = models.Polynomial2D(bp_order, **bp_polycoef)

        # See SIP definition paper for definition of u, v, f, g
        SIP_forward = (models.Shift(-(crpix[0]-1)) & models.Shift(-(crpix[1]-1)) | # Calculate u and v
             models.Mapping((0, 1, 0, 1, 0, 1)) | a_poly & b_poly & models.Identity(2) |
             models.Mapping((0, 2, 1, 3)) | models.math.AddUfunc() & models.math.AddUfunc() |
             models.AffineTransformation2D(matrix=cdmat) | models.Pix2Sky_TAN() |
             models.RotateNative2Celestial(crval[0], crval[1], 180))

        SIP_backward = (models.RotateCelestial2Native(crval[0], crval[1], 180) |
            models.Sky2Pix_TAN() | models.AffineTransformation2D(matrix=cdmat).inverse |
            models.Mapping((0, 1, 0, 1, 0, 1)) | ap_poly & bp_poly & models.Identity(2) |
            models.Mapping((0, 2, 1, 3)) | models.math.AddUfunc() & models.math.AddUfunc() |
            models.Shift((crpix[0]-1)) & models.Shift((crpix[1]-1)))

        full_distortion_model = SIP_forward & models.Identity(2)
        full_distortion_model.inverse = SIP_backward & models.Identity(2)

        imagepipe = []

        det_frame = cf.Frame2D(name="detector")
        imagepipe.append((det_frame, full_distortion_model))

        world_frame = cf.CelestialFrame(name="world", unit = (u.Unit("deg"), u.Unit("deg")),
                             axes_names=('lon', 'lat'), axes_order=(0, 1),
                             reference_frame="ICRS")
        imagepipe.append((world_frame, None))

        grism_pipeline.extend(imagepipe)

        self.geometric_transforms = gwcs.WCS(grism_pipeline)
