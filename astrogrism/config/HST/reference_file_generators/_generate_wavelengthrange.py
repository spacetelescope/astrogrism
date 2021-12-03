import datetime

from asdf.tags.core import Software, HistoryEntry
from astropy import units as u
try:
    from jwst.datamodels.wcs_ref_models import WavelengthrangeModel
except ModuleNotFoundError as e:
    raise ModuleNotFoundError("The jwst package, though not required for "
                              "Astrogrism, is required to generate "
                              "reference files for HST. Try `pip install "
                              "jwst` and try again!") from e

from ._common import common_reference_file_keywords


def create_tsgrism_wavelengthrange(outname="wfc3_tsgrism_wavelengthrange.asdf",
                                   history="WFC3 TSGrism wavelengthrange",
                                   author="STScI",
                                   wavelengthrange=None,
                                   extract_orders=None):
    """Create a wavelengthrange reference file for WFC3 TSGRISM mode.

    Parameters
    ----------
    outname: str
        The output name of the file
    history: str
        History information about it's creation
    author: str
        Person or entity making the file
    wavelengthrange: list(tuples)
        A list of tuples that set the order, filter, and
        wavelength range min and max
    extract_orders: list[list]
        A list of lists that specify

    """
    ref_kw = common_reference_file_keywords(reftype="wavelengthrange",
                                            title="WFC3 TSGRISM reference file",
                                            description="WFC3 Grism-Filter Wavelength Ranges",
                                            exp_type="WFC3_TSGRISM",
                                            author=author,
                                            pupil="ANY",
                                            model_type="WavelengthrangeModel",
                                            filename=outname,
                                            )

    if wavelengthrange is None:
        # This is a list of tuples that specify the
        # order, filter, wave min, wave max
        wavelengthrange = [(-1, 'G102', 0.7488958984375, 1.1496958984375),
                           (0, 'G102', 0.7401, 1.2297),
                           (1, 'G102', 0.7496, 1.1979),
                           (2, 'G102', 0.7401, 1.1897),
                           (3, 'G102', 0.7571, 0.9878),
                           (-1, 'G141', 1.031, 1.7845),
                           (0, 'G141', 1.0402, 1.6998),
                           (1, 'G141', 0.9953, 1.7697),
                           (2, 'G141', 0.9702, 1.5903),
                           (3, 'G141', 1.0068, 1.3875),
                           (4, 'G141', 1.031, 1.7845),
                           ]

    # array of integers of unique orders
    orders = sorted(set((x[0] for x in wavelengthrange)))
    filters = sorted(set((x[1] for x in wavelengthrange)))

    # Nircam has not specified any limitation on the orders
    # that should be extracted by default yet so all are
    # included.
    if extract_orders is None:
        extract_orders = [('G102', [1]),
                          ('G141', [1]),
                          ]

    ref = WavelengthrangeModel()
    ref.meta.update(ref_kw)
    ref.meta.exposure.p_exptype = "WFC3_TSGRISM"
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.wavelengthrange = wavelengthrange
    ref.extract_orders = extract_orders
    ref.order = orders
    ref.waverange_selector = filters

    history = HistoryEntry({'description': history,
                            'time': datetime.datetime.utcnow()})
    software = Software({'name': 'wfc3_reftools.py',
                         'author': author,
                         'homepage': 'https://github.com/spacetelescope/astrogrism_sandbox',
                         'version': '0.0.0'})
    history['software'] = software
    ref.history = [history]
    ref.validate()
    ref.to_asdf(outname)
