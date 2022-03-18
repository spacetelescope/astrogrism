import csv
import datetime
import pkg_resources

from asdf.tags.core import Software, HistoryEntry
from astropy import units as u
try:
    from jwst.datamodels.wcs_ref_models import WavelengthrangeModel
except ModuleNotFoundError as e:
    raise ModuleNotFoundError("The jwst package, though not required for "
                              "Astrogrism, is required to generate "
                              "reference files for HST. Try `pip install "
                              "jwst` and try again!") from e

from .common import common_reference_file_keywords


def create_grism_wavelengthrange(grism, chip=None, author="STScI",
                                 wavelengthrange=None, extract_orders=None):
    """Create a wavelengthrange reference file for specified grism.

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
    history = f"{grism} wavelengthrange"

    acs_order_dict = {"A": 1, "B": 0, "C": 2, "D": 3, "E": -1, "F": -2, "G": -3}

    if chip is None and grism in ("G800L", "G280"):
        raise ValueError("Must specify chip 1 or 2")

    if grism == "G800L":
        chip_str = f"CHIP{chip}"
        outname = f"ACS_G800L_CCD{chip}_wavelengthrange.asdf"
    elif grism == "G280":
        chip_str = f"UVIS{chip}"
        outname = f"WFC3_G280_CCD{chip}_wavelengthrange.asdf"
    elif grism in ("G102", "G141"):
        chip_str = "IR"
        outname = f"WFC3_{grism}_wavelengthrange.asdf"
    else:
        raise ValueError(f"Grism {grism} is not currently supported.")

    ref_kw = common_reference_file_keywords(reftype="wavelengthrange",
                                            title="Grism reference file",
                                            description=f"{grism}{chip} Wavelength Ranges",
                                            exp_type="GRISM",
                                            author=author,
                                            pupil="ANY",
                                            model_type="WavelengthrangeModel",
                                            filename=outname,
                                            )

    if wavelengthrange is None:
        # This is a list of tuples that specify the
        # order, filter, wave min, wave max
        wavelengthrange = []
        waverange_tab = "config/common/WFSS_Wavelength_Range_V1.1_with_ACS.tab"
        waverange_fname = pkg_resources.resource_filename("astrogrism", waverange_tab)

        with open(waverange_fname, "r") as f:
            freader = csv.reader(f, delimiter=" ")
            for row in freader:
                if row[2] != grism or row[3] != chip_str:
                    continue
                if grism == "G800L":
                    order = acs_order_dict[row[5]]
                else:
                    order = int(row[5])

                wavelengthrange.append((order, grism, float(row[6]), float(row[7])))

    # array of integers of unique orders
    orders = sorted(set((x[0] for x in wavelengthrange)))
    filters = sorted(set((x[1] for x in wavelengthrange)))

    # Default list of orders to extract for each grism
    if extract_orders is None:
        extract_orders = [('G102', [1]),
                          ('G141', [1]),
                          ('G280', [1]),
                          ('G800L', [1])
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
                         'homepage': 'https://github.com/spacetelescope/astrogrism',
                         'version': '0.0.0'})
    history['software'] = software
    ref.history = [history]
    ref.validate()
    ref.to_asdf(outname)
