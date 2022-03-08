from math import floor
from os import environ
from pathlib import Path
from shutil import copy
from tarfile import open as tar_open
from tempfile import gettempdir
from warnings import warn

from astropy.utils.data import download_file
from astropy.io.fits import ImageHDU
from synphot import Observation
from astropy.wcs import WCS
import numpy as np
import progressbar

from astrogrism import GrismObs


def download_stsynphot_files(download_path):

    # Prepare environment required for stsynphot to be imported
    environ['PYSYN_CDBS'] = str(download_path / 'grp' / 'redcat' / 'trds')

    # Download HST Instrument Data Archive
    download_path.mkdir(parents=True, exist_ok=True)
    hst_data_files_archive = Path(download_file('https://ssb.stsci.edu/'
                                                'trds/tarfiles/synphot1.tar.gz', cache=True))
    with tar_open(hst_data_files_archive) as tar:
        # Only extract the files that are missing
        for file in tar:
            if not (download_path / Path(file.name)).exists():
                tar.extract(file, path=download_path)
    # Download Vega CALSPEC Reference Atlas
    vega_reference_atlas_path = download_path / 'grp/redcat/trds/calspec/alpha_lyr_stis_010.fits'
    # Check if it exists first before trying to download it
    if not vega_reference_atlas_path.exists():
        vega_reference_atlas_path.parent.mkdir(parents=True, exist_ok=True)
        archive_url = ('https://archive.stsci.edu/hlsps/reference-atlases/cdbs/'
                       'current_calspec/alpha_lyr_stis_010.fits')
        temp_download = Path(download_file(archive_url, cache=True))
        copy(str(temp_download), str(vega_reference_atlas_path))


def generate_synthetic_spectrum(grism, detector=None, temp_path=gettempdir(), verbose=False):
    """
    Initializes and uses STSynphot to generate a Vega spectrum within the bandpass of a given grism

    Parameters
    ----------
    grism : str
        String representation of one of the four supported HST Grisms
        Valid grisms: G141, G102, G280, G800L

    detector : int
        For detectors with multiple chips, specifies which chip to simulate
        Only useful for G280 and G800L Grisms

    temp_path : str
        Path to download necessary files for STSynphot. Fallsback to Python's
        default temporary folder location

    """
    if detector not in (1, 2, None):
        raise ValueError("Invalid detector argument. Please choose 1 or 2")

    SIM_DATA_DIR = Path(temp_path) / "astrogrism_simulation_files"

    download_stsynphot_files(SIM_DATA_DIR)

    # Now that we have all our reference files, we can import stsynphot
    # (This is why it's not a top-line import)
    from stsynphot import Vega, band
    if grism == 'G141':
        if detector and verbose:
            warn("WFC3's G141 grism does not have multiple detectors. Ignoring detector argument",
                 RuntimeWarning)
        bandpass = band('wfc3,ir,g141')
    elif grism == 'G102':
        if detector and verbose:
            warn("WFC3's G102 grism does not have multiple detectors. Ignoring detector argument",
                 RuntimeWarning)
        bandpass = band('wfc3,ir,g102')
    elif grism == 'G280':
        bandpass = band(f'wfc3,uvis{detector},g280')
    elif grism == 'G800L':
        bandpass = band(f'acs,wfc{detector},g800l')
    else:
        raise ValueError(f"Unrecognized grism: {grism}. Valid grisms: G141, G102, G280, G800L")

    spectrum = Observation(Vega, bandpass, binset=bandpass.binset).to_spectrum1d()

    # Find the first value with a non-zero "flux"
    for i in range(len(spectrum.flux.value)):
        value = spectrum.flux.value[i]
        if value != 0.0:
            min_slice = i
            break

    # Find the last value with a non-zero "flux"
    for i in reversed(range(len(spectrum.flux.value))):
        value = spectrum.flux.value[i]
        if value != 0.0:
            max_slice = i
            break
    return spectrum[min_slice:max_slice]


def disperse_spectrum_on_image(grism_file, wide_field_image, spectrum, detector=None):

    grismobs = GrismObs(grism_file)
    if detector:
        image2grism = grismobs.geometric_transforms[f'CCD{detector}'].get_transform('detector', 'grism_detector')
    else:
    image2grism = grismobs.geometric_transforms.get_transform('detector', 'grism_detector')

    if type(wide_field_image) is ImageHDU:
        data = wide_field_image.data
    else:
        data = wide_field_image
    
    shape = data.shape
    simulated_data = np.zeros(shape)

    # Start progressbar
    bar = progressbar.ProgressBar(maxval=shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    # For each pixel in the science image, we need to disperse it's spectrum
    for x_pixel in range(0, shape[0]):
        bar.update(x_pixel)
        for y_pixel in range(0, shape[1]):
            # Get the flux of the science pixel; we'll need to scale the spectrum to this brightness
            data_flux = wide_field_image[1].data[x_pixel][y_pixel]
            if np.isnan(data_flux):
                continue

            # For each Wavelength in the spectrum, calculate where, in pixels, that wavelength would fall on the detector
            dispersed_coords = image2grism(x_pixel, y_pixel, spectrum.wavelength, 1)
            for step in range(0, len(spectrum.wavelength)):
                spectrum_flux = spectrum.flux[step]

                dispersed_x = dispersed_coords[0][step]
                dispersed_y = dispersed_coords[1][step]

                # If the dispersed position of the wavelength is inside the bounds of the image, write the spectrum
                if (0 < dispersed_x < simulated_data.shape[0]) and (0 < dispersed_y < simulated_data.shape[1]):
                    # Scale the flux of the spectrum to the brightness of the original pixel
                    # NOTE: Is floor the right approach to determine which pixel to write to? Maybe sufficient until we accomplish the "drizzling" part of the simulation?
                    simulated_data[floor(dispersed_x)][floor(dispersed_y)] = simulated_data[floor(dispersed_x)][floor(dispersed_y)] + (data_flux * spectrum_flux).value
    bar.finish()
    return simulated_data


def simulate_grism(grism, wide_field_image):
    spectrum = generate_simulation_spectrum(grism)
    disperse_spectrum_on_image(grism, wide_field_image, spectrum)
