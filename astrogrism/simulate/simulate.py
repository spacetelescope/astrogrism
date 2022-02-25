from os import environ
from pathlib import Path
from tarfile import open as tar_open
from tempfile import mkdtemp
from warnings import warn

from astropy.utils.data import download_file

SIM_DATA_DIR = Path(mkdtemp()) / "astrogrism_simulation_files"


def generate_simulation_spectrum(grism, detector=None):
    if detector not in (1, 2, None):
        raise ValueError("Invalid detector argument. Please choose 1 or 2")

    # Prepare environment required for stsynphot to be imported
    environ['PYSYN_CDBS'] = str(SIM_DATA_DIR / 'grp' / 'redcat' / 'trds')

    # Always download our data files fresh. Ensures changes from redcat team will be picked up
    # Download HST Instrument Data Files
    SIM_DATA_DIR.mkdir(parents=True, exist_ok=True)
    hst_data_files_archive = Path(download_file('https://ssb.stsci.edu/trds/tarfiles/synphot1.tar.gz')) # noqa
    with tar_open(hst_data_files_archive) as tar:
        tar.extractall(path=SIM_DATA_DIR)
    # Download Vega CALSPEC Reference Atlas
    calspec_dir = SIM_DATA_DIR / 'grp' / 'redcat' / 'trds' / 'calspec'
    calspec_dir.mkdir(parents=True, exist_ok=True)
    vega_reference_atlas = Path(download_file('https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/alpha_lyr_stis_010.fits')) # noqa
    vega_reference_atlas.replace(calspec_dir / 'alpha_lyr_stis_010.fits')

    # Now that we have all our reference files, we can import stsynphot
    # (This is why it's not a top-line import)
    from stsynphot import Vega, band # noqa
    if grism == 'G141':
        if detector:
            warn("WFC3's G141 grism does not have multiple detectors. Ignoring detector argument",
                 RuntimeWarning)
        bandpass = band('wfc3,ir,g141')
    elif grism == 'G102':
        if detector:
            warn("WFC3's G102 grism does not have multiple detectors. Ignoring detector argument",
                 RuntimeWarning)
        bandpass = band('wfc3,ir,g102')
    elif grism == 'G280':
        bandpass = band(f'wfc3,uvis{detector},g280')
    elif grism == 'G800L':
        bandpass = band(f'acs,wfc{detector},g800l')
    else:
        raise ValueError(f"Unrecognized grism: {grism}. Valid grisms: G141, G102, G280, G800L")

    spectrum = (bandpass * Vega).to_spectrum1d()

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
