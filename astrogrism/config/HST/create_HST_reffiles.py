from pathlib import Path

from astrogrism.config.HST.reference_file_generators._generate_specwcs import create_grism_specwcs  # noqa
from astrogrism.config.HST.reference_file_generators._generate_wavelengthrange import create_tsgrism_wavelengthrange  # noqa
# from astrogrism.config.HST.reference_file_generators.reference_file_generators._generate_distortion import create_distortion  # noqa


def create_reference_files(conffile, hst_grism, outpath=Path.cwd(), outbasename=None):
    """
    Generates the HST reference files required for the GrismObs class.
    Currently, wavelengthrange and specwcs are supported.
    Distortion addition is stil TODO

    Parameters
    ----------
    conffile : str or pathlib.PATH
        The path to the grism's specific GRISMCONF configuration file
    hst_grism : str
        The specific grism to generate the configuration files for
    outpath : str or pathlib.PATH
        Directory where the files should be output
    outbasename : str or pathlib.PATH
        The filename basename of the reference files to be named
    """

    if outbasename is None:
        outbasename = Path(conffile).name

    wavelengthrange_filename = str(Path(outpath) / (str(outbasename) + "_wavelengthrange.asdf"))
    create_tsgrism_wavelengthrange(outname=str(wavelengthrange_filename))

    specwcs_filename = str(Path(outpath) / (str(outbasename) + "_specwcs.asdf"))
    create_grism_specwcs(conffile=str(conffile), pupil=hst_grism, outname=str(specwcs_filename))

    # TODO: Implement distortion generation (non conf generator)
    # create_distortion(detector, apname, outname, subarr, exp_type)
