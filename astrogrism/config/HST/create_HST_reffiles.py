from pathlib import Path

from astrogrism.config.HST.reference_file_generators._generate_specwcs import create_grism_specwcs  # noqa
from astrogrism.config.HST.reference_file_generators._generate_wavelengthrange import create_tsgrism_wavelengthrange  # noqa
# from astrogrism.config.HST.reference_file_generators.reference_file_generators._generate_distortion import create_distortion  # noqa


def create_reference_files(conffile, hst_grism, outpath=Path.cwd(), outbasename="HST"):
    create_tsgrism_wavelengthrange(outname=str(Path(outpath) / (
                                                                str(outbasename) +
                                                                "_wavelengthrange.asdf")
                                               )
                                   )
    create_grism_specwcs(conffile=str(conffile),
                         pupil=hst_grism,
                         outname=str(Path(outpath) / (str(outbasename) + "_specwcs.asdf")))
    # create_distortion(detector, apname, outname, subarr, exp_type)
