from pathlib import Path

from .reference_file_generators._generate_specwcs import create_grism_specwcs
from .reference_file_generators._generate_wavelengthrange import create_tsgrism_wavelengthrange
#from .reference_file_generators._generate_distortion import create_distortion

def create_reference_files(conffile, hst_grism, outpath=Path.cwd()):
    create_tsgrism_wavelengthrange(outname=str(Path(outpath) / "HST_wavelengthrange.asdf"))
    create_grism_specwcs(conffile=str(conffile), pupil=hst_grism, outname=str(Path(outpath) / "HST_specwcs.asdf"))
    #create_distortion(detector, apname, outname, subarr, exp_type)
