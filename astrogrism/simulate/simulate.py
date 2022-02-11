from math import floor
from pathlib import Path

from astrogrism import GrismObs
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import numpy as np

from specutils import Spectrum1D
import stsynphot as stsyn


def _fake_spectrum():
    spec_axis = np.linspace(6000, 8000, 10) * u.AA
    flux = (np.random.randn(len(spec_axis.value)) +
            10*np.exp(-0.001*(spec_axis.value-6563)**2) +
            spec_axis.value/500) * u.Jy

    return Spectrum1D(spectral_axis=spec_axis, flux=flux)


def _generate_simulation_spectrum(grism):
    return _fake_spectrum()


def _create_simulation_cube(spectrum, shape, wcs):
    data = np.tile(spectrum.flux.value * spectrum.flux.unit, (shape[0], shape[1], 1))
    print(data.shape)
    return Spectrum1D(flux=data, spectral_axis=spectrum.spectral_axis, wcs=wcs)


def _disperse_simulation_cube(grism, wide_field_image, spectrum_cube):
    #if Path(wide_field_image).is_file:

    grismobs = GrismObs(wide_field_image)

    shape = wide_field_image['SCI'].data.shape
    simulated_data = np.zeroes(shape)
    wcs = WCS(wide_field_image['SCI'].header)
    # For each pixel in the science image, we need to disperse it's spectrum
    for horizontal in range(0, shape[0]):
        for vertical in range(0, shape[1]):
            # Get the flux of the science pixel; we'll need to scale the spectrum to this brightness
            data_flux = wide_field_image['SCI'].data[horizontal][vertical]
            # Calculate the pixels to the RA/Dec sky coordinates
            sky_coord = wcs.pixel_to_world(horizontal, vertical)
            # Lookup the associated RA/Dec's spectrum in the simulated spectrum cube
            spectrum_coord = spectrum_cube.wcs.world_to_pixel(sky_coord.ra, sky_coord.dec)
            #TBF: Slice the clube using the RA/Dec's pixel equivalents to get the actual spectrum
            spectrum = spectrum_cube[spectrum_coord[0]][spectrum_coord[1]]

            # For each Wavelength in the spectrum, calculate where, in pixels, that wavelength would fall on the detector
            image2grism = grismobs.geometric_transforms.get_transform('detector', 'grism_detector')
            for wavelength in spectrum:
                spectrum_flux = spectrum[wavelength]
                #TBF: What is xcenter/ycenter in this context??
                dispersion = image2grism.evaluate(x_center, y_center, wavelength, 1)
                x = (dispersion[0])
                y = (dispersion[1])

                # If the dispersed position of the wavelength is inside the bounds of the image, write the spectrum
                if x in range(0, simulated_data.shape[0]) and y in range(0, simulated_data.shape[1]):
                    # Scale the flux of the spectrum to the brightness of the original pixel
                    # NOTE: Is floor the right approach to determine which pixel to write to?
                    # TBF: How to properly index a 3D numpy array?
                    simulated_data[(floor(x),floor(y))] = data_flux * spectrum_flux

    return simulated_data


def simulate_grism(grism, wide_field_image):
    spectrum = _generate_simulation_spectrum(grism)
    tiled_spectrum = _create_simulation_cube(spectrum,
                                             wide_field_image['SCI'].data.shape,
                                             WCS(wide_field_image['SCI'].header))
    _disperse_simulation_cube(grism, wide_field_image, tiled_spectrum)
