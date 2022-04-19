
.. _transforms:

Geometric Transforms
====================

The main functionality provided by Astrogrism is currently the geometric
transforms to convert between the three reference frames relevant to grism
observations: the world (sky) frame, the direct (undispersed) image frame,
and the dispersed grism image frame. The inputs and outputs going to/from
each frame are listed below - note that the input needed to go from e.g. the
world to the direct detector frame is the same as the output of the transform
in the other direction, e.g. from the direct detector to the world frame.

Astropy ``Quantity`` input is currently accepted for wavelength, which will 
be converted automatically internally to the necessary units in that case. 
Currently, sky coordinates (right ascension and declination) must be input 
as decimal degrees. 

Note that in the following descriptions, "Input/Output" describes both the
inputs to the geometric transform when transforming from that frame to another,
and the outputs when transforming from another frame to the one being described.


World Coordinate Frame
----------------------

Input/Output: ``(right ascension, declination, wavelength, order)``

The ``right ascension`` and ``declination`` are specified in degrees. The 
``wavelength`` units depend on instrument: Micron for ACS and WFC3 IR, 
Angstrom for WFC3 UVIS. If the wavelength is input as an astropy ``Quantity``, 
it will be automatically converted to the correct units, otherwise (if a float
is input) it must be in the correct units. The dispersion ``order`` is a 
unitless integer. 

Direct Detector Frame
---------------------

Input/Output: ``(x, y, wavelength, order)``

``x`` and ``y`` are unitless pixel positions on the detector; ``wavelength`` 
and ``order`` are the same as for the World frame.

Grism Detector Frame
--------------------

Input/Output: ``(x, y, x0, y0, order)``

All variables are unitless. ``x0`` and ``y0`` are the pixel positions on the direct detector.
``x`` and ``y`` are the corresponding dispersed pixel locations on the grism detector. 


Array Inputs
------------

Currently, the geometric transforms can take array inputs for either the 
spatial coordinates or wavelength values, but not both at the same time. For
example, the following would be valid uses::

    from astrogrism import GrismObs
    import astropy.units as u
    g_obs = GrismObs("sample_grism_flt.fits")
    detector_to_grism = g_obs.geometric_transforms["CCD1"].get_transform("detector", "grism_detector")

    # Array inputs for spatial (pixel) coordinates
    detector_to_grism([1024, 1030, 1036], [2048.0, 2050, 2052], .7*u.um, 1.0)

    # Array input for wavelength
    detector_to_grism(1024, 2048.0, [.7, .75, .8]*u.um, 1.0)

However, the following does not currently work::

    detector_to_grism([1024, 1030, 1036], [2048.0, 2050, 2052], [.7, .75, .8]*u.um, 1.0)

Using array inputs for wavelength is currently recommended over using array
inputs for the spatial coordinates, because the wavelength handling is properly
vectorized on the backend and has better performance compared to the array
handling for spatial coordinates. Using array inputs for spatial coordinates
currently only gives a factor of 2 increase in performance over looping, rather
than an order of magnitude or more as might be expected.