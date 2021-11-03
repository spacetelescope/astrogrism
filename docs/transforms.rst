
.. _transforms:

Geometric Transforms
====================

The main functionality provided by Astrogrism is currently the geometric
transforms to go between the three reference frames relevant to grism
observations: the world (sky) frame, the direct (undispersed) image frame,
and the dispersed grism image frame. The inputs and outputs going to/from
each frame are listed below - note that the input needed to go from e.g. the
world to the direct detector frame is the same as the output of the transform
in the other direction, e.g. from the direct detector to the world frame.

Astropy ``Quantity`` input is currently accepted for wavelength, which will 
be converted automatically internally to the necessary units in that case. 
Currently, sky coordinates (right ascension and declination) must be input 
as decimal degrees. 


World
-----

Input/Output: ``(right ascension, declination, wavelength, order)``

The ``right ascension`` and ``declination`` are specified in degrees. The 
``wavelength`` units depend on instrument: Micron for ACS and WFC3 IR, 
Angstrom for WFC3 UVIS. If the wavelength is input as an astropy ``Quantity``, 
it will be automatically converted to the correct units, otherwise (if a float
is input) it must be in the correct units. The dispersion ``order`` is a 
unitless integer. 

Direct Detector
---------------

Input/Output: ``(x, y, wavelength, order)``

``x`` and ``y`` are unitless pixel positions on the detector; ``wavelength`` 
and ``order`` are the same as for the World frame.

Grism Detector
--------------

Input/Output: ``(x, y, x0, y0, order)``

All variables are unitless. ``x0`` and ``y0`` are the pixel positions on the direct detector.
``x`` and ``y`` are the corresponding dispersed pixel locations on the grism detector. 
