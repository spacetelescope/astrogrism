
.. _quickstart:

Quickstart
==========

Once installed, ``astrogrism`` can be imported and used in an interactive Python
session or in a Jupyter notebook. The main class used to read in grism
observation FITS files and expose the functionality is the ``astrogrism.GrismObs`` 
class.

To load your data, simply run the following in a Python interpreter or Jupyter
notebook cell::
    
    from astrogrism import GrismObs
    g_obs = GrismObs("sample_grism_flt.fits")

This object makes available the tranforms between the ``world``, ``detector`` 
(i.e. the undispersed direct image), and ``grism_detector`` frames. To get 
the transform between, for example, the world and grism detector frames, you 
can call::

    world_to_grism = g_obs.geometric_transforms.get_transform("world", "grism_detector")

For instruments with two chips (WFC3 UVIS and ACS), you must specify the chip for
which you want the transform, e.g.::
   
    world_to_grism = g_obs.geometric_transforms["CCD1"].get_transform("world", "grism_detector")

This transform would then allow you to calculate the coordinates on the dispered image
given a right ascension and declination (currently required to be in degrees), the
wavelength of interest, and the order of the dispersed spectrum::

    world_to_grism(264.10183, -32.90802, 0.7, 1.0)

This returns a tuple of values ``(x, y, x0, y0, order)``, where ``(x0, y0)`` denotes 
coordinates on the direct image, and ``(x, y)`` denotes the dispersed coordinates 
on the grism image. For more detail about the inputs and outputs of the geometric
transforms, see :ref:`transforms`.

Note that, in addition to the geometric transforms, the ``GrismObs`` object 
stores the contents of the input FITS file as an ``astropy`` HDUList in the 
``grism_image`` attribute. The direct image can also be loaded, e.g::

    g_obs = GrismObs("sample_grism_flt.fits", direct_image="sample_direct.fits")

and is stored in the ``direct_image`` attribute in the same way as the ``grism_image``::

    GrismObs.direct_image
    >>> [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x7f9240ca8880>, 
         <astropy.io.fits.hdu.image.ImageHDU object at 0x7f9201bcfd60>, ...]
    
    GrismObs.grism_image
    >>> [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x7f9240ca8880>,
         <astropy.io.fits.hdu.image.ImageHDU object at 0x7f9201bcfd60>, ...]
