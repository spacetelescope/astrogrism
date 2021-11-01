
.. _quickstart:

Quickstart
==========

Once installed, ``astrogrism`` can be imported and used in an interactive Python
session or in a Jupyter notebook. The main class exposing used to read in grism
observation FITS files and expose the functionality is the ``astrogrism.GrismObs`` 
class.

To load your data, simply run the following in a Python interpreter or Jupyter
notebook cell::
    
    from astrogrism.grism_observation import GrismObs
    g_obs = GrismObs("sample_file.fits")

This object makes available the tranforms between the ``world``, ``detector`` 
(i.e. the undispersed direct image), and ``grism_detector`` frames. To get 
the transform between, for example, the world and grism detector frames, you 
can call::

    world_to_grism = g_obs.geometric_transforms.get_transform("world", "grism_detector")


