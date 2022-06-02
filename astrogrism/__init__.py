# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

# Expose subpackage API at package level.
from astrogrism.grism_observation import GrismObs # noqa
from astrogrism.extraction import * # noqa
