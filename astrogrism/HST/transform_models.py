import warnings

import numpy as np
import astropy.units as u
from astropy.modeling.core import Model
from astropy.modeling.models import (Rotation2D, Identity, Mapping, Tabular1D, Const1D)
from astropy.modeling.models import math as astmath


class AstrogrismForwardGrismDispersion(Model):
    """Return the transform from grism to image for the given spectral order.

    Parameters
    ----------
    orders : list [int]
        List of orders which are available

    lmodels : list [astropy.modeling.Model]
        List of models which govern the wavelength solutions for each order

    xmodels : list [astropy.modeling.Model]
        List of models which govern the x solutions for each order

    ymodels : list [astropy.modeling.Model]
        List of models which givern the y solutions for each order

    Returns
    -------
    x, y, wavelength, order in the grism image for the pixel at x0,y0 that was
    specified as input using the input delta pix for the specified order

    Notes
    -----
    Based on the JWST NIRISS transform model code, see
    https://github.com/spacetelescope/jwst/blob/master/jwst/transforms/models.py
    """

    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 5
    n_outputs = 4

    def __init__(self, orders, lmodels=None, xmodels=None, ymodels=None,
                 theta=0., name=None, meta=None, l_unit="micron",
                 x_offset=0, y_offset=0):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.xmodels = xmodels
        self.ymodels = ymodels
        self.lmodels = lmodels
        self.theta = theta
        self.orders = orders
        self.l_unit = l_unit
        self.x_offset = x_offset
        self.y_offset = y_offset
        meta = {"orders": orders}
        if name is None:
            name = 'astrogrism_forward_grism_dispersion'
        super(AstrogrismForwardGrismDispersion, self).__init__(name=name,
                                                               meta=meta)
        # starts with the backwards pixel and calculates the forward pixel
        self.inputs = ("x", "y", "x0", "y0", "order")
        self.outputs = ("x", "y", "wavelength", "order")

    def evaluate(self, x, y, x0, y0, order):
        """Return the valid pixel(s) and wavelengths given x0, y0, x, y, order
        Parameters
        ----------
        x0: int,float,list
            Source object x-center
        y0: int,float,list
            Source object y-center
        x :  int,float,list
            Input x location
        y :  int,float,list
            Input y location
        order : int
            Spectral order to use

        Returns
        -------
        x0, y0, lambda, order in the direct image for the pixel that was
        specified as input using the wavelength l and spectral order

        Notes
        -----
        I kept the possibility of having a rotation like NIRISS, although I
        don't know if there is a use case for it for HST instruments.

        The two `flatten` lines may need to be uncommented if we want to use
        this for array input.
        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except AttributeError:
            iorder = self._order_mapping[order]
        except KeyError:
            raise ValueError("Specified order is not available")

        t = np.linspace(0, 1, 10)

        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]

        # For subarrays, we need to get use the location on the full array
        # to calculate offsets
        dx = xmodel.evaluate(x0 + self.x_offset, y0+self.y_offset, t, t_op="outer")
        dy = ymodel.evaluate(x0 + self.x_offset, y0+self.y_offset, t, t_op="outer")

        if len(dx.shape) == 2:
            dx = dx[:, 0]
            dy = dy[:, 0]

        if self.theta != 0.0:
            rotate = Rotation2D(self.theta)
            dx, dy = rotate(dx, dy)

        so = np.argsort(dx)
        tab = Tabular1D(dx[so], t[so], bounds_error=False, fill_value=None)

        dxr = astmath.SubtractUfunc()

        # Need to build this compound model differently depending on lmodel inputs
        if lmodel.n_inputs == 1:
            wavelength = dxr | tab | lmodel
            model = Mapping((2, 3, 0, 2, 4)) | (Const1D(x0) & Const1D(y0) &
                                                wavelength & Const1D(order))
        elif lmodel.n_inputs == 3:
            wavelength = Identity(2) & dxr | Identity(2) & tab | lmodel
            model = Mapping((2, 3, 0, 1, 0, 2, 4)) | (Const1D(x0) & Const1D(y0) &
                                                      wavelength & Const1D(order))

        x_out, y_out, l_out, order_out = model(x, y, x0, y0, order)

        # In this case we had a single wavelength with multiple pixel locations
        if isinstance(l_out, np.ndarray):
            if x0.shape == l_out.shape or y0.shape == l_out.shape:
                l_out = l_out[0]

        return x_out, y_out, l_out*u.Unit(self.l_unit), order_out


class AstrogrismBackwardGrismDispersion(Model):
    """Return the dispersed pixel(s) given center x, y, lambda, and order

    Parameters
    ----------
    xmodels : list[tuple]
        The list of tuple(models) for the polynomial model in x
    ymodels : list[tuple]
        The list of tuple(models) for the polynomial model in y
    lmodels : list
        The list of models for the polynomial model in l
    orders : list
        The list of orders which are available to the model
    theta : float
        Angle [deg] - defines the NIRISS filter wheel position

    Returns
    -------
    x, y, x0, y0, order in the grism image for the pixel at x0,y0 that was
    specified as input using the wavelength l for the specified order

    Notes
    -----
    Based on the JWST transform model code, see
    https://github.com/spacetelescope/jwst/blob/master/jwst/transforms/models.py
    """
    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 4
    n_outputs = 5

    def __init__(self, orders, lmodels=None, xmodels=None, ymodels=None,
                 theta=None, name=None, meta=None, interpolate_t=False,
                 l_unit="micron", x_offset=0, y_offset=0):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.xmodels = xmodels
        # TODO: Raise a warning if no inverse transform is possible (for example
        # the current state of UVIS)
        self.ymodels = ymodels
        self.lmodels = lmodels
        self.orders = orders
        self.theta = theta
        self.interpolate_t = interpolate_t
        self.l_unit = l_unit
        self.x_offset = x_offset
        self.y_offset = y_offset
        meta = {"orders": orders}
        if name is None:
            name = 'astrogrism_backward_grism_dispersion'
        super(AstrogrismBackwardGrismDispersion, self).__init__(name=name,
                                                                meta=meta)
        self.inputs = ("x", "y", "wavelength", "order")
        self.outputs = ("x", "y", "x0", "y0", "order")

    def evaluate(self, x, y, wavelength, order):
        """Return the dispersed pixel(s) given center x, y, lam and order
        Parameters
        ----------
        x :  int,float
            Input x location on the direct image
        y :  int,float
            Input y location on the direct image
        wavelength : float
            Wavelength to disperse
        order : list
            The order to use

        Returns
        -------
        x, y in the grism image for the pixel at x0, y0 that was
        specified as input using the wavelength and order specified

        Notes
        -----
        I kept the potential for rotation from NIRISS, unsure if it's actually
        needed/useful for HST instruments. Original note:

        There's spatial dependence for NIRISS so the forward transform
        dependes on x,y as well as the filter wheel rotation. Theta is
        usu. taken to be the different between fwcpos_ref in the specwcs
        reference file and fwcpos from the input image.
        """
        if self.ymodels is None:
            return None
        if np.any(wavelength < 0):
            raise ValueError("Wavelength should be greater than zero")

        # Check if wavelength is in unit expected (if Quantity)
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(self.l_unit).value
        else:
            warnings.warn(f"Assuming input wavelength is in {self.l_unit}. To "
                          "specify wavelength unit, input an astropy Quantity.")
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except AttributeError:
            iorder = self._order_mapping[order]
        except KeyError:
            raise ValueError("Specified order is not available")

        if self.interpolate_t:
            # If the displ coefficients are too complex to invert, have to interpolate t
            t = np.linspace(0, 1, 1000)
            if self.lmodels[iorder].n_inputs == 1:
                l = self.lmodels[iorder].evaluate(t)
            elif self.lmodels[iorder].n_inputs == 3:
                l = self.lmodels[iorder].evaluate(x+self.x_offset,
                                                  y+self.y_offset,
                                                  t, t_op="outer")

            # Loop to account for arrays
            t_fit = []
            for i in range(l.shape[1]):
                so = np.argsort(l[:, i])
                tab = Tabular1D(l[so, i], t[so], bounds_error=False, fill_value=None)
                t_fit.append(tab(wavelength))
            # tab = Tabular1D(l, t, bounds_error=False, fill_value=None)
            # t = tab(wavelength)
            t = np.array(t_fit).flatten()
        else:
            if self.lmodels[iorder].n_inputs == 1:
                t = self.lmodels[iorder](wavelength)
            elif self.lmodels[iorder].n_inputs == 3:
                t = self.lmodels[iorder](x+self.x_offset, y+self.y_offset, wavelength)

        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]

        dx = xmodel.evaluate(x+self.x_offset, y+self.y_offset, t)
        dy = ymodel.evaluate(x+self.x_offset, y+self.y_offset, t)

        # Rotate by theta. Do we need to account for subarray offsets here?
        if self.theta != 0.0:
            rotate = Rotation2D(self.theta)
            dx, dy = rotate(dx, dy)

        return (x+dx, y+dy, x, y, order)
