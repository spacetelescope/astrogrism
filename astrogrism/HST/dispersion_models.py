import numpy as np
from astropy.modeling import Model
from asdf.extension import Converter


class DISPXY_Model(Model):
    n_inputs = 3
    n_outputs = 1
    inputs = ("x", "y", "t")
    outputs = ("variable",)

    _tag = "tag:stsci.edu:grismstuff/dispxy_model-1.0.0"
    _name = "DISPXY_Model"

    def __init__(self, ematrix, offset, inv=False):
        self.ematrix = np.array(ematrix)
        self.inv = inv
        self.offset = offset

        if self.ematrix.shape == (2,):
            # Reshape to add second dimension for ematrix with only two values
            self.ematrix = np.reshape(self.ematrix, [2, 1])

        if len(self.ematrix.shape) > 1:
            if self.inv and self.ematrix.shape[0] > 2:
                # Can't invert these here, need to interpolate from the other direction
                raise ValueError("Can't invert higher order coefficient matrices")

        # This seems to be needed to avoid an error in calling the model
        self._model_set_axis = False

    # Note that in the inverse case, input "t" here is actually dx or dy
    def evaluate(self, x, y, t, t_op = "multiply"):

        e = self.ematrix
        offset = self.offset
        reshape_output = False

        # Handle reshaping of x and y to handle arrays if needed
        if x.ndim != y.ndim:
            raise ValueError("Input x and y arrays must have same dimensionality."
                             "2D arrays will be used as-is, 1D arrays will be broadcast"
                             "together. See documentation for further detail.")

        if x.ndim == 2:
            if x.shape != y.shape:
                raise ValueError("If x and y inputs are 2D their shapes must match")

        # If x and y are 1D arrays, assume we need to create a meshgrid
        elif x.ndim == 1 and x.shape != (1,) and y.shape != (1,):
            #mesh = np.meshgrid(x, y)
            #x = mesh[0]
            #y = mesh[1]
            if x.shape != y.shape:
                raise ValueError("If x and y are 1D arrays they must have the"
                                 "same shape. See documentation for instructions"
                                 "for dispersing a 2D region.")

        elif x.shape == (1,):
            if y.shape == (1,):
                # Explicitly include this just to not forget about this case
                pass
            if y.shape != (1,):
                x = np.full(y.shape, x[0])

        elif y.shape == (1,) and x.shape != (1,):
            y = np.full(x.shape, y[0])

        else:
            raise ValueError("Array input for x and y can only be 1 or 2 dimensional")


        # x and y should be the same shape at this point if at least one was an array
        if x.ndim == 2:
            reshape_output = True
            output_shape = x.shape
            x = x.flatten()
            y = y.flatten()

        const = np.full(x.shape, 1)

        coeffs = {1: np.array([1]),
                  6: np.array([const, x, y, x**2, x*y, y**2])}

        t_order = e.shape[0]
        if len(e.shape) == 1:
            c_order = 1
        else:
            c_order = e.shape[1]

        f = 0
        t = np.array(t)

        #print(f"e[0,:]: {e[0,:]}")
        #print(f"coeffs: {coeffs[c_order]}")
        #print(f"t: {t}")

        if self.inv:
            f = ((t + offset - np.dot(e[0, :], coeffs[c_order])) /
                  np.dot(e[1, :], coeffs[c_order]))
        else:
            for i in range(0, t_order):
                if t_op == "outer":
                    f += np.outer(t**i, (np.dot(e[i, :], coeffs[c_order])))
                elif t_op == "multiply":
                    f += np.multiply(t**i, (np.dot(e[i, :], coeffs[c_order])))

        if reshape_output:
            f = np.reshape(f, output_shape)

        return f


class DISPXY_ModelConverter(Converter):
    tags = ["tag:stsci.edu:grismstuff/dispxy_model-*"]
    types = [DISPXY_Model]

    def to_yaml_tree(self, obj, tags, ctx):
        # ASDF will know how to turn the nested lists into yaml properly
        return {"ematrix": obj.ematrix, "inverse_flag": obj.inv}

    def from_yaml_tree(self, node, tags, ctx):
        ematrix = node['ematrix']
        inverse_flag = node['inverse_flag']
        return DISPXY_Model(ematrix, inverse_flag)


class DISPXY_Extension():
    extension_uri = "asdf://stsci.edu/grismstuff/extensions/extension-1.0"
    converters = [DISPXY_ModelConverter()]
    tags = ["tag:stsci.edu:grismstuff/dispxy_model-1.0.0"]
