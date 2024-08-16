

from astropy.modeling.functional_models import Fittable2DModel
from astropy.modeling.parameters import Parameter
from astropy.units import UnitsError

import numpy as np

class MoffatElliptical2D(Fittable2DModel):
    """
    Two dimensional Moffat model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    y_0 : float
        y position of the maximum of the Moffat model.
    gamma1 : float
        Core width in 1st arbitrary(?) direction of the Moffat model
    gamma2 : float
        Core width in 2nd arbitrary(?) direction of the Moffat model
    phi : float
        Rotation angle of the Moffat model (rad)
    alpha : float
        Power index of the Moffat model.

    See Also
    --------
    Gaussian2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1, description="Amplitude (peak value) of the model")
    x_0 = Parameter(
        default=0, description="X position of the maximum of the Moffat model"
    )
    y_0 = Parameter(
        default=0, description="Y position of the maximum of the Moffat model"
    )
    gamma1 = Parameter(default=1, description="Core width in 1st arbitrary(?) direction of the Moffat model")
    gamma2 = Parameter(default=1, description="Core width in 1st arbitrary(?) direction of the Moffat model")
    phi = Parameter(default=0, description="Rotation angle of the Moffat model (rad)")
    alpha = Parameter(default=1, description="Power index of the Moffat model")

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        `this notebook by Yoonsoo Bach
        <https://nbviewer.jupyter.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat>`_.
        """
        return 2.0 * np.abs(self.gamma) * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, gamma1, gamma2, phi, alpha):
        """Two dimensional Moffat model function."""
        cosr = np.cos(phi)
        sinr = np.sin(phi)
        A = (cosr/gamma1)**2 + (sinr/gamma2)**2
        B = (sinr/gamma1)**2 + (cosr/gamma2)**2
        C = 2 * sinr * cosr * (gamma1**-2 - gamma2**-2)
        dx = (x - x_0)
        dy = (y - y_0)
        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        return amplitude / D**alpha

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma1, gamma2, phi, alpha):
        """Two dimensional Moffat model derivative with respect to parameters."""
        cosr = np.cos(phi)
        sinr = np.sin(phi)

        a1 = (cosr/gamma1)**2
        a2 = (sinr/gamma2)**2
        A = a1 + a2
        dA_dg1 = -2*a1/gamma1
        dA_dg2 = -2*a2/gamma2

        b1 = (sinr/gamma1)**2
        b2 = (cosr/gamma2)**2
        B = b1 + b2
        dB_dg1 = -2*b1/gamma1
        dB_dg2 = -2*b2/gamma2

        C = 2 * sinr * cosr * (gamma1**-2 - gamma2**-2)
        dC_dg1 = -4 * sinr * cosr / gamma1**3
        dC_dg2 = 4 * sinr * cosr / gamma2**3

        dA_dphi = -C
        dB_dphi = C
        dC_dphi = C * (cosr / sinr - sinr / cosr)  # = 2 * C / tan(2 * phi)

        dx = (x.ravel() - x_0)
        dy = (y.ravel() - y_0)

        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        dD_dx0 = -2*A*dx - C
        dD_dy0 = -2*B*dy - C
        dD_dg1 = dA_dg1*dx**2 + dB_dg1*dy**2 + dC_dg1*dx*dy
        dD_dg2 = dA_dg2*dx**2 + dB_dg2*dy**2 + dC_dg2*dx*dy
        dD_dphi = dA_dphi*dx**2 + dB_dphi*dy**2 + dC_dphi*dx*dy

        I = amplitude / D**alpha
        f = -alpha / D

        return [f * I * dD_dx0,
                f * I * dD_dy0,
                I/amplitude,
                f * I * dD_dg1,
                f * I * dD_dg2,
                f * I * dD_dphi,
                -I * np.log(D)]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        else:
            return {
                self.inputs[0]: self.x_0.input_unit,
                self.inputs[1]: self.y_0.input_unit,
            }

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "y_0": inputs_unit[self.inputs[0]],
            "gamma1": inputs_unit[self.inputs[0]],
            "gamma2": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }
        
        
        
class Moffat2D(Fittable2DModel):
    """
    Two dimensional Moffat model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    y_0 : float
        y position of the maximum of the Moffat model.
    gamma : float
        Core width of the Moffat model.
    alpha : float
        Power index of the Moffat model.

    See Also
    --------
    Gaussian2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1, description="Amplitude (peak value) of the model")
    x_0 = Parameter(
        default=0, description="X position of the maximum of the Moffat model"
    )
    y_0 = Parameter(
        default=0, description="Y position of the maximum of the Moffat model"
    )
    gamma = Parameter(default=1, description="Core width of the Moffat model")
    alpha = Parameter(default=1, description="Power index of the Moffat model")

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        `this notebook by Yoonsoo Bach
        <https://nbviewer.jupyter.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat>`_.
        """
        return 2.0 * np.abs(self.gamma) * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Moffat model function."""
        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma**2
        return amplitude * (1 + rr_gg) ** (-alpha)

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Moffat model derivative with respect to parameters."""
        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma**2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = 2 * amplitude * alpha * d_A * (x - x_0) / (gamma**2 * (1 + rr_gg))
        d_y_0 = 2 * amplitude * alpha * d_A * (y - y_0) / (gamma**2 * (1 + rr_gg))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * rr_gg / (gamma * (1 + rr_gg))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        else:
            return {
                self.inputs[0]: self.x_0.input_unit,
                self.inputs[1]: self.y_0.input_unit,
            }

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "y_0": inputs_unit[self.inputs[0]],
            "gamma": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }