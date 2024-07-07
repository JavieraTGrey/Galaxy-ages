import numpy as np
import scipy


class IMF_salpeter():
    default_mmin = 0.3
    default_mmax = 120

    def __init__(self, alpha=2.35, mmin=None, mmax=None):

        self._mmin = self.default_mmin if mmin is None else mmin
        self._mmax = self.default_mmax if mmax is None else mmax

        self.alpha = alpha

        self.slope = -alpha
        assert (self._mmin < self._mmax)
        assert (self._mmin > 0)
        assert (self._mmin != -1)

    def pdf(self, x):

        is_is_range = (x >= self._mmin) * (x <= self._mmax)

        return x**self.slope * is_is_range

    def generate_masses(self, delta, sfr):
        mass_range = np.arange(self._mmin, self._mmax + delta, delta)
        masses = mass_range * sfr
        return masses

    def __call__(self, m):
        integral, _ = scipy.integrate.quad(self.pdf, self._mmin, self._mmax)
        assert (np.round(integral, 5) == 1.)

        return self.pdf(m)
