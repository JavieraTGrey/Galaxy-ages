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


# UPDATED IMF CLASS

class IMF_salpeter:
    default_mmin = 0.3
    default_mmax = 120

    def __init__(self, alpha=2.35, mmin=None, mmax=None):
  
        self._mmin = self.default_mmin if mmin is None else mmin
        self._mmax = self.default_mmax if mmax is None else mmax
        self.alpha = alpha
        self.slope = -alpha

        assert self._mmin < self._mmax, "Minimum mass must be less than maximum mass."
        assert self._mmin > 0, "Minimum mass must be greater than zero."
        assert self._mmin != -1, "Minimum mass cannot be -1."

    def pdf(self, x):

        is_in_range = (x >= self._mmin) & (x <= self._mmax)
        return x**self.slope * is_in_range

    def __call__(self, m):

        integral, _ = scipy.integrate.quad(self.pdf, self._mmin, self._mmax)
        return self.pdf(m) / integral

    def generate_masses(self, min_mass, max_mass, e=1):
        #salpeter
        cte = 1.35
        return ((min_mass**-cte) - (max_mass ** -cte)) * (e / cte)

    def generate_mass_range(self, delta):
        return np.arange(self._mmin, self._mmax + delta, delta)

    def imf_fraction(self, mass_list, e=1):

        tot = self.generate_masses(mass_list[0], mass_list[-1], e)
        frac = np.zeros_like(mass_list[:-1])
        n = len(mass_list) - 1

        for i in range(n):
            frac[i] = self.generate_masses(mass_list[i], mass_list[i + 1], e)

        return frac / tot

    def generate_imf(self, delta, e=1):

        mass_range = self.generate_mass_range(delta)
        fractions = self.imf_fraction(mass_range, e)
        return mass_range, fractions

