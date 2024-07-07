from astropy.constants import k_B
import numpy as np


class STAR:
    def __init__(self, mass):
        self.mass = mass
        self.properties()

    def properties(self):

        # Using Mass-Luminosity relation
        if self.mass < 0.43:
            self.luminosity = 0.23 * (self.mass ** 2.3)
        elif 0.43 <= self.mass < 2.0:
            self.luminosity = self.mass ** 4.0
        elif 2.0 <= self.mass < 55:
            self.luminosity = 1.4 * (self.mass ** 3.5)
        elif 55 <= self.mass:
            self.luminosity = 32000 * self.mass


        # Using stefan-boltzmann law
        self.temperature = (self.luminosity / 4 * np.pi * k_B) ** 0.25

        if self.temperature > 30000:
            self.spectral_type = 'O'
        elif self.temperature > 10000:
            self.spectral_type = 'B'
        elif self.temperature > 7500:
            self.spectral_type = 'A'
        elif self.temperature > 6000:
            self.spectral_type = 'F'
        elif self.temperature > 5200:
            self.spectral_type = 'G'
        elif self.temperature > 3700:
            self.spectral_type = 'K'
        else:
            self.spectral_type = 'M'
