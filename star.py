import numpy as np
from astropy.constants import sigma_sb, L_sun, M_sun, R_sun
import astropy.units as u


class STAR:

    def __init__(self, mass, born):

        self.mass = mass*M_sun
        self.born = born * u.Gyr
        self.stage = self.born
        self.branch = 'MS'
        self.properties()

    def properties(self):

        # Using Mass-Luminosity relation
        if self.mass/M_sun < 0.43:
            self.luminosity = 0.23 * ((self.mass/M_sun) ** 2.3) * L_sun
        elif 0.43 <= self.mass/M_sun < 2.0:
            self.luminosity = ((self.mass/M_sun) ** 4.0) * L_sun
        elif 2.0 <= self.mass/M_sun < 55:
            self.luminosity = 1.4 * L_sun * ((self.mass/M_sun) ** 3.5)
        elif 55 <= self.mass/M_sun:
            self.luminosity = 32000 * (self.mass/M_sun) * L_sun

        # Using stafan-boltzmann law
        sigma = sigma_sb.to(u.W / (u.m**2 * u.K**4))
        self.radii = R_sun * (self.mass / M_sun)**0.8
        denominator = (4 * np.pi * (self.radii**2) * sigma)
        self.temperature = ((self.luminosity / denominator) ** 0.25)
        
        if self.temperature.value >= 33000:
            self.spectral_type = 'O'
        elif self.temperature.value >= 10000:
            self.spectral_type = 'B'
        elif self.temperature.value >= 7300:
            self.spectral_type = 'A'
        elif self.temperature.value >= 6000:
            self.spectral_type = 'F'
        elif self.temperature.value >= 5300:
            self.spectral_type = 'G'
        elif self.temperature.value >= 3900:
            self.spectral_type = 'K'
        else:
            self.spectral_type = 'M'

    def update(self, t):
        # T es edad del universo!
        cons = (10 * u.Gyr)*(M_sun**3)
        self.t_ms = cons / (self.mass**3)

        life = t - self.stage

        if (self.mass/M_sun) < 5:
            if life > self.t_ms:
                if self.branch == 'MS':
                    self.branch == 'RG'
        else:
            if life > self.t_ms:
                self.branch = "Dead"
