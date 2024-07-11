from astropy.constants import sigma_sb, L_sun, M_sun, R_sun
from astropy.io import fits
import astropy.units as u
from io import BytesIO
import numpy as np


class STAR:

    def __init__(self, mass, born, spectra):

        self.spectra = spectra
        self.mass = mass*M_sun
        self.born = born * u.Gyr
        self.stage = self.born
        self.branch = 'MS'
        self.properties()
        self.get_star_spectrum()

        self.wavelength = self.wavelengths[self.wavelengths < 10000]
        self.spectrum = self.spectrum[self.wavelengths < 10000]
        self.flux = self.flux[self.wavelengths < 10000]

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

        if (self.mass/M_sun) > 16:
            self.spectral_type = 'O'
        elif (self.mass/M_sun) > 2.1:
            self.spectral_type = 'B'
        elif (self.mass/M_sun) > 1.4:
            self.spectral_type = 'A'
        elif (self.mass/M_sun) > 1.04:
            self.spectral_type = 'F'
        elif (self.mass/M_sun) > 0.8:
            self.spectral_type = 'G'
        elif (self.mass/M_sun) > 0.45:
            self.spectral_type = 'K'
        else:
            self.spectral_type = 'M'

    def read_spectrum_from_fits(self, data):
        with fits.open(BytesIO(data)) as hdul:
            header = hdul[0].header
            spectrum_data = hdul[0].data

            wavelengths = np.arange(header['CRVAL1'],
                                    header['CRVAL1'] + header['CDELT1'] * len(spectrum_data),
                                    header['CDELT1'])
            return wavelengths[:len(spectrum_data)], spectrum_data

    def get_star_spectrum(self):
        new_type = [type for type in list(self.spectra) if type.startswith(self.spectral_type)]
        new_spectral_type = np.random.choice(new_type)
        if new_spectral_type in self.spectra:

            data = self.spectra[new_spectral_type]
            self.wavelengths, self.spectrum = self.read_spectrum_from_fits(data)

            self.total_lum = np.trapz(self.spectrum, x=self.wavelengths)

            self.flux = self.spectrum / self.total_lum

        else:
            raise ValueError(f"Spectral type {new_spectral_type} not available in spectra.")

    def plot(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.wavelength,
                    self.flux,
                    color='mediumslateblue')
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Normalized Flux (ergs^-1Å^-1)') # ergs−1A−1
        plt.title(f'Normalized Spectral Data - Type {self.spectral_type}')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.show()

    def update(self, t):
        # T es edad del universo!

        life = t - self.born
        HB_life = 1 *u.Gyr + 1*u.Myr

        if life.value < 0:
            self.branch = 'Unborn Star'

        elif (self.mass/M_sun) < 5:
            if life > self.t_ms:
                if life - self.t_ms < 1 * u.Gyr:
                    if self.branch == 'MS':
                        self.branch = 'RG'

                elif life - self.t_ms < HB_life:
                    self.branch = 'HB'
            
                else:
                    self.branch = 'Dead'
                    self.spectrum = np.zeros(self.spectrum.shape)
        else:
            if life > self.t_ms:
                self.branch = "Dead"
                self.spectrum = np.zeros(self.spectrum.shape)  
