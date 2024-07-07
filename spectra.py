from astropy.io import fits
from ftplib import FTP
from io import BytesIO
import numpy as np


# Function to download spectrum from FTP
def download_spectrum_ftp(url):
    from urllib.parse import urlparse

    parsed_url = urlparse(url)
    ftp_host = parsed_url.hostname
    ftp_file_path = parsed_url.path

    ftp = FTP(ftp_host)
    ftp.login()

    file_data = BytesIO()
    ftp.retrbinary(f"RETR {ftp_file_path}", file_data.write)
    ftp.quit()

    file_data.seek(0)
    return file_data.read()


# Function to open FITS data
def open_fits_from_data(data):
    with fits.open(BytesIO(data)) as hdul:
        header = hdul[0].header
        data = hdul[0].data
        return header, data


# Function to calculate the total luminosity of the spectrum
def calculate_total_luminosity(wavelength, spectrum):
    luminosity = np.trapz(spectrum, x=wavelength)
    return luminosity


# Function to normalize a stellar spectrum relative to the total luminosity
def normalize_spectrum_to_luminosity(header, data):
    crpix1 = header['CRPIX1']  # Reference pixel
    crval1 = header['CRVAL1']  # Value of the reference pixel
    cdelt1 = header['CDELT1']  # Increment per pixel

    # Calculate the wavelength range
    wavelength = crval1 + (np.arange(len(data)) - crpix1 + 1) * cdelt1

    # Calculate the total luminosity of the spectrum
    total_luminosity = calculate_total_luminosity(wavelength, data)

    # Normalize the spectrum by dividing by the total luminosity
    normalized_spectrum = data / total_luminosity

    return wavelength, normalized_spectrum, total_luminosity
