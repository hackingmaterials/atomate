# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
from argparse import ArgumentParser

import numpy as np


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def raman_intensity(R):
    """
    Scattering intensity
    
    Args:
        R (3x3 numpy array): Susceptibility tensor

    Returns:
        float: scattering intensity according to Eq 7, Phys Rev B 73, 104304 2006
    """
    return(4 * (R[0, 0]**2 + R[1, 1]**2 + R[2, 2]**2)
           + 7 * (R[0, 1]**2 + R[0, 2]**2 + R[1, 2]**2)
           + R[0, 0]*R[1, 1] + R[0, 0]*R[2, 2] + R[1, 1]*R[2, 2])


def lorentzian(w, wj, delta_w=4.0):
    """ 
    delta function approximation
    
    Args:
        w (float/np.ndarray): frequency
        wj (float): frequency where the delta function will be centered
        delta_w (float): peak width, default= 4 cm^-1(from the above paper)

    Returns:
        broadened delta function
    """
    # (Eq 8, Phys Rev B 73, 104304 2006)
    return 4.0 / np.pi * (w**2 * delta_w) / ((w**2 - wj**2)**2 + 4*delta_w**2 * w**2)
    # wolfram
    #return 1./np.pi * (delta_w/2)/ ( (w - wj)**2 + (delta_w/2)**2 )

    
def get_freq_intensities(input_file):
    """
    Return frequencies and normalized intensities from the json input file

    Args:
        input_file (str): path to the json file containing the raman tensors

    Returns:
        (numpy array, numpy array): frequencies, intensities
    """
    raman = json.load(open(input_file))
    raman_tensor_dict = raman["raman_tensor"]
    # positive(and zero freq acoustic modes) eigenvals are unstable, so skipped
    intensities = []
    freq = []
    for i, f in enumerate(raman["normalmodes"]["eigenvals"]):
        if f < 0:
            freq.append(raman["frequencies"][i])
            intensities.append(raman_intensity(np.array(raman_tensor_dict[str(i)])))
    return np.array(freq), np.array(intensities) / np.max(intensities)  # normalize


def plot_spectrum(freq, intensities, lorentzian_width):
    """
    Plot the Raman spectrum.

    Args:
       freq (numpy array): frequencies
       intensities (numpy array): intensities computed from the raman tensors
       lorentzian_width (float): width of the lorentzian delta function
    """
    import matplotlib.pyplot as plt
    
    # x-axis: wavenumber
    # multiply the wavenumbers by 2*pi    
    wavenumbers = np.linspace(0, np.max(freq)+10, 100)*2*np.pi
    # y-axis: intensity modulated by the lorentzian delta function
    broadened_intensities = np.zeros_like(wavenumbers)
    for i, f in enumerate(freq):
        broadened_intensities += intensities[i] * lorentzian(wavenumbers, f * 2 * np.pi, lorentzian_width)
    plt.plot(wavenumbers, broadened_intensities)
    plt.savefig("Raman_spectrum.pdf")


def run():
    parser = ArgumentParser(description='Example: plot_raman.py -i raman.json -w 4')
    parser.add_argument('-i', '--input_file', help='path to json input file')
    parser.add_argument('-w', '--width', help='lorentzian delta width in cm^-1', type=float)
    args = parser.parse_args()

    # input_file = "FeO_raman.json"
    # lorentzian_width = 4 # cm^-1, Phys Rev B 73, 104304 2006
    frequencies, intensities = get_freq_intensities(args.input_file)
    print("Frequencies(cm^-1) and Normalized intensities: \n {} \n{}".format(frequencies, intensities))
    plot_spectrum(frequencies, intensities, args.width)


if __name__ == "__main__":
    run()
