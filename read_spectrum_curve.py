import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.constants import gas_constant
from shutil import copyfile


class Spectrum:
    def __init__(self, name, path='spectra_csv/absorbtivity', data=None):
        self.full_path = path
        self.data = data
        if path is not None:
            self.full_path = path + '/' + name
            self.data = self._read_data()
        self.x = self.data.x
        self.y = self.data.y

    def __mul__(self, other):
        new_data = pd.DataFrame({'x': self.data.x, 'y': other * self.data.y})
        return Spectrum(name=None, path=None, data=new_data)

    def __rmul__(self, other):
        return self.__mul__(other)

    def _read_data(self):
        data = pd.read_csv(self.full_path, ';')
        data.columns = ['x', 'y']
        data.astype('float32')
        return data

    # def _get_full_width_at_half_maximum(self, max_point):

    def to_integrated_mol_absorptivity(self, freq, width=0.05, temperature=296):
        p = 1.013e5  # normal pressure in Pascals
        constant = gas_constant * temperature / p * 1e10
        corrected_intencity = self.data.y * constant

        spectrum_function = interp1d(self.data.x, corrected_intencity)

        delta = width * freq
        nu1 = freq - delta
        nu2 = freq + delta
        x_min = max(min(self.data.x), nu1)
        x_max = min(max(self.data.x), nu2)

        absorptivity = quad(spectrum_function, x_min, x_max)[0] * 1e-5

        return absorptivity

    def to_integrated_mol_absorptivity_full_spectra(self, width=0.05, number_of_points=100, temperature=296):
        step = int(len(self.x) / number_of_points)
        points = [self.x[step * i] for i in range(number_of_points)]
        absorptivity = [self.to_integrated_mol_absorptivity(point, temperature=temperature, width=width) for point in points]
        return Spectrum(name=None, path=None, data=pd.DataFrame({'x': points, 'y': absorptivity}))

    def plot(self, ax=plt):
        ax.plot(self.x, self.y)

    def show(self, ax=plt):
        self.plot(ax=plt)
        plt.show()


def read_spectrum_curve(name, path='spectra_csv'):
    return Spectrum(name, path)


def convert_to_absorbtivity(name, length, p1, p_total, path='spectra_csv/initial_spectra', units='absorbance'):
    initial_spectrum = Spectrum(name, path=path)
    if units == 'absorbance':
        absorbance = initial_spectrum.y
    elif units == 'transmitance':
        absorbance = - np.log10(initial_spectrum.y)
    else:
        print('Unknown units', units)
        return
    c = p1 / p_total * 1e6
    absorbtivity = absorbance / c / length / 1e-2

    new_path = 'spectra_csv/absorbtivity/' + name
    with open(new_path, 'w') as input_file:
        for i in range(len(absorbtivity)):
            if absorbtivity[i] != float('inf'):
                string = str(initial_spectrum.x[i]) + ';' + str(absorbtivity[i]) + str('\n')
                input_file.write(string)

    print('New spectrum was placed at', new_path)
    return new_path


if __name__ == '__main__':
    Spectrum('salicylic_acid.CSV').to_integrated_mol_absorptivity_full_spectra(width=0.05).show()
    # Spectrum('anthraflavic_acid.CSV', path='spectra_csv/initial_spectra').show()

    # convert_to_absorbtivity('salicylic_acid.CSV', length=1, p1=30, p_total=600, units='absorbance')
