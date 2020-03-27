from find_function_max import find_function_max
from load_structures import load_structures
from read_spectrum_curve import Spectrum
import matplotlib.pyplot as plt


def calculate_intencity_scaling_factor(structure_name, bounds=(500, 4000)):
    structure = load_structures([structure_name])[structure_name]
    calc_maxes = find_function_max(structure.pm.get_spectrum_function().calculate, structure.freqs, bounds=bounds)
    spectrum = Spectrum(structure_name + '.csv')
    exp_int = []
    for freq_max in calc_maxes.peak_position:
        exp_int.append(spectrum.to_integrated_mol_absorptivity(freq_max, width=0.2))

    calc_maxes['exp_int'] = exp_int
    calc_maxes['sf'] = calc_maxes.intencity / calc_maxes.exp_int

    functional_groups = [' '.join(structure.get_fg_by_freq(freq)) for freq in calc_maxes.calc_freq]
    calc_maxes['functional_group'] = functional_groups
    return calc_maxes


if __name__ == '__main__':
    int_sf_table = calculate_intencity_scaling_factor('acetone')
    plt.plot(int_sf_table.calc_freq, int_sf_table.sf)

