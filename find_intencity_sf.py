from load_structures import load_structures
from find_function_max import find_function_max
import lmfit
import numpy as np
import rampy as rp
import pandas as pd
from scipy.constants import gas_constant
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('TkAgg')


def generate_parameters(norm_constant, intensities, mus, sigmas=None):
    params = lmfit.Parameters()
    for i in range(len(intensities)):
#                         (Name, Value, Vary, Min, Max, Expr)
        intensity_param = ('i' + str(i), intensities[i] / norm_constant, True, 0, None, None)
        mu_param = ('mu' + str(i), mus[i], True, mus[i] - 56, mus[i] + 56, None)
        sigma_param = ('sigma' + str(i), 30, True, 10, 80, None)
        params.add_many(intensity_param, mu_param, sigma_param)
    return params


def residual(pars, x, data=None, eps=None):  # Function definition
    # unpack parameters, extract .value attribute for each parameter
    # Using the Gaussian model function from rampy
    peaks = []
    for i in range(int(len(pars) / 3)):
        i_i = pars['i' + str(i)].value
        mu_i = pars['mu' + str(i)].value
        sigma_i = pars['sigma' + str(i)].value
        peak_i = rp.gaussian(x, i_i, mu_i, sigma_i)
        peaks.append(peak_i)
    model = sum(peaks)  # The global model is the sum of the Gaussian peaks

    if data is None:  # if we don't have data, the function only returns the direct calculation
        return [model, *peaks]
    if eps is None:  # without errors, no ponderation
        return (model - data)


def find_intencity_sf(structure_name, tol=1e-7, sigma=10, bounds=(900, 4000), max_iter=None, all_freqs=False):
    structure = load_structures([structure_name])[structure_name]
    if not all_freqs:
        calc_maxes = find_function_max(structure.pm.get_spectrum_function(sigma=sigma).calculate, structure.freqs, bounds=bounds)
        calc_maxes = calc_maxes[['peak_position', 'intencity']].drop_duplicates()
    else:
        calc_maxes = pd.DataFrame({'peak_position': structure.freqs, 'intencity': structure.intencities})

    intencity_norm_constant = np.max(calc_maxes.intencity)
    params = generate_parameters(intencity_norm_constant, calc_maxes.intencity.values, calc_maxes.peak_position.values)

    spectrum = np.genfromtxt('spectra_csv/absorbtivity/' + structure_name + '.CSV', delimiter=';')

    algo = 'nelder'
    x_fit = spectrum[:, 0]
    y_fit = spectrum[:, 1]

    # units_constant = gas_constant * 296 / 101.3e3 * 1e10 * 1e-5
    y_smo = rp.smooth(x_fit, y_fit, method="bartlett", window_length=31)
    norm_constant = np.max(y_smo) / 10
    y_smo /= norm_constant  # normalise spectra to maximum intensity, easier to handle
    result = lmfit.minimize(residual, params, method=algo, args=(x_fit, y_smo), tol=tol, bounds=bounds,
                            options={'maxiter': max_iter})  # fit data with  nelder model from scipy
    peaks = residual(result.params, x_fit)
    yout = peaks[0]

    psis = []
    sigmas = []
    mus = []
    units_constant = gas_constant * 296 / 101.3e3 * 1e10
    for i in range(len(calc_maxes.intencity.values)):
        i_i = result.params['i' + str(i)].value * units_constant * norm_constant * 1e-5
        sigma_i = result.params['sigma' + str(i)].value
        mu_i = result.params['mu' + str(i)].value

        psi = i_i * sigma_i * np.sqrt(2 * np.pi)
        psis.append(psi)
        sigmas.append(sigma_i)
        mus.append(mu_i)

    calc_maxes['psi'] = psis
    calc_maxes['sigma'] = sigmas
    calc_maxes['mu'] = mus
    calc_maxes['sf'] = calc_maxes.psi / calc_maxes.intencity

    fg_types = find_function_max(structure.pm.get_spectrum_function(sigma=sigma).calculate, structure.freqs, bounds=bounds)
    functional_groups = [' '.join(structure.get_fg_by_freq(freq)) for freq in fg_types.calc_freq]
    fg_types['functional_group'] = functional_groups

    plt.figure(figsize=(20, 10))
    plt.plot(x_fit, y_smo, 'k-', label='experimental spectrum')
    plt.plot(x_fit, yout, 'r-', label='gaussian approximation')
    plt.legend()
    plt.show()

    return calc_maxes, fg_types


if __name__ == '__main__':
    print(find_intencity_sf('isobutanol', tol=1, sigma=30, all_freqs=True))
