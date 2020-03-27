import matplotlib.pyplot as plt
import math
import numpy as np
from load_structures import load_structures
import utilities.functions as func
from constants import FREQ_RANGE, func_group_to_coords, fuc_group_colors, atoms_to_hb, atoms_to_epox

def make_bars(structure, column_names=None, color='#000000', width=7, label=None, coeff=1):
    if label is not None:
        label = label + ' ' + structure.name

    if column_names is None:
        x = [coeff * x for x in structure.freqs]
        y = structure.intencities
        plt.bar(x, y, width, color=color, label=label)

    else:
        for column_name in column_names:
            x = [coeff * x for x in structure.freqs]
            y = structure.intencities * structure.fg_coordinates_assigned_matrix[column_name]
            plt.bar(x, y, width, color=color, label=label)


def make_sum_spectrum(structures_, colorized_groups=None, full_spectrum=False):
    if colorized_groups is None:
        colorized_groups = []
    width = len(colorized_groups) + 5

    if full_spectrum:
        for structure in structures_:
            make_bars(structure, color='#d3d3d3', label='calc. spectrum')

    for group in colorized_groups:
        for structure in structures_:
            coords_names = func_group_to_coords[group]
            for k, coords_name in enumerate(coords_names):
                if k > 0:
                    label = None
                else:
                    label = group
                if coords_name in structure.fg_coordinates_assigned_matrix.columns:
                    make_bars(structure,
                              column_names=[coords_name],
                              color=fuc_group_colors[group],
                              label=label,
                              width=width)
                width -= 1


def read_gauss(filename='peaks.txt'):
    gauss_func_list = []
    with open(filename, 'r') as input_file:
        input_file.readline()
        lines = [line.strip('\n').split() for line in input_file.readlines()]
        for line in lines:
            if line[1] == 'Gaussian':
                mu = float(line[2])
                sigma = float(line[4]) / 2 / math.sqrt(2 * math.log(2, math.e))
                k = float(line[3]) * sigma * math.sqrt(2 * math.pi)
                gauss_func_list.append(k * func.GaussFunction(mu, sigma, FREQ_RANGE))

    return gauss_func_list


def make_func_group_gauss(structures=None, func_groups=None, sigma=30):
    if structures is None:
        structures = []
    if func_groups is None:
        func_groups = []

    for struct in structures:
        for func_group in func_groups:
            coords = func_group_to_coords[func_group]
            if coords[0] in struct.fg_coordinates_assigned_matrix.columns:
                for coord in coords:
                    mu = [x for x in np.array(struct.freqs) * struct.fg_coordinates_assigned_matrix[coord].to_numpy()
                          if x > 0]
                    intensity = [x for x in
                                 np.array(struct.intencities) * struct.fg_coordinates_assigned_matrix[coord].to_numpy()
                                 if x > 0]
                    for i in range(len(mu)):
                        gauss = intensity[i] * sigma * math.sqrt(2 * math.pi) * \
                                func.GaussFunction(mu_=mu[i], sigma_=sigma, x_=FREQ_RANGE)

                        color = fuc_group_colors[func_group]

                        plt.plot(gauss.x, gauss.y, color=color, label=func_group + ' ' + struct.name)


def make_func_group_gauss_sum(structures=None, func_groups=None, sigma=30, color=None, axis=plt):
    if structures is None:
        structures = []
    if func_groups is None:
        func_groups = []

    result = []

    for struct in structures:
        for func_group in func_groups:
            coords = func_group_to_coords[func_group]
            if coords[0] in struct.fg_coordinates_assigned_matrix.columns:
                sum_gauss = 0 * func.GaussFunction(mu_=1, sigma_=sigma, x_=FREQ_RANGE)
                added_freq = []
                for coord in coords:
                    check_list = [False] * len(struct.freqs)
                    for i, freq in enumerate(struct.freqs):
                        if freq not in added_freq:
                            check_list[i] = True
                    if coord in struct.fg_coordinates_assigned_matrix.columns:
                        mu = [x for x in np.array(struct.freqs) *
                              struct.fg_coordinates_assigned_matrix[coord].to_numpy() if x > 0]
                        intensity = [x for x in np.array(struct.intencities) *
                                     struct.fg_coordinates_assigned_matrix[coord].to_numpy() if x > 0]

                        for i in range(len(mu)):
                            added_freq.append(mu[i])

                            sum_gauss += sigma * math.sqrt(2 * math.pi) * intensity[i] * func.GaussFunction(
                                mu_=mu[i], sigma_=sigma, x_=FREQ_RANGE)

                if color is None:
                    new_color = fuc_group_colors[func_group]

                elif color == 'rand':
                    new_color = np.random.rand(3, )

                else:
                    new_color = color

                result.append((struct.name, func_group, sum_gauss))

                axis.plot(sum_gauss.x, sum_gauss.y, color=new_color, label=func_group + ' ' + struct.name)
                plt.xlabel('Frequency, cm-1')
                plt.ylabel('Intensity')
    return result


def get_group_peaks():
    result = dict()
    for fg in func_group_to_coords.keys():
        for struct in load_structures():
            if func_group_to_coords[fg][0] in struct.matrix_frame.columns:
                fg_curve = make_func_group_gauss_sum([fg], [struct])[0][1]
                prev_point = fg_curve.y[0]
                for i, point in enumerate(fg_curve.y[1:]):
                    if i != len(fg_curve.y[1:]) - 1:
                        next_point = fg_curve.y[i + 2]
                        if point > next_point and point > prev_point:
                            if fg not in result:
                                result[fg] = []
                            result[fg].append(fg_curve.x[i + 1])
    return result


def make_fg_type_bars(structure, group_name='hydroxide', type_='HB', axis=plt):
    if type_ == 'HB':
        new_farme = structure.fg_coordinates_assigned_matrix
        new_farme['freq'] = structure.freqs
        new_farme['intensity'] = structure.intencities
        new_farme['old_freq'] = structure.old_freqs

        new_farme['hb'] = [False] * new_farme.shape[0]
        new_farme['epox_hb'] = [False] * new_farme.shape[0]
        new_farme['type'] = ['STRE'] * new_farme.shape[0]
        for k, freq in enumerate(new_farme['old_freq']):
            for j in range(1):
                new_farme['type'][k] = structure.ped_manager.data[freq][j][1].type_
                for coord in structure.ped_manager.data[freq][j][1].coordinates:
                    for atom in coord.atoms:
                        if atom in atoms_to_hb[structure.name] and \
                                atom in atoms_to_hb[structure.name]:
                            new_farme['hb'][k] = True
                        if atom in atoms_to_epox[structure.name] and \
                                atom in atoms_to_epox[structure.name]:
                            new_farme['epox_hb'][k] = True
        new_farme['hydroxide'] = new_farme[[*(func_group_to_coords['hydroxide'])]].sum(axis=1)
        #
        hb_freqs = new_farme['freq'].loc[(new_farme['hydroxide'] > 0) & (new_farme['hb'])].values
        hb_int = new_farme['intensity'].loc[(new_farme['hydroxide'] > 0) & (new_farme['hb'])].values

        epox_hb_freqs = new_farme['freq'].loc[(new_farme['hydroxide'] > 0) & (new_farme['epox_hb'])].values
        epox_hb_int = new_farme['intensity'].loc[(new_farme['hydroxide'] > 0) & (new_farme['epox_hb'])].values

        non_hb_freqs = new_farme['freq'].loc[(new_farme['hydroxide'] > 0) & (new_farme['hb'] != True)
                                             & (new_farme['epox_hb'] != True)].values
        non_hb_int = new_farme['intensity'].loc[(new_farme['hydroxide'] > 0) & (new_farme['hb'] != True)
                                                & (new_farme['epox_hb'] != True)].values

        if len(hb_freqs) != 0:
            label = 'OH-OH hydrogen bonds'
            axis.bar(hb_freqs, hb_int, 5, color=fuc_group_colors['hydroxide'], label=label)

        if len(epox_hb_freqs) != 0:
            label = 'OH-epoxide hydrogen bonds'
            axis.bar(epox_hb_freqs, epox_hb_int, 5, color=fuc_group_colors['epoxide'], label=label)

        label = 'without hydrogen bonds'
        axis.bar(non_hb_freqs, non_hb_int, 5, color='r', label=label)