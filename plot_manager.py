import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utilities.functions as func
from constants import *
from find_function_max import find_function_max


class PlotManager:
    def __init__(self, structure):
        self.structure = structure

    def getattr(self, attr_name):
        return getattr(self, attr_name)

    def spectrum_bar(self, axis=plt, color=None):
        if color is None:
            color = '#808080'
        elif color == 'rand':
            color = np.random.rand(3, )

        x = self.structure.freqs
        y = self.structure.intencities
        label = 'calc. spectrum ' + self.structure.name
        axis.bar(x, y, 10, color=color, label=label)

    def get_spectrum_function(self, fg=None, sigma=30, curve_type='gauss'):
        if fg is None:
            freqs = self.structure.freqs
            intencities = self.structure.intencities
        else:
            freqs = (pd.DataFrame(self.structure.freqs).loc[self.structure.fg_assigned_matrix[fg] == 1]).T.values[0]
            intencities = \
                (pd.DataFrame(self.structure.intencities).loc[self.structure.fg_assigned_matrix[fg] == 1]).T.values[0]

        x = FREQ_RANGE
        y = 0 * func.GaussFunction(mu_=0, sigma_=sigma, x_=x)

        for i, freq in enumerate(freqs):
            if curve_type == 'gauss':
                y += intencities[i] * sigma * \
                     math.sqrt(2 * math.pi) * func.GaussFunction(mu_=freq, sigma_=sigma, x_=x)
            elif curve_type == 'lorenz':
                y += intencities['intensity'][i] * func.LorenzFunction(x0_=freq, gamma_=sigma, x_=x)
        return y

    def annotate(self, function, points, axis=plt):
        annotation_points = find_function_max(function, points, bounds=(900, 4000)).drop(['calc_freq'], axis=1).drop_duplicates()
        for point in annotation_points.iterrows():
            text = str(int(point[1].peak_position))
            axis.annotate(text, (point[1].peak_position, point[1].intencity), clip_on=True)

    def spectrum_curve(self, axis=plt, color=None, sigma=30, curve_type='gauss', annotate=True):
        if color is None:
            color = '#808080'
        elif color == 'rand':
            color = np.random.rand(3, )

        if curve_type == 'gauss':
            param_name = 'sigma'
        else:
            param_name = 'gamma'

        y = self.get_spectrum_function(sigma=sigma, curve_type=curve_type)

        label = 'spectrum curve ' + self.structure.name + ', ' + param_name + '=' + str(sigma)
        axis.plot(y.x, y.y, label=label, color=color)
        if annotate:
            self.annotate(y.calculate, self.structure.freqs, axis)
        return y

    def fg_bar(self, fg_list, axis=plt, color_list=None):
        if not isinstance(fg_list, list):
            fg_list = [fg_list]

        for i, fg in enumerate(fg_list):
            if fg in self.structure.fg_assigned_matrix.columns:
                if color_list is None:
                    color = fuc_group_colors[fg]
                elif color_list == 'rand':
                    color = np.random.rand(3, )
                else:
                    color = color_list[i]

                x = (pd.DataFrame(self.structure.freqs).loc[self.structure.fg_assigned_matrix[fg] == 1]).T.values[0]
                y = (pd.DataFrame(self.structure.intencities).loc[self.structure.fg_assigned_matrix[fg] == 1]).T.values[0]

                label = fg + ' ' + self.structure.name
                axis.bar(x, y, width=10, color=color, label=label)

    def fg_sum_curve(self, fg_list, axis=plt, color_list=None, sigma=30, curve_type='gauss', annotate=True):
        if not isinstance(fg_list, list):
            fg_list = [fg_list]

        curves = []
        for i, fg in enumerate(fg_list):
            if fg in self.structure.fg_assigned_matrix.columns:
                if color_list is None:
                    color = fuc_group_colors[fg]
                elif color_list == 'rand':
                    color = np.random.rand(3, )
                else:
                    color = color_list[i]

                if curve_type == 'gauss':
                    param_name = 'sigma'
                else:
                    param_name = 'gamma'

                y = self.get_spectrum_function(fg=fg, sigma=sigma, curve_type=curve_type)

                label = fg + ' ' + self.structure.name + ', ' + param_name + '=' + str(sigma)
                axis.plot(y.x, y.y, label=label, color=color)
                curves.append(y)
                if annotate:
                    self.annotate(y.calculate, self.structure.get_fg_freqs(fg, corrected=True), axis)
        return curves


