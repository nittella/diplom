from scipy.optimize import minimize
import numpy as np
import pandas as pd


# Ищет ближайшие максимумы функуции function для точек points в диапазоне bounds
def find_function_max(function, points=None, bounds=(0, 1)):
    if points is None:
        points = np.array([0.5])

    points = points[(points > bounds[0]) & (points < bounds[1])]
    result = pd.DataFrame(columns=['calc_freq', 'peak_position', 'intencity'])
    for point in points:
        res = minimize(lambda x: -function(x[0]), np.array([int(point)]), method='CG')
        peak_position = res.x[0]
        peak_value = -res.fun
        result = result.append({'calc_freq': point,
                                'peak_position': round(peak_position, 0),
                                'intencity': round(peak_value, 0)},
                               ignore_index=True)

    return result
