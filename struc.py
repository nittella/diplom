import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from constants.func_groups import *
from utilities import functions as func
from constants.constants import *
from internal_coordinate import ComplexCoordinate, Coordinate
from plot_manager import PlotManager
from load_structures import load_structures
from generate_coord_set import generate_coord_set
from atom import read_atoms

pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_rows', None)


# scaling factor определен при помоще статитики сравнения расчитанных небольших молекул методом B3LYP с базисом 6-31+G**
# погрешность для диапазона (1000, 2000) +-0.0344, для диапазона (2900, 4000) +- 0.0298
def sf(freq):
    if freq < 2600:
        return -8e-7 + 0.9788
    return 0.9574


def get_number_of_atoms(name):
    path = 'structures/' + name + '/coords'
    with open(path, 'r') as input_file:
        return len(input_file.readlines())


class Structure:
    def __init__(self, name, ped_manager, sf_function=sf):
        self.name = name
        self.path = 'structures/' + self.name + '/'  # путь к папке с файлами структуры

        self.old_freqs = np.array(self.get_freqs(sf_function=lambda x: 1))  # исходные частоты
        self.freqs = np.array(self.get_freqs(sf_function=sf_function))  # откорректированные частоты

        self.atoms = read_atoms(path=self.path)
        self.number_of_atoms = get_number_of_atoms(self.name)

        self.veda_matrix = self.get_veda_matrix()   # матрица разложения потенциальной энергии по внутренним координатам
        # внутренние координаты, по которым проведено разложение внутренней энергии:
        self.internal_coordinates = self.get_internal_coordinates()
        #  разложение потенциальной энергии в координатах функциональных групп
        self.fg_coordinates_matrix = self.generate_fg_coordinates_matrix()
        # Матрица из 0 и 1. 1 там где координата вносит значительный вклад в потенциальну энергию
        self.fg_coordinates_assigned_matrix = self.choose_freqs_by_group()
        # Матрица из 0 и 1. 1 там где функциональная группа вносит значительный вклад в потенциальну энергию
        self.fg_assigned_matrix = self.assign_coords_to_fg()

        self.ped_manager = ped_manager  # хранит и сортирует координаты
        self.freq_range = (0, 4000)  # ограничение по частотам
        self.intencities = self.read_intencities()  # расчетные интенсивности
        self.pm = PlotManager(self)  # строит графики

# Считывает частоты из файла и масштабирует их согласно sf_function
    def get_freqs(self, sf_function):
        if sf_function is None:
            sf_function = sf

        freqs = []
        freq_path = self.path + 'skra.ved'
        with open(freq_path, encoding='utf8') as input_file:
            line = input_file.readline().strip('\n').split()
            while len(line) == 0 or line[0] != 'PED:':
                if len(line) > 0 and '.' in line[0]:
                    for f in line:
                        f = sf_function(float(f)) * float(f)
                        freqs.append(f)
                line = input_file.readline().strip('\n').split()
        return freqs

# Считывает матрицу разложения потенциальной энергии по внутренним координатам из файла
    def get_veda_matrix(self):
        matrix_path = self.path + 'skra.ved'
        with open(matrix_path) as input_file:
            lines = [x.strip('\n').split() for x in input_file.readlines()]
        start = 0
        for start, line in enumerate(lines):
            if line != [] and line[0] == 'PED:':
                start += 2
                break
        end = 0
        for end, line in enumerate(lines):
            if end > start and line == []:
                end -= 1
                break

        size = len(lines[start]) - 2
        matrix = np.zeros((size, size))
        for freq_id, line in enumerate(lines[start:end]):
            for coord_id, elem in enumerate(line[1:-1]):
                matrix[freq_id, coord_id] = int(elem)
        table = pd.DataFrame(matrix)
        table['freq'] = self.old_freqs
        return table

# Считывает внутренние координаты из файла
    def get_internal_coordinates(self):
        coord_path = self.path + 'skra.dd2'
        with open(coord_path) as input_file:
            lines = [x.strip('\n').split() for x in input_file.readlines()]
        start = 2
        end = 0
        for end, line in enumerate(lines):
            if line != [] and line[0] == '****':
                break

        coordinates = []
        for line in lines[start:end]:
            if line[0] in ('s', 'k', 'v'):
                # complex coordinate parameters
                coord_type = line[3]
                number = int(line[1])

                # coordinate parameters
                coefficient = line[2]
                if coord_type == 'STRE':
                    atoms = [int(x) for x in line[4:6]]
                    name = line[6]
                    value = float(line[7])
                elif coord_type == 'BEND':
                    atoms = [int(x) for x in line[4:7]]
                    name = line[7]
                    value = float(line[8])
                else:
                    atoms = [int(x) for x in line[4:8]]
                    name = line[8]
                    value = float(line[9])
                coordinate = Coordinate(coefficient, atoms, name, value)
                coordinates.append(ComplexCoordinate(coord_type=coord_type, number=number, coordinates=[coordinate]))
            else:
                coord_type = coordinates[-1].type
                # coordinate parameters
                coefficient = line[0]
                i = 1
                if line[i] in ('STRE', 'BEND', 'TORS', 'OUT'):
                    i = 2
                if coord_type == 'STRE':
                    atoms = [int(x) for x in line[i:i + 2]]
                    name = line[i + 2]
                    value = float(line[i + 3])
                elif coord_type == 'BEND':
                    atoms = [int(x) for x in line[i:i + 3]]
                    name = line[i + 3]
                    value = float(line[i + 4])
                else:
                    atoms = [int(x) for x in line[i:i + 4]]
                    name = line[i + 4]
                    value = float(line[i + 5])
                coordinate = Coordinate(coefficient, atoms, name, value)
                coordinates[-1].add_coordinate(coordinate)
        return coordinates

# Считывает расчетные интенсивности из файла
    def read_intencities(self, file_name=None):
        intencities = []
        if file_name is None:
            file_name = 'structures/' + self.name + '/intencities'
        with open(file_name, 'r') as input_file:
            lines = [line.strip('\n').split() for line in input_file.readlines()]
            lines.sort(key=lambda x: -float(x[0]))
            for line in lines[:]:
                if self.freq_range[0] < float(line[0]) < self.freq_range[1]:
                    intencities.append(float(line[1]))
        return intencities

# По номеру координаты возвращает имя соответвующей координаты функциональной группы
    def get_fg_coord_name(self, complex_coord_number):
        name = ''

        for coordinate in self.internal_coordinates[complex_coord_number].coordinates:
            fg_list = []
            for atom in coordinate.atoms:
                if self.atoms[atom].functional_group not in fg_list:
                    fg_list.append(self.atoms[atom].functional_group)
            fg_list.sort()
            name = '_'.join(fg_list)
        if self.internal_coordinates[complex_coord_number].type == 'STRE':
            name += '_STRE'
        else:
            name += '_ANG'
        return name

# Матрица частоты - координаты функциональных групп
    def generate_fg_coordinates_matrix(self):
        coords = list(generate_coord_set(func_groups={atom.functional_group for atom in self.atoms}))
        matrix = np.zeros((len(self.freqs), len(coords)))
        for i, freq_line in self.veda_matrix.iterrows():
            for ind in freq_line[:-1].to_numpy().nonzero()[0]:
                name = self.get_fg_coord_name(ind)
                j = coords.index(name)
                matrix[i, j] += abs(freq_line[ind])
        frame = pd.DataFrame(matrix)
        frame.columns = coords
        frame = frame.loc[:, (frame != 0).any(axis=0)]
        return frame

# Относит частоты к той или иной координате по некому правилу
    def choose_freqs_by_group(self):
        result = pd.DataFrame()
        abs_fg_coordinates_matrix = self.fg_coordinates_matrix.abs()
        for column in abs_fg_coordinates_matrix.columns[:]:
            if max(abs_fg_coordinates_matrix[column]) > 33:
                threshold = 0.1 * max(abs_fg_coordinates_matrix[column])
            else:
                q_25 = abs_fg_coordinates_matrix[column].loc[abs_fg_coordinates_matrix[column] > 0].quantile(q=0.25)
                q_75 = abs_fg_coordinates_matrix[column].loc[abs_fg_coordinates_matrix[column] > 0].quantile(q=0.75)
                threshold = q_75 + 1.5 * (q_75 - q_25)
            result[column + '_assignment'] = np.where(abs_fg_coordinates_matrix[column] > threshold, 1, 0)
        return result

# Относит координаты к функциональной группе
    def assign_coords_to_fg(self):
        result = pd.DataFrame(data=np.zeros((len(self.freqs), len(func_group_to_coords.keys()))),
                              columns=func_group_to_coords.keys())
        for i, freq in enumerate(self.freqs):
            for fg, list_coords in func_group_to_coords.items():
                for coord in list_coords:
                    if coord in self.fg_coordinates_assigned_matrix.columns and \
                            self.fg_coordinates_assigned_matrix[coord][i] == 1:
                        result[fg][i] = 1
                        break
        result = result.loc[:, (result != 0).any(axis=0)]
        return result

    # Возвращает частоты относящиеся только к указанной функциональной группе. Если corrected=True, то частоты
    # откорректированные
    def get_fg_freqs(self, fg, corrected=False):
        if corrected:
            return self.freqs[self.fg_assigned_matrix[fg] == 1]
        else:
            return self.old_freqs[self.fg_assigned_matrix[fg] == 1]

    # По указанной частоте freq возвращает индексы координат, отсортированные по убыванию вкладов в потенциальную энергию, так чтобы
    # сумма их вкладов была не меньше чем min_sum или минимальное значение больше threshold.
    def _get_coordinates_by_freq(self, freq, min_sum=80, threshold=3):
        sum_ = 0
        indexes = []
        values = []
        coordinates = []
        freq_series_copy = self.veda_matrix.loc[self.veda_matrix.freq == freq].drop(['freq'], axis=1).abs()[:]
        value = threshold
        while sum_ < min_sum and value >= threshold:
            coord_ind = freq_series_copy.idxmax(axis=1).values[0]
            value = freq_series_copy[coord_ind].values[0]
            coord = self.internal_coordinates[coord_ind]

            sum_ += value

            if value > threshold:
                indexes.append(coord_ind)
                values.append(value)
                coordinates.append(coord)

            freq_series_copy[coord_ind].values[0] = 0

        return indexes, values, coordinates

    # Возвращает функциональные группы, вносящие основной вклад в колебание с данной частотой
    def get_fg_by_freq(self, freq, corrected=True):
        freqs = self.freqs
        if not corrected:
            freqs = self.old_freqs

        freq_line = (self.fg_assigned_matrix[abs(freqs - freq) < 1e-4] == 1).values[0]
        return self.fg_assigned_matrix.columns[freq_line].values

# Возвращает таблицу с частотами, характерными для указанной функиональной группы в диапазне freq_range, вклады
# внутренних координат в потенциоальную энергию и сами внутренние координаты
    def get_report(self, fg, freq_range=(900, 4000), corrected=False):
        report_table = pd.DataFrame(columns=['freq, cm-1', 'values, %', 'internal_coordinates'])

        fg_freqs = self.get_fg_freqs(fg, corrected=False)
        fg_freqs_in_range = fg_freqs[(fg_freqs < freq_range[1]) & (fg_freqs > freq_range[0])]

        corrected_fg_freqs = self.get_fg_freqs(fg, corrected=True)
        corrected_fg_freqs_in_range = corrected_fg_freqs[(fg_freqs < freq_range[1]) & (fg_freqs > freq_range[0])]

        for i, freq in enumerate(fg_freqs_in_range):
            indexes, values, coordinates = self._get_coordinates_by_freq(freq, min_sum=100)
            values_str = '+'.join(map(lambda x: str(int(x)), values))
            coords_str = '+'.join([x.coordinates[0].name + '(' + '-'.join(map(str, x.coordinates[0].atoms)) + ')'
                                   for x in coordinates])
            if corrected:
                freq_for_printing = corrected_fg_freqs_in_range[i]
            else:
                freq_for_printing = freq

            report_table = report_table.append({'freq, cm-1': freq_for_printing,
                                                'values, %': values_str,
                                                'internal_coordinates': coords_str},
                                               ignore_index=True)
        return report_table