import dataclasses
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cosine

@dataclasses.dataclass
class Coordinate:
    coefficient = 0
    atoms = [1, 2]
    name = 'CC'
    type = 'STRE'
    value = 1.0

    def __str__(self):
        return ' '.join([
            str(self.coefficient),
            self.name,
            self.type,
            str(self.atoms)
        ])


class ComplexCoordinate:
    def __init__(self):
        self.type = 'STRE'
        self.name = 'CC'
        self.coordinates = [Coordinate()]
        self.number = 0
        self.o_coordinate = False
        self.o_share = 0
        self.is_epoxyde = False
        self.is_hydroxyde = False

    def __str__(self):
        result = ' '.join(['№' + str(self.number),
                           self.type,
                           self.name,
                           'O-share: ' + str(round(self.o_share, 2))])
        return result + '\n' + '\n'.join(map(str, self.coordinates)) + '\n' + '*' * 50


class PEDManager:
    def __init__(self, matrix_, freqs_, coordinates_, sorting_type='o'):
        self.matrix = matrix_
        self.coordinates = coordinates_
        self.data = dict() # key is frequency keep sorting by value * o_share pair value + coordinate
        self.data2 = dict() # key is coordunate number, keep sorted by value pair value + freq
        self.sorted_by_value = dict() # key is frequency keep sorting by value pair value + coordinate
        for key_freq, matrix_line in zip(freqs_, matrix_):
            self.data[key_freq] = [x for x in zip(matrix_line, coordinates_)]
            self.sorted_by_value[key_freq] = [x for x in zip(matrix_line, coordinates_)]
            if sorting_type == 'o':
                self.data[key_freq].sort \
                    (key=lambda x: (-abs(x[0] * x[1].o_share), -abs(x[0])))
            else:
                self.data[key_freq].sort(
                    key=lambda x: -abs(x[0]))
            self.sorted_by_value[key_freq].sort(key=lambda x: -abs(x[0]))
            for coordinate, elem in zip(coordinates_, matrix_line):
                if coordinate.number not in self.data2:
                    self.data2[coordinate.number] = []
                self.data2[coordinate.number].append((elem, key_freq))
        for coord_number in self.data2:
            self.data2[coord_number].sort(key=lambda x: -abs(x[0]))

    def __str__(self):
        result = ''
        for freq, line in self.data.items():
            if line[0][1].o_coordinate:
                result += '*'
            result += str(freq) + ': '
            for value, coordinate in line:
                if abs(value) > -1:
                    result += '(' + str(value) + ' ' \
                              + str(coordinate) + ')' + ', '
            result += '\n'
        return result

    def choose_freq(self, max_o_value=8.):
        # print(self.get_first_o())
        result = [y for x, y in zip(self.get_first_and_half(), self.get_freqs()) if x >= max_o_value]
        return result

    def print_data(self, o_marker=False, only_o=False,
                   abs_value_greater_than=-1):
        for freq, line in self.data.items():
            if line[0][1].o_coordinate and o_marker:
                print('*', end='')
            print(freq, end=': ')
            for value, coordinate in line:
                if abs(value) > abs_value_greater_than:
                    print('(' + str(value), str(coordinate) + ')', end=', ')
            print('\n', end='')

    def get_freq_by_coord(self, coordinate_number):
        return self.data2[coordinate_number]

    def get_coord_by_freq(self, freq, sort=None):
        if sort == 'by_value':
            return self.sorted_by_value[freq]
        return self.data[freq]

    def get_max_values(self):
        max_values = []
        for freq, line in self.data.items():
            max_values.append(abs(line[0][0]))
        return max_values

    def get_freqs(self):
        return [i for i in self.data.keys()]

    def get_o_sum(self, n=-1):
        o_sums = []
        for i, (freq, line) in enumerate(self.data.items()):
            o_sum = 0
            for value, coordinate in line:
                if coordinate.o_coordinate:
                    o_sum += abs(value) * coordinate.o_share
            o_sums.append(o_sum)
        return o_sums

    def make_o_sum_plot(self, freq):
        o_sums = [0]
        for value, coordinate in self.data[freq]:
            new_sum = o_sums[-1] + abs(value) * coordinate.o_share
            o_sums.append(new_sum)
        plt.plot(o_sums, label=str(freq))
        return o_sums

    def get_mean_for_freq(self, only_non_zero=False):
        abs_matrix = abs(self.matrix)
        if only_non_zero:
            for x in abs_matrix:
                print(len(*x.nonzero()))
            print([len(*x.nonzero()) for x in abs_matrix])
            return abs_matrix.sum(axis=1) / [len(*x.nonzero()) for x in
                                             abs_matrix]
        return abs_matrix.mean(axis=1)

    def get_first_o(self, n=0):
        result = []
        for freq, line in self.data.items():
            sum_n = 0
            for i in range(n + 1):
                sum_n += abs(line[i][0] * line[i][1].o_share)
            result.append(sum_n)
        return result

    def get_first_and_half(self):
        result = []
        for freq, line in self.data.items():
            sum_n = abs(line[0][0] * line[0][1].o_share)
            first = abs(line[0][0] * line[0][1].o_share)
            i = 1
            while i < len(line) and (
                    first / 2 < abs(line[i][0] * line[i][1].o_share) or \
                    abs(line[i][0] * line[i][1].o_share) == first / 2):
                sum_n += abs(line[i][0] * line[i][1].o_share)
                i += 1
            result.append(sum_n)
        return result

    def get_o_coordinates(self):
        result = []
        for coordinate in self.coordinates:
            if coordinate.o_coordinate:
                result.append(coordinate)
                print(coordinate)
        return result

    def make_plot(self, x1, x2, kind='scatter', n=0, normalize=False):
        axis = [None, None]
        labels = [None, None]
        for i, x_i in enumerate((x1, x2)):
            labels[i] = x_i
            if x_i == 'freq':
                axis[i] = self.get_freqs()
            elif x_i == 'o_sum':
                axis[i] = self.get_o_sum()
            elif x_i == 'max_value':
                axis[i] = self.get_max_values()
            elif x_i == 'av_freq':
                axis[i] = self.get_mean_for_freq()
            elif x_i == 'av_non_zero_freq':
                axis[i] = self.get_mean_for_freq(only_non_zero=True)
            elif x_i == 'first_o':
                axis[i] = self.get_first_o(n)
                # if normalize:
                #     o_sum = np.array(self.get_o_sum())
                #     o_sum[o_sum == 0] = 1
                #     axis[i] = np.array(self.get_first_o()) / np.array(o_sum)
            elif x_i == 'first_o_and_half':
                axis[i] = self.get_first_and_half()
        title = labels[1] + '(' + labels[0] + ')'
        if kind == 'scatter':
            plt.scatter(axis[0], axis[1], label=title)
        elif kind == 'plot':
            plt.plot(axis[0], axis[1], label=title)
        elif kind == 'hist':
            plt.hist(axis[0], axis[1], bins=10)
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])

    def cut_freq_range(self, freq_start, freq_stop):
        freq_for_del = [freq for freq in self.data.keys()
                        if freq_start > freq or freq > freq_stop]
        for freq in freq_for_del:
            del self.data[freq]
            del self.sorted_by_value[freq]


def read_freq(file_name_='skra.ved'):
    freq = []
    with open(file_name_, encoding='utf8') as input_file:
        line = input_file.readline().strip('\n').split()
        while len(line) == 0 or line[0] != 'PED:':
            if len(line) > 0 and '.' in line[0]:
                for f in line:
                    freq.append(float(f))
            line = input_file.readline().strip('\n').split()
    return freq


def read_veda_matrix(file_name_='skra.ved'):
    with open(file_name_, encoding='utf8') as input_file:
        line = input_file.readline().strip('\n').split()
        while len(line) == 0 or line[0] != 'PED:':
            line = input_file.readline().strip('\n').split()
        line = input_file.readline().strip('\n').split()
        matrix_size = len(line)
        result = np.ndarray((matrix_size, matrix_size))
        for i in range(matrix_size):
            line = input_file.readline().strip('\n').split()
            for j, elem in enumerate(line[1:-1]):
                result[i, j] = elem
    return result


def read_coord(n_, file_name_='skra.dd2'):
    result = []
    with open(file_name_, encoding='utf8') as input_file:
        input_file.readline()
        input_file.readline()
        line = input_file.readline().strip('\n').split()
        for i in range(n_):
            o_coordinates_count = 0
            complex_coordinate = ComplexCoordinate()
            complex_coordinate.coordinates = []
            coordinate = Coordinate()
            complex_coordinate.number = int(line[1])
            coordinate.coefficient = float(line[2])
            complex_coordinate.type = line[3]
            coordinate.type = line[3]
            j = 4
            coordinate.atoms = []
            if not line[j].isdigit():
                j += 1
            while line[j].isdigit():
                if int(line[j]) in [72, 73, 74]:
                    complex_coordinate.is_epoxyde = True
                elif int(line[j]) in [76, 77, 78, 79, 80, 81]:
                    complex_coordinate.is_hydroxyde = True
                coordinate.atoms.append(int(line[j]))
                j += 1
            coordinate.name = line[j]
            if 'O' in coordinate.name:
                complex_coordinate.o_coordinate = True
                o_coordinates_count += 1
            complex_coordinate.name = line[j]
            coordinate.value = float(line[j + 1])
            complex_coordinate.coordinates.append(coordinate)
            line = input_file.readline().strip('\n').split()
            while line[0] != 's' and line != ['****']:
                coordinate = Coordinate()
                coordinate.coefficient = float(line[0])
                j = 1
                coordinate.atoms = []
                coordinate.type = complex_coordinate.type
                if not line[j].isdigit():
                    coordinate.type = line[j]
                    j += 1
                while line[j].isdigit():
                    if int(line[j]) in [72, 73, 74]:
                        complex_coordinate.is_epoxyde = True
                    elif int(line[j]) in [76, 77, 78, 79, 80, 81]:
                        complex_coordinate.is_hydroxyde = True
                    coordinate.atoms.append(int(line[j]))
                    j += 1
                coordinate.name = line[j]
                if 'O' in coordinate.name:
                    complex_coordinate.o_coordinate = True
                    o_coordinates_count += 1
                coordinate.value = float(line[j + 1])
                complex_coordinate.coordinates.append(coordinate)
                line = input_file.readline().strip('\n').split()
            complex_coordinate.o_share = o_coordinates_count / len(
                complex_coordinate.coordinates)
            result.append(complex_coordinate)
    return result


def assign_freq(internal_ped_manager, external_ped_manager):
    assigned_freqs = []
    deltas = []
    not_assigned = set(external_ped_manager.get_freqs())
    for i, internal_freq in enumerate(internal_ped_manager.get_freqs()):
        deltas.append(4000)
        assigned_freqs.append([internal_freq, external_ped_manager.get_freqs()[0]])
        for j, external_freq in enumerate(external_ped_manager.get_freqs()):
            if abs(external_freq - internal_freq) < deltas[i]:
                assigned_freqs[i][1] = external_freq
                deltas[i] = external_freq - internal_freq
    for freqs in assigned_freqs:
        if freqs[1] in not_assigned:
            not_assigned.remove(freqs[1])
    return assigned_freqs, deltas, sorted([float(x) for x in not_assigned])
    # print('max delta:', max(deltas))
    # print(*zip(assigned_freqs, coords), sep='\n')


def analyze_mode_type(coordinates_1):
    all_types = {'TORS', 'BEND', 'OUT', 'STRE', 'TORS_O', 'BEND_O', 'OUT_O', 'STRE_O'}
    vector_1_ = dict()
    for coord_type in all_types:
        vector_1_[coord_type] = 0
    for complex_coordinate in coordinates_1:
        for coordinate in complex_coordinate[1].coordinates:
            coord_type = coordinate.type
            if 'O' in coordinate.name:
                coord_type += '_O'
            vector_1_[coord_type] += abs(complex_coordinate[0]) / len(complex_coordinate[1].coordinates)
    # print([x for x in vector_1.values()])
    # for key in ['TORS', 'BEND', 'OUT', 'STRE', 'TORS_O', 'BEND_O', 'OUT_O', 'STRE_O']:
    #     vector_1_.pop(key)
    return [x for x in vector_1_.keys()], [x for x in vector_1_.values()]


if __name__ == '__main__':
    ped_dict = dict()
    for structure_name, coords_number in [('104', 219),
                                          ('104oh6', 255),
                                          ('g104w3', 246)]:
        source_name = 'structures/' + structure_name + '/skra.'
        ped_dict[structure_name] = PEDManager(read_veda_matrix(
            file_name_=source_name + 'ved'),
            read_freq(file_name_=source_name + 'ved'),
            read_coord(coords_number, file_name_=source_name + 'dd2')
        )
        ped_dict[structure_name].cut_freq_range(900, 1800)

    freqs, deltas, _ = assign_freq(ped_dict['104'], ped_dict['g104w3'])


    labels = []
    vectors = []
    for freq in freqs:
        labels, vector = analyze_mode_type(ped_dict['104'].get_coord_by_freq(freq[0]))
        vectors.append(vector)
    data = np.array(vectors)

    labels_2 = []
    vectors_2 = []
    for freq in freqs:
        labels_2, vector_2 = analyze_mode_type(ped_dict['g104w3'].get_coord_by_freq(freq[1]))
        vectors_2.append(vector_2)
    data_2 = np.array(vectors_2)

    distances = []
    for i, vector_1 in enumerate(data):
        distances.append(cosine(vector_1, data_2[i]))
    plt.hist(distances)
    plt.show()

    for freq in freqs:
        labels_2, vector_2 = analyze_mode_type(ped_dict['g104w3'].get_coord_by_freq(freq[1]))
        vectors_2.append(vector_2)
    data_2 = np.array(vectors_2)
