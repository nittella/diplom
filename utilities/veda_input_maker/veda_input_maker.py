class GaussianCoord:
    def __init__(self, x_, y_, z_):
        self.x = x_
        self.y = y_
        self.z = z_

    def __str__(self):
        return '   '.join(map(
            str,
            [round(self.x, 6),
             round(self.y, 6),
             round(self.z, 6)]
        )
        )


def read_force_matrix(file_name):
    matrix_dict = dict()
    with open(file_name) as input_file:

        for line in input_file.readlines():
            line.strip('\n')
            list_line = line.split()
            if len(list_line) > 1 and '.' in list_line[1]:
                if list_line[0] not in matrix_dict:
                    matrix_dict[list_line[0]] = []
                matrix_dict[list_line[0]] += [float(x) for x in list_line[1:]]
    return matrix_dict


def read_cartesian_coordinates(file_name):
    with open(file_name) as input_file:
        atomic_data = dict()
        atomic_data['center_number'] = []
        atomic_data['atomic_number'] = []
        atomic_data['atomic_type'] = []
        atomic_data['coordinates'] = []
        center_number = 1
        for line in input_file.readlines():
            line.strip('\n')
            list_line = line.split()
            atomic_number = int(list_line[0])
            atomic_data['center_number'].append(center_number)
            atomic_data['atomic_number'].append(atomic_number)
            atomic_data['atomic_type'].append(0)
            x = float(list_line[1])
            y = float(list_line[2])
            z = float(list_line[3])
            atomic_data['coordinates'].append([x, y, z])
            center_number += 1
    return atomic_data


def generate_dmatrix_size(n):
    sum_ = 0
    for i in range(n * 3 + 1):
        sum_ += i
    return sum_


def write_in_columns(list_, file_, columns=5):
    counter = 1
    for line in list_:
        for elem in line:
            if isinstance(elem, int):
                file_.write(' ' * 11)
                file_.write(str(elem))
            elif isinstance(elem, float):
                if elem >= 0:
                    file_.write(' ' * 2)
                else:
                    file_.write(' ')
                file_.write(str("%10.8E"% elem))
            if counter == columns:
                counter = 1
                file_.write('\n')
            else:
                counter += 1
    file_.write('\n')


def create_veda_output(force_matrix_, atomic_data_, file_name_='output.fmu'):
    with open(file_name_, 'w') as output_file:
        output_file.write(
            'Atomic numbers                             I   N=          '
            + str(len(atomic_data_['atomic_number']))
            + '\n')
        write_in_columns([atomic_data_['atomic_number']], output_file)
        output_file.write(
            'Current cartesian coordinates              R   N=          '
            + str(len(atomic_data_['coordinates']) * 3)
            + '\n')
        write_in_columns(atomic_data_['coordinates'], output_file)
        output_file.write(
            'Int Atom Types                             I   N=          '
            + str(len(atomic_data_['atomic_type']))
            + '\n')
        write_in_columns([atomic_data_['atomic_type']], output_file)
        output_file.write(
            'Cartesian Force Constants                  R   N=         '
            + str(generate_dmatrix_size(len(atomic_data_['atomic_type'])))
            + '\n')
        write_in_columns(force_matrix_.values(), output_file)
    print('Output have been successfully created!')


if __name__ == '__main__':
    coords = read_cartesian_coordinates('coords')
    force_matrix = read_force_matrix('matrix.mol')
    create_veda_output(force_matrix, coords)
