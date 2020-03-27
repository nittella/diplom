import numpy as np

with open('force_matrix.mol') as input_file:
    matrix_dict = dict()
    for line in input_file.readlines():
        line.strip('\n')
        list_line = line.split()
        if '.' in list_line[1]:
            if list_line[0] not in matrix_dict:
                matrix_dict[list_line[0]] = []
            matrix_dict[list_line[0]] += [float(x) for x in list_line[1:]]
result = ''
counter = 0
for i, (key, value) in enumerate(matrix_dict.items()):
    for j, elem in enumerate(value):
        if elem >= 0:
            result += '  '
        else:
            result += ' '
        result += str("%10.8E"% (elem))
        if counter == 4:
            counter = 0
            result += '\n'
        else:
            counter += 1
with open('gauss_force_matrix', 'w') as output_file:
    output_file.write(result)

