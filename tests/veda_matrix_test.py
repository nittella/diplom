from load_structures import load_structures
from create_structures import create_structures
import numpy as np


def veda_matrix_test():
    # test on water
    create_structures(['water1'])
    struct_water = load_structures(['water1'])['water1']
    water_control_matrix = np.array([[100, 0, 0],
                                     [0, 100, 0],
                                     [0, 0, 100]])
    read_matrix = struct_water.veda_matrix.iloc[:, :-1].to_numpy()
    if np.array_equal(read_matrix, water_control_matrix):
        print('Water test OK!')
    else:
        print('Water test error')
        print(read_matrix)

#     water trimer test
    create_structures(['water3'])
    struct_water3 = load_structures(['water3'])['water3']
    water_control_matrix = np.array([[34, 28, 0, 12, 0, -1, 0, 0, 0, 0, 0, 0, 22, 0, 0, 3, 0, 0, 0, 0, 0],
                                     [0, -26, 1, 71, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [-23, 41, -4, 11, -2, 1, 0, 0, 0, 0, 0, 0, -15, 0, 0, -2, 0, 0, 0, 0, 0],
                                     [0, -1, -25, 3, 70, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [2, -1, -26, 0, -3, 64, 0, 0, -2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 2, 43, 1, 21, 31, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [5, 0, 0, 0, 0, 0, 1, 0, -6, 1, 10, 1, -26, 1, 9, -1, 16, 4, -11, -4, 4],
                                     [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -12, 1, 1, -1, 19, 0, -15, 0, 8, 39, -1],
                                     [11, 0, 0, 0, 0, 0, 0, 0, -5, 0, -5, -2, -25, 0, -8, -2, 9, 11, 15, 7, 0],
                                     [-1, 0, 0, 0, 0, 0, 1, -1, -16, 1, 18, -4, 0, -4, 23, 1, -9, -4, 7, 3, -7],
                                     [0, 0, 0, 0, 0, 0, -1, 0, 12, -4, 17, -9, 1, 9, 1, -2, 3, 10, 6, 10, 14],
                                     [0, 0, 0, 0, 0, 0, -3, 0, 2, -3, -17, 0, 0, 4, 34, -2, 1, 10, -1, -10, 13],
                                     [10, 0, 0, 0, 0, 0, -1, -1, 10, -9, 7, 0, 0, -3, 1, 4, -5, -16, -6, -3, -24],
                                     [4, 0, 0, 0, 0, 0, -1, 1, 3, -3, -11, -15, -3, -2, 2, 7, -9, -12, 8, 6, -12],
                                     [0, 0, 0, 0, 0, 0, 0, 0, 12, -1, 0, 1, 0, -14, 2, 0, -11, -15, -15, -16, 9],
                                     [6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 41, 1, -1, -11, -15, 5, 1, -13],
                                     [0, 0, 0, 0, 0, -1, 47, 44, 1, -7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 1, -2, -2, 17, -1, 47, -10, -5, -6, -1, 0],
                                     [1, 0, 0, 0, 0, 0, -1, -2, 22, 18, 1, 33, -1, -1, 0, 8, 0, 0, 11, 0, -1],
                                     [0, 0, 0, 0, 0, 0, -35, 43, 2, -2, 0, 0, 0, -5, -1, -11, 1, 1, 0, 0, 1],
                                     [0, 0, 0, 0, 0, 0, -2, -1, -6, 50, 1, -28, 0, 0, 0, -7, 0, 0, -4, -1, 0]])
    read_matrix = struct_water3.veda_matrix.iloc[:, :-1].to_numpy()
    if np.array_equal(read_matrix, water_control_matrix):
        print('Water3 test OK!')
    else:
        print('Water3 test error')
        print(read_matrix)


if __name__ == '__main__':
    veda_matrix_test()
