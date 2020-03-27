class InternalCoordinatesGenerator:
    def __init__(self):
        self._circle_connection = dict([
            (1, 0),
            (2, 1), (3, 1),
            (4, 2), (5, 3),
            (6, 4), (7, 4), (8, 5),
            (9, 6), (10, 7), (11, 8)
        ])
        self._circle_dict = dict([
            (0, [21, 75, 31, 32, 33, 22]),
            (1, [12, 20, 30, 42, 34, 23]),
            (2, [11, 19, 40, 43, 35, 14]),
            (3, [13, 10, 29, 41, 44, 24]),
            (4, [3, 18, 39, 51, 36, 15]),
            (5, [5, 9, 28, 49, 45, 25]),
            (6, [4, 8, 38, 50, 46, 16]),
            (7, [2, 17, 47, 52, 37, 7]),
            (8, [6, 1, 27, 48, 53, 26]),
            (9, [66, 63, 60, 57, 54, 69]),
            (10, [65, 62, 59, 56, 71, 68]),
            (11, [67, 64, 61, 58, 55, 70])
        ])

    def generate_hh_bonds(self):
        for (i, circle) in self._circle_dict.items():
            for j, atom in enumerate(circle):
                prev_atom = circle[(j - 1) % 6]
                next_atom = circle[(j + 1) % 6]
                if i > 8:
                    print(' 3 . Ru ' + str(atom) + ' : C ' +
                          str(self._circle_dict[self._circle_connection[i]][j])
                          + ' , Ru ' + str(prev_atom)
                          + ' ,  Ru ' + str(next_atom))
                elif i > 0:
                    print(' 3 . C ' + str(atom) + ' : C ' +
                          str(self._circle_dict[self._circle_connection[i]][j])
                          + ' , C ' + str(prev_atom)
                          + ' ,  C ' + str(next_atom))
                else:
                    print(' 3 . C ' + str(atom) + ' : C ' +
                          str(self._circle_dict[1][j]) + ' , C ' + str(
                        prev_atom) + ' ,  C ' + str(next_atom))


InternalCoordinatesGenerator().generate_hh_bonds()
