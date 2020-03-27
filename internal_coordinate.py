class Coordinate:
    def __init__(self, coefficient, atoms, name, value):
        self.coefficient = coefficient
        self.atoms = atoms
        self.name = name
        self.value = value

    def __str__(self):
        return 'coef.: ' + str(self.coefficient) + ' atoms: ' + str(self.atoms) \
               + ' name: ' + self.name + ' val.: ' + str(self.value) + '\n'


class ComplexCoordinate:
    def __init__(self, coord_type, number, coordinates):
        self.type = coord_type
        self.coordinates = coordinates
        self.number = number

    def __str__(self):
        return '№' + str(self.number) + ' type: ' + self.type + ' coord.:\n' + '\n'.join(map(str, self.coordinates))

    def add_coordinate(self, coordinate):
        self.coordinates.append(coordinate)