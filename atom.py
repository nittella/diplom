from abc import ABC, abstractmethod


class Atom:
    def __init__(self, bound_atoms=None, functional_group='unknown'):
        if bound_atoms is None:
            bound_atoms = []

        self.functional_group = functional_group
        self.bound_atoms = bound_atoms  # Хранит ссылки на связанные атомы

    def update_simple_functional_group(self):
        pass

    def update_complex_functional_group(self):
        pass

    def update_functional_group(self):
        pass


class CarbonAtom(Atom):
    def update_simple_functional_group(self):
        self.functional_group = 'carbon'

    def update_complex_functional_group(self):
        oxygen_atoms = [bound_atom for bound_atom in self.bound_atoms if isinstance(bound_atom, OxygenAtom)]
        if len(oxygen_atoms) == 2:

            rules = [
                (('ketone', 'hydroxide'), 'carboxyl'),
                (('ketone', 'epoxide'), 'lactone'),
                (('ketone', 'epoxide'), 'lactol')
            ]

            for rule, complex_functional_group in rules:
                if set([oxygen_atom.functional_group for oxygen_atom in oxygen_atoms]) == set(rule):
                    for atom in [self, *oxygen_atoms]:
                        atom.functional_group = complex_functional_group
                    break


def check_type(atoms, types):
    return set([type(atom) for atom in atoms]) == set(types)


class OxygenAtom(Atom):
    def update_simple_functional_group(self):
        rules = [
            (len(self.bound_atoms) == 1, 'ketone'),
            (check_type(self.bound_atoms, [HydrogenAtom, HydrogenAtom]), 'water'),
            (check_type(self.bound_atoms, [CarbonAtom, CarbonAtom]), 'epoxide'),
            (check_type(self.bound_atoms, [HydrogenAtom, CarbonAtom]), 'hydroxide')
        ]

        for rule, simple_functional_group in rules:
            if rule:
                self.functional_group = simple_functional_group
                break


class HydrogenAtom(Atom):
    def update_functional_group(self):
        if self.bound_atoms[0].functional_group != 'carbon':
            self.functional_group = self.bound_atoms[0].functional_group
        else:
            self.functional_group = 'hydrogen'


atom_number_to_class = {
    '0': Atom,
    '1': HydrogenAtom,
    '6': CarbonAtom,
    '8': OxygenAtom
}


def create_atom(atom_number, bound_atoms=None, functional_group='unknown'):
    return atom_number_to_class[atom_number](bound_atoms, functional_group)


def update_atoms_fg_type_all(atom_list):
    #  Oпределяем к какой простой функциональной группе принадлежит атом кислорода (ketone, water, epoxide,
    #  hydroxide). Все атомы углерода относим к классу carbon, все атомы водорода к классу hydrogen.
    for atom in atom_list:
        atom.update_simple_functional_group()

    # Определяем сложные функциональные группы (carboxyl, lactone, lactol)
    for atom in atom_list:
        atom.update_complex_functional_group()

    # Относим атомы водорода к тому же классу, что и связанный с ними атом.
    for atom in atom_list:
        atom.update_functional_group()


# Читает aтомные номера из файла
def get_atoms_number(path):
    path = path + 'coords'
    with open(path, 'r') as input_file:
        result = ['0']
        lines = [line.strip('\n').split() for line in input_file.readlines()]
        for line in lines:
            result.append(line[0])
        return result


# Читает связанные атомы из файла
def get_bound_atoms(path, number_of_atoms):
    path = path + 'bonding_data'
    with open(path) as input_file:
        lines = [line.strip('\n').split() for line in input_file.readlines()]
        bonds = [[] for _ in range(number_of_atoms + 1)]
        for line in lines:
            if len(line) < 3 or line[2] != 'H':
                atom1 = int(line[0])
                atom2 = int(line[1])
                bonds[atom1].append(atom2)
                bonds[atom2].append(atom1)
        return bonds


# Создает лист из связанных между собой атомов(по сути граф)
def read_atoms(path):
    # Сохраняем нумерацию атомов, поэтому первый атом пустой
    atoms = []
    atom_numbers = get_atoms_number(path)
    number_of_atoms = len(atom_numbers)
    bound_atoms = get_bound_atoms(path, number_of_atoms)
    # atoms_fg_types = ['unknown'] * (number_of_atoms + 1)
    # создаем атомы и добавляем их в список атомов
    for i in range(number_of_atoms):
        new_atom = create_atom(atom_numbers[i])
        atoms.append(new_atom)

    # создаем ссылки в atom.bound_atoms на связанные атомы
    for i in range(number_of_atoms):
        atoms[i].bound_atoms = [atoms[k] for k in bound_atoms[i]]

    update_atoms_fg_type_all(atoms)

    return atoms

