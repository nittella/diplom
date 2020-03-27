from constants import func_groups

path = '../../structures/'


def get_atoms_type(structure_name):
    atoms = [('carbon', 'C')] * (func_groups.number_of_atoms[structure_name] + 1)
    for i in range(54, 72):
        atoms[i] = ('hydrogen', 'H')
    for i, atom in enumerate(atoms):
        if i in func_groups.func_group_in_structures[structure_name]:
            atoms[i] = func_groups.func_group_in_structures[structure_name][i]
    return atoms


def generate_mpo(structure_name):
    bonds = [[i] for i in range(func_groups.number_of_atoms[structure_name] + 1)]
    bonds_path = path + '/' + structure_name + '/bonding_data'
    with open(bonds_path) as input_file:
        lines = [line.strip('\n').split() for line in input_file.readlines()]
        for line in lines:
            if len(line) < 3 or line[2] != 'H':
                atom1 = int(line[0])
                atom2 = int(line[1])
                bonds[atom1].append(atom2)
                bonds[atom2].append(atom1)
    bonds.sort(key=lambda x: -len(x))

    atom_types = get_atoms_type(structure_name)
    with open('output.mpo', 'w') as output_file:
        for line in bonds[:-1]:
            output_file.write(' ')
            output_file.write(str(len(line) - 1))
            output_file.write(' . ')
            output_file.write(str(atom_types[line[0]][1]) + ' ' + str(line[0]))
            output_file.write(' : ')
            tmp = [str(atom_types[elem][1]) + ' ' + str(elem) for elem in line[1:]]
            output_file.write(' , '.join(tmp))
            output_file.write('\n')
    print('for struture is', structure_name, 'done!')


if __name__ == '__main__':
    generate_mpo('gketoohw1')