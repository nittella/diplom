from constants import structure_list
from struc import Structure, get_number_of_atoms
from ped_manager import PEDManager, read_veda_matrix, read_freq, read_coord
import pickle


# Создает или пересоздает структуру и сохраняет ее в pickle файл в папаке структуры
def create_structures(names_of_structures=None, sf_function=None):
    ped_dict = dict()
    if names_of_structures is None:
        names_of_structures = structure_list
    for structure_name in names_of_structures:
        # for structure_name, coords_number in number_of_atoms.items():
        coords_number = get_number_of_atoms(structure_name)
        source_name = 'structures/' + structure_name + '/skra.'
        ped_dict[structure_name] = PEDManager(read_veda_matrix(
            file_name_=source_name + 'ved'),
            read_freq(file_name_=source_name + 'ved'),
            read_coord(coords_number, file_name_=source_name + 'dd2')
        )
    for structure_name in names_of_structures:
        with open('structures/' + structure_name + '/pickle', 'wb') as pickle_file:
            pickle.dump(Structure(structure_name, ped_dict[structure_name], sf_function=sf_function), pickle_file)

    print('New structures have been successfully created!')


if __name__ == '__main__':
    create_structures()
