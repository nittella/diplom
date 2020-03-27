from constants import structure_list
import pickle


# Загружает структуру из pickle файла. Чтобы создать pickle файл необходимо запустить функцию create_structures
# из create_structures.py
def load_structures(structure_names=None):
    structs = []

    if structure_names is None:
        structure_names = structure_list

    for struct_name in structure_names:
        with open('structures/' + struct_name + '/pickle', 'rb') as pickle_file:
            structs.append(pickle.load(pickle_file))
    print('Structure(s)', ' '.join(structure_names), 'have been successfully loaded!' )
    return {s.name: s for s in structs}
