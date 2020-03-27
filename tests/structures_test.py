import unittest
import warnings
from utilities.freq_comparator import freq_comparator
from create_structures import create_structures
from load_structures import load_structures

warnings.filterwarnings("ignore")

# все структуры, можно только добавлять новые
# structure_list = ['104', '104oh6', 'g104w3', '104oh6w6', 'goh3', 'gef3', 'ghole1', 'gep3oh6', 'gketoohw12',
#     'gketoohw1', 'water1', 'water3', 'lactone_lactol_w3', 'lactone_lactol_def', 'lactone_aromatic', 'ethylene_oxide',
#     'valerolactone', 'acetone', 'benzene', 'benzoic_acid', 'cyclohexanone', 'acetic_acid', 'isobutanol']
#
# оставлять только интересующие
structure_list = ['104']

create_structures(structure_list)


class TestStructure(unittest.TestCase):

    def test_create_structures(self):
        create_structures(structure_list)

    def test_load_structures(self):
        load_structures(structure_list)

    def test_get_report(self):
        structs = load_structures(structure_list)
        for structure in structs.values():
            for functional_group in structure.fg_assigned_matrix.columns:
                structure.get_report(functional_group)

    def test_freqs(self):
        for structure in structure_list:
            self.assertTrue(freq_comparator.compare_freq(structure, path='structures/', print_log=False) < 10,
                            msg=structure)


class TestPlotManager(unittest.TestCase):

    def test_spectrum_bar(self):
        structs = load_structures(structure_list)
        for structure in structs.values():
            structure.pm.spectrum_bar()

    def test_spectrum_curve(self):
        structs = load_structures(structure_list)
        for structure in structs.values():
            structure.pm.spectrum_curve()

    def test_fg_bar(self):
        structs = load_structures(structure_list)
        for structure in structs.values():
            for functional_group in structure.fg_assigned_matrix.columns:
                structure.pm.fg_bar(functional_group)

    def test_fg_sum_curve(self):
        structs = load_structures(structure_list)
        for structure in structs.values():
            for functional_group in structure.fg_assigned_matrix.columns:
                structure.pm.fg_sum_curve(functional_group)


if __name__ == '__main__':
    # create_structures(structure_list)
    unittest.main()
