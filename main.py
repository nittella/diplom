from create_structures import create_structures
from load_structures import load_structures
from read_spectrum_curve import Spectrum

# matplotlib.use('TkAgg')
# pd.set_option('display.max_columns', 1000)  # or 1000
# pd.set_option('display.max_rows', 1000)  # or 1000
# pd.set_option('display.max_colwidth', -1)  # or 199


if __name__ == '__main__':
    struct_name = 'salicylic_acid'
    create_structures([struct_name])
    # create_structures()
    structs = load_structures([struct_name])
    structs[struct_name].pm.spectrum_curve(sigma=20)
    structs[struct_name].pm.spectrum_bar()
    structs[struct_name].pm.fg_bar(['ketone'])
    Spectrum('salicylic_acid.CSV').to_integrated_mol_absorptivity_full_spectra().show()
