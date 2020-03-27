from matplotlib import pyplot as plt
import seaborn as sns
from create_structures import create_structures
from load_structures import load_structures


# Строит на одном графике
def many_plots(structures_list,  func, fg_list=None, axis=plt, color_palette=None):
    if fg_list is None:
        fg_number = 1
    else:
        fg_number = len(fg_list)

    if color_palette is None:
        length = len(structures_list) * fg_number
        color_palette = sns.color_palette("husl", length)

    for structure in structures_list.values():
        if fg_list is not None:
            # np.random.RandomState(seed=0)
            structure.pm.getattr(func)(fg_list, axis, color_palette[:fg_number])
            color_palette = color_palette[fg_number:]
        else:
            # np.random.RandomState(seed=0)
            structure.pm.getattr(func)(axis, color_palette[0])
            color_palette = color_palette[1:]


if __name__ == '__main__':
    structure_list = ['2_benzoyl_benzoic_acid', '3_hydroxy_benzoic_acid']
    create_structures(structure_list)
    structs = load_structures(structure_list)
    many_plots(structs, 'spectrum_curve')
    plt.legend()
    plt.show()
