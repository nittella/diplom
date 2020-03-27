from itertools import combinations


def generate_coord_set(func_groups):
    coords = set()
    for i in range(1, 5):
        for coord in combinations(sorted(func_groups), i):
            for coord_type in ['STRE', 'ANG']:
                coords.add('_'.join(coord) + '_' + coord_type)
    return coords
