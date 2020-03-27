# Какие координаты к каким функциональным круппам относятся
func_group_to_coords = {
    'hydroxide': [
        'hydroxide_STRE_assignment',
        'carbon_hydroxide_STRE_assignment',
        'carbon_hydroxide_ANG_assignment'
    ],
    'water': [
        'water_STRE_assignment',
        'water_ANG_assignment',

    ],
    'epoxide': [
        'carbon_epoxide_STRE_assignment',
        'carbon_epoxide_ANG_assignment'
    ],
    'carbon': [
        'carbon_STRE_assignment',
        'carbon_ANG_assignment'
    ],
    'ketone': [
        'carbon_ketone_STRE_assignment',
        'carbon_ketone_ANG_assignment'
    ],
    'carboxyl': [
        # 'carbon_carboxyl_STRE_assignment',
        # 'carbon_carboxyl_ANG_assignment',
        'carboxyl_STRE_assignment'
    ],
    'carboxyl_hydroxide': [
        'carbon_carboxyl_hydroxide_ANG_assignment'
    ],
    'hydrogen': [
        'carbon_hydrogen_STRE_assignment',
        'carbon_hydrogen_ANG_assignment'
    ],
    'hydroxide_water': [
        'hydroxide_water_ANG_assignment',
        'hydroxide_water_STRE_assignment',
        'carbon_hydroxide_water_ANG_assignment'
    ],
    'lactone': [
        'lactone_ANG_assignment',
        'lactone_STRE_assignment'
    ],
    'lactol': [
        'lactol_ANG_assignment',
        'lactol_STRE_assignment'
    ]
}

# Цвета определяют каким цыетом будет отрисован спектр для соответвующей функционльной группы
fuc_group_colors = {
    'hydroxide': '#2b2be5',
    'water': '#1fb69b',
    'epoxide': '#db8d16',
    'carbon': '#000000',
    'ketone': '#0fd831',
    'carboxyl': '#ff0000',
    'carboxyl_hydroxide': '#850fd8',
    'hydrogen': '#89dfec',
    'hydroxide_water': '#33C9FF',
    'lactone': '#f65ef0',
    'lactol': '#CC66FF'
}

atoms_to_hb = {
    'goh3': {
        77: True
    },
    '104oh6w6': {
        78: True,
        80: True,
        76: True,
        77: True,
        71: True
    },
    '104oh6': {
        77: True,
        86: True,
        81: True
    }
}

atoms_to_epox = {
    'goh3': {},
    '104oh6w6': {},
    '104oh6': {
        78: True,
        80: True,
        76: True,
        82: True,
        84: True,
        87: True
    }
}