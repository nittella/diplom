import unittest
from generate_coord_set import generate_coord_set


class TestGenerateCoordsSet(unittest.TestCase):
    def test(self):
        functional_groups = ['oh', 'cooh', 'co', 'coc']
        result = {'co_ANG', 'co_STRE', 'co_coc_ANG', 'co_coc_STRE', 'co_coc_cooh_ANG', 'co_coc_cooh_STRE', 'co_coc_cooh_oh_ANG', 'co_coc_cooh_oh_STRE', 'co_coc_oh_ANG', 'co_coc_oh_STRE', 'co_cooh_ANG', 'co_cooh_STRE', 'co_cooh_oh_ANG', 'co_cooh_oh_STRE', 'co_oh_ANG', 'co_oh_STRE', 'coc_ANG', 'coc_STRE', 'coc_cooh_ANG', 'coc_cooh_STRE', 'coc_cooh_oh_ANG', 'coc_cooh_oh_STRE', 'coc_oh_ANG', 'coc_oh_STRE', 'cooh_ANG', 'cooh_STRE', 'cooh_oh_ANG', 'cooh_oh_STRE', 'oh_ANG', 'oh_STRE'}
        self.assertEqual(generate_coord_set(functional_groups), result)


if __name__ == '__main__':
    unittest.main()
