glycerol_atoms = {'C': 3, 'H': 8, 'O': 3}
water_atoms = {'H': 2, 'O': 1}
pc_atoms = {'C': 5, 'H': 14, 'N': 1, 'P': 1, 'O': 4}
pe_atoms = {'C': 2, 'H': 8, 'N': 1, 'P': 1, 'O': 4}
ps_atoms = {'C': 3, 'H': 8, 'N': 1, 'P': 1, 'O': 6}
pi_atoms = {'C': 6, 'H': 13, 'P': 1, 'O': 9}


class sph_base:
    def __init__(self):
        self.formula = None

    def from_chain(self, chain):
        import re
        from loguru import logger
        match = re.fullmatch(r'([mdt])(\d+):(\d+)', chain)
        if match:
            hydroxyl_tag, carbon_atoms, db_nums = match.group(1), int(match.group(2)), int(match.group(3))
        else:
            logger.error('Invalid chain format: {}', chain)
            return None

        self.formula = {'C': carbon_atoms, 'H': 2 * carbon_atoms + 3 - 2 * db_nums, 'O': 0, 'N': 1}
        if hydroxyl_tag == 'm':
            self.formula['O'] += 1
        elif hydroxyl_tag == 'd':
            self.formula['O'] += 2
        elif hydroxyl_tag == 't':
            self.formula['O'] += 3

        # C18H39NO2
        return self.formula


sph = sph_base()

__all__ = ['glycerol_atoms', 'water_atoms', 'pc_atoms', 'pe_atoms', 'ps_atoms', 'pi_atoms', 'sph']
