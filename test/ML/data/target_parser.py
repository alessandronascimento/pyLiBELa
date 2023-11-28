from pathlib import Path
from itertools import islice

class Parser():

    @staticmethod
    def parse_index_file(path : str | Path) -> dict[str, float] :
        out_dict = {}
        replace_map = (('mM', 'uM', 'nM', 'pM', 'fM'),
                       ('E-3', 'E-6', 'E-9', 'E-12', 'E-15'))

        with open(path, 'r') as file: 
            for line in islice(file, 6, None):
                name, k = line.split()[0], line.split()[3]
                if 'IC50' in k: continue
                val_str = k.strip("Kdi<>=")
                out_dict[name] = sum(map(lambda unit, value: eval(val_str.strip(unit) + value)  if unit in val_str else 0,  replace_map[0], replace_map[1]))

        return out_dict




