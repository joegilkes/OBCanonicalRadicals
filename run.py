import sys
from openbabel import pybel
from openbabel import openbabel as ob

from OBCR import *

def main():
    molspec = sys.argv[1]
    if molspec.endswith('xyz'):
        mol = get_molecules(molspec)[0]
        initial_smi = pbmol_to_smi(mol)
    else:
        initial_smi = molspec
        mol = pybel.readstring('can', molspec)
        mol.addh()

    if is_radical(initial_smi):
        print('Type  idx  rads')
        for atom in ob.OBMolAtomIter(mol.OBMol):
            print(f'{atom.GetType()}    {atom.GetIdx()}    {get_radical_state(atom)}')
        print()
        mol = fix_radicals(mol)
        # Force re-parsing of structure to ensure aromaticity is detected.
        mol.addh()
        final_smi = pbmol_to_smi(mol)
    else:
        final_smi = initial_smi

    print(f'Initial: {initial_smi}')
    print(f'Final: {final_smi}')

if __name__ == '__main__':
    main()