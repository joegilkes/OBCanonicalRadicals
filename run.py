import sys
from openbabel import pybel
from openbabel import openbabel as ob

from OBCR import get_molecules, pbmol_to_smi, obmol_to_smi, is_radical, get_radical_state, fix_radicals

def main():
    molspec = sys.argv[1]
    if molspec.endswith('xyz'):
        mol = get_molecules(molspec)[0]
        initial_smi = pbmol_to_smi(mol)
    else:
        initial_smi = molspec
        mol = pybel.readstring('can', molspec)
        mol.addh()

    if '.' in initial_smi:
        species = mol.OBMol.Separate()
    else:
        species = [mol.OBMol]

    final_smi = []
    for i, spec in enumerate(species):
        smi = obmol_to_smi(spec)
        if is_radical(smi):
            print(f'Species {i+1}')
            print('Type  idx  rads')
            for atom in ob.OBMolAtomIter(mol.OBMol):
                print(f'{atom.GetType()}    {atom.GetIdx()}    {get_radical_state(atom)}')
            print()
            spec_mol = pybel.Molecule(spec)
            spec_mol = fix_radicals(spec_mol)
            # Force re-parsing of structure to ensure aromaticity is detected.
            spec_mol.addh()
            final_smi.append(pbmol_to_smi(spec_mol))
        else:
            final_smi.append(smi)

    if len(final_smi) == 1:
        final_smi = final_smi[0]
    else:
        final_smi = '.'.join(final_smi)

    print(f'Initial: {initial_smi}')
    print(f'Final: {final_smi}')

if __name__ == '__main__':
    main()