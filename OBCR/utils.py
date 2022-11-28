from typing import List

from openbabel import pybel
from openbabel import openbabel as ob

def get_molecules(path: str) -> List[pybel.Molecule]:
    '''Reads in molecule(s) from an xyz file.
    
    Assumes multiple molecules may be in the same xyz file, therefore
    returns a list of pybel Molecule objects.
    '''
    gen = pybel.readfile('xyz', path)
    mol = []
    gen_stat = True
    while gen_stat:
        try:
            next_mol = next(gen)
            mol.append(next_mol)
        except StopIteration:
            gen_stat = False
    return mol


def pbmol_to_smi(pbmol: pybel.Molecule) -> str:
    '''Creates the Canonical SMILES representation of a given pybel Molecule.'''
    smi = pbmol.write('can').split()[0].strip()
    return smi


def obmol_to_smi(obmol: ob.OBMol) -> str:
    '''Creates the Canonical SMILES representation of a given OpenBabel OBMol.'''
    pbmol = pybel.Molecule(obmol)
    return pbmol_to_smi(pbmol)


def is_radical(smi: str) -> bool:
    '''Determines whether a given SMILES string contains radicals.'''
    hydrogens = ['[H]', '[H][H]']
    if smi not in hydrogens and ('[' in smi):
        return True
    else:
        return False


def get_radical_state(obatom: ob.OBAtom) -> int:
    '''Gets the radical state of a given OBAtom.'''
    typical_val = ob.GetTypicalValence(
        obatom.GetAtomicNum(), 
        obatom.GetTotalValence(),
        obatom.GetFormalCharge()
    )
    current_val = obatom.GetTotalValence()
    return typical_val - current_val