from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np
import sys

def get_molecules(path):
    rgen = pybel.readfile('xyz', path)
    rmol = []
    gen_stat = True
    while gen_stat:
        try:
            rm = next(rgen)
            rmol.append(rm)
        except StopIteration:
            gen_stat = False
    return rmol

def pbmol_to_smi(pbmol):
    smi = pbmol.write('can').split()[0].strip()
    return smi

def is_radical(smi):
    if ('[' in smi):
        return True
    else:
        return False

def get_radical_state(obatom):
    typical_val = ob.GetTypicalValence(
        obatom.GetAtomicNum(), 
        obatom.GetTotalValence(),
        obatom.GetFormalCharge()
    )
    current_val = obatom.GetTotalValence()
    return typical_val - current_val


def get_hydrogenation(obatom, max_depth, curr_hydrog, curr_depth, prev_idx):
    # Check we aren't over-iterating.
    curr_depth += 1
    print(f'Current atom idx = {obatom.GetIdx()}')
    print(f'Current depth = {curr_depth}')
    if curr_depth >= max_depth:
        print(f'Max depth reached, hydrogenation = {curr_hydrog}')
        return curr_hydrog

    # Recurse into neighbours.
    obneighbours = []
    for neigh in ob.OBAtomAtomIter(obatom):
        # Ignore previously visited neighbours.
        if neigh.GetIdx() == prev_idx:
            continue
        # Increment on hydrogen neighbours.
        elif neigh.GetType() == 'H':
            print('Found hydrogen.')
            curr_hydrog += 1
        # Recurse for any continuations of the chain.
        else:
            print(f'Found non-hydrogen neighbour at index {neigh.GetIdx()}')
            obneighbours.append(neigh)

    if len(obneighbours) == 0:
        print('Reached end of chain.')
        return curr_hydrog
    else:
        for neigh in obneighbours:
            print(f'Recursing into neighbour at index {neigh.GetIdx()}')
            curr_hydrog = get_hydrogenation(neigh, max_depth, curr_hydrog, curr_depth, obatom.GetIdx())

    return curr_hydrog


def find_starting_radical(targets, obmol):
    start_found = False
    curr_depth = 1
    max_depth = 7
    while not start_found:
        curr_depth += 2
        hyd = [0 for _ in range(len(targets))]
        for i, targ_idx in enumerate(targets):
            print(f'### Exploring around atom with index {targ_idx} ###')
            hyd[i] = get_hydrogenation(obmol.GetAtom(targ_idx), curr_depth, 0, 0, targ_idx)
        
        # If there is only one maximum, the start point has been found.
        if hyd.count(max(hyd)) == 1:
            start_found = True
            start_atom = targets[np.argmax(hyd)]
            print(f'Start at atom with index {start_atom}\n')
            return start_atom
        # If there are multiple maxima, check against max depth.
        else:
            if curr_depth == max_depth:
                print('Max depth reached, starting from any of the targets should be valid.\n')
                return targets[0]
            # Shortlist targets with joint max hydrogenation and iterate again.
            else:
                targets = [targets[i] for i in np.argwhere(hyd == np.amax(hyd)).flatten()]
                print(f'Multiple targets found: {targets}')
                print('Retrying at greater depth.\n')


def resolve_radicals(obatom, prev_idx, start_idx, bonds_changed, start_direction=None):
    obneighidx = []
    for neigh in ob.OBAtomAtomIter(obatom):
        if neigh.GetType() == 'H':
            continue
        elif neigh.GetIdx() == prev_idx:
            continue
        elif neigh.GetIdx() == start_idx:
            continue
        else:
            obneighidx.append(neigh.GetIdx())

    if len(obneighidx) == 0:
        print('Finished exploring chain.')
        return bonds_changed

    if start_direction is not None:
        obneighidx.pop(np.where(np.array(obneighidx) == start_direction)[0][0])
        obneighidx.insert(0, start_direction)
    
    obneighbours = [obatom.GetParent().GetAtom(idx) for idx in obneighidx]
    
    for neigh in obneighbours:
        if (get_radical_state(obatom) > 0) and (get_radical_state(neigh) > 0):
            bond = obatom.GetBond(neigh)
            bond.SetBondOrder(bond.GetBondOrder()+1)
            bonds_changed = True
        bonds_changed = resolve_radicals(neigh, obatom.GetIdx(), start_idx, bonds_changed)

    return bonds_changed


def get_best_resolution(pbmols):
    obmols = [pbmol.OBMol for pbmol in pbmols]

    # Check for invalid structures, e.g. allenes within rings.
    invalid_mols = []
    for i, obmol in enumerate(obmols):
        rings = pbmols[i].sssr
        n_rings = len(rings)
        if n_rings == 0:
            continue
        else:
            allene_in_ring = False
            for ring in rings:
                if allene_in_ring: continue
                ringpath = ring._path
                ringsize = len(ringpath)
                ringpath = [idx for idx in ringpath]
                ringpath.extend([ringpath[0], ringpath[1]])

                prev_double = False
                for j in range(ringsize+1):
                    bond = obmol.GetBond(ringpath[j], ringpath[j+1])
                    if bond.GetBondOrder() > 1:
                        if prev_double == True:
                            allene_in_ring = True
                            break
                        else:
                            prev_double = True
                    else:
                        prev_double = False

            if allene_in_ring:
                print(f'Allene found within ring in structure {i+1}')
                invalid_mols.append(i)


    # Remove invalid structures.
    for i in sorted(invalid_mols, reverse=True):
        pbmols.pop(i)
        obmols.pop(i)

    if len(obmols) == 0:
        raise RuntimeError('No valid structure found for radical resolution of molecule.')
    elif len(obmols) == 1:
        return pbmols[0]
    # If there is still a choice of valid configurations, choose
    # the one with the greatest overall bond order.
    else:
        total_bond_orders = [0 for _ in range(len(obmols))]
        for i, obmol in enumerate(obmols):
            for bond in ob.OBMolBondIter(obmol):
                total_bond_orders[i] += bond.GetBondOrder()

        if total_bond_orders.count(max(total_bond_orders)) == 1:
            best_pbmol = pbmols[np.argmax(total_bond_orders)]
        else:
            best_pbmol = pbmols[0]

        return best_pbmol
        



def fix_radicals(pbmol):
    obmol = pbmol.OBMol

    # Find radical states of all atoms.
    radical_states = {}
    for atom in ob.OBMolAtomIter(obmol):
        atom_idx = atom.GetIdx()
        atom_rad_state = get_radical_state(atom)
        radical_states[atom_idx] = atom_rad_state

    # Select starting point.
    max_radical_state = max(radical_states.values())
    highest_radical_atoms = [k for k, v in radical_states.items() if v == max_radical_state]
    # If multiple starting points, find the canonical start point.
    if len(highest_radical_atoms) > 1:
        start_idx = find_starting_radical(highest_radical_atoms, obmol)
    else:
        start_idx = highest_radical_atoms[0]

    print(f'Starting radical resolution from atom with index {start_idx}')
    obneighidx = []
    for neigh in ob.OBAtomAtomIter(obmol.GetAtom(start_idx)):
        if neigh.GetType() == 'H':
            continue
        else:
            obneighidx.append(neigh.GetIdx())
    n_neighbours = len(obneighidx)

    if n_neighbours > 1:
        print(f'Starting atom has {n_neighbours} neighbours, resolving best direction.')
        pbmols = [pbmol.clone for _ in range(n_neighbours)]
        obmols = [mol.OBMol for mol in pbmols]
        for i, obmol in enumerate(obmols):
            start_atom = obmol.GetAtom(start_idx)
            bonds_changed = True
            while bonds_changed:
                bonds_changed = resolve_radicals(start_atom, start_idx, start_idx, False, start_direction=obneighidx[i])

        pbmol = get_best_resolution(pbmols)

    else:
        start_atom = obmol.GetAtom(start_idx)
        bonds_changed = True
        while bonds_changed:
            bonds_changed = resolve_radicals(start_atom, start_idx, start_idx, False)

    return pbmol