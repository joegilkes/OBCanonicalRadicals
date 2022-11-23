from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np

def get_molecules(path):
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


def pbmol_to_smi(pbmol):
    '''Creates the Canonical SMILES representation of a given pybel Molecule.'''
    smi = pbmol.write('can').split()[0].strip()
    return smi


def is_radical(smi):
    '''Determines whether a given SMILES string contains radicals.'''
    if ('[' in smi):
        return True
    else:
        return False


def get_radical_state(obatom):
    '''Gets the radical state of a given OBAtom.'''
    typical_val = ob.GetTypicalValence(
        obatom.GetAtomicNum(), 
        obatom.GetTotalValence(),
        obatom.GetFormalCharge()
    )
    current_val = obatom.GetTotalValence()
    return typical_val - current_val

    
def get_hydrogenation(obatom, max_depth, curr_hydrog, curr_depth, prev_idx):
    '''Recurses through a molecule, evaluating hydrogenation from a given starting atom.
    
    Starting from the atom defined by `obatom`, recurses through a molecule
    via each atom's neighbours, up to a given depth. Evaluates overall
    hydrogenation of the molecule around this atom.

    Arguments:
        obatom: OBAtom object to evaluate hydrogenation around.
        max_depth: Maximum number of bonds to traverse out from.
        curr_hydrog: Current hydrogenation (initialise from 0).
        curr_depth: Current depth (initialise from 0).
        prev_idx: Index of previous atom in recursive call (initialise from obatom.GetIdx()).
    '''
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
    '''Find the radical atom to start canonical radical resolution from.
    
    In order to properly canonicalise the radical resolution process,
    resolution must start from the same atom each time. If multiple
    atoms share the same maximum radical state of the molecule, the
    canonical starting point is determined by which atom has the
    greatest neighbouring hydrogenation.

    This can be found be recursively exploring the neighbours of the
    target atoms. In cases where neighbouring hydrogenation is equal,
    the 'depth' of the search is increased, up to a maximum depth 
    where it can be safely concluded that the proposed start points
    have sufficiently similar environments that there is some symmetry,
    and the position of the start point will not matter.

    Arguments:
        targets: List of OBAtom indeces for potential starting atoms.
        obmol: OBMol object. 
    '''
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
    '''Recurses through a molecule, resolving dangling radical bonds.
    
    Starting from the atom defined by `obatom`, recurses through a
    molecule's bonds. If a bond is present between two atoms with
    free radicals, increments the bond order of this bond by one,
    implicitly removing these radicals.

    Returns a boolean indicating if any bonds have been changed
    by one full recursion over the molecule, such that the function
    can be continuously called until no more bonds are being changed.

    Arguments:
        obatom: OBAtom object to start resolving radicals around.
        prev_idx: Index of previous atom in recursive call (initialise from obatom.GetIdx()).
        start_idx: Index of atom which the function was originally called from (initialise from obatom.GetIdx()).
        bonds_changed: Boolean flag to keep track of bond change status through recursion (initialise as False).
        start_direction (Default = None): Optional argument, denoting which neighbour to recurse down first.
    '''
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
    '''Determines the best radical resolution from a given list.
    
    Given multiple radical resolutions of the same molecule (e.g. from
    different starting directions or atoms), throws out any invalid
    structures and determines the best by whichever has the highest
    overall bond order.

    Arguments:
        pbmols: List of Pybel Molecule objects.
    '''
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
    '''Canonicalise radical molecules within OpenBabel.
    
    OpenBabel sometimes struggles to parse structures with neighbouring
    radicals consistently, leading to multiple interpretations of radical
    structure coming from very similar geometries.

    This tries to fix the issue by enforcing that all neighbouring radicals
    should join together to form a bond, i.e. by transforming all radicals
    into their most stable state.

    Examples:
        * The species [CH2][C]C and C[C]=C have equal likelihood of being
        detected from the same geometry (via `pybel.readfile()`). These
        species have the same atoms, but differ in their bonding with the
        former having a radical CH2 and a diradical C, and the latter
        having only a radical C. This function resolves the discrepancy
        by forming a bond between the neighbouring radical carbons in the
        former, transforming it into the latter and thus canonicalising
        the radical structure.
        * If the cyclic molecule CC1C[CH][CH][C]1 were read in, it could be
        simplified to either CC1C[CH]C=[C]1 or CC1CC=C[C]1, depending on where
        the simplification of the radicals starts from and which direction
        around the ring it goes. This enforces a set of rules such that only
        CC1C[CH]C=[C]1 will be output every time as the canonical radical
        structure from that geometry.

    '''
    obmol = pbmol.OBMol

    # Find radical states of all atoms.
    radical_states = {}
    for atom in ob.OBMolAtomIter(obmol):
        atom_idx = atom.GetIdx()
        atom_rad_state = get_radical_state(atom)
        radical_states[atom_idx] = atom_rad_state

    # Check if radical resolution is worth doing.
    # If there is only one (or less) radical atom, no new bonds will
    # be formed so just return the original Molecule object.
    rad_vals = list(radical_states.values())
    if np.count_nonzero(rad_vals) <= 1:
        return pbmol

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