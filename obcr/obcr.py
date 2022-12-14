from typing import Union
import logging as log

from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np

from obcr.utils import *


class HydrogenationResolver:
    def __init__(self, obatom: ob.OBAtom, max_depth: int):
        '''Recurses through a molecule, evaluating hydrogenation from a given starting atom.
    
        Starting from the atom defined by `obatom`, recurses through a molecule
        via each atom's neighbours, up to a given depth. Evaluates overall
        hydrogenation of the molecule around this atom.

        Arguments:
            obatom: OBAtom object to evaluate hydrogenation around.
            max_depth: Maximum number of bonds to traverse out from.
        '''
        self.obatom = obatom
        self.max_depth = max_depth
        self.curr_hydrog = 0
        self.prev_idx = self.obatom.GetIdx()
        self.prev_idxs = [self.prev_idx]

    def __call__(self, obatom: Union[ob.OBAtom, None]=None,
                 prev_idx: Union[int, None]=None, 
                 curr_depth: Union[int, None]=None):
        if obatom is None:
            obatom = self.obatom
        if prev_idx is None:
            prev_idx = self.prev_idx
        else:
            self.prev_idxs.append(prev_idx)
        if curr_depth is None:
            curr_depth = 0

        # Check we aren't over-iterating.
        curr_depth += 1
        log.info(f'Current atom idx = {obatom.GetIdx()}')
        log.info(f'Current depth = {curr_depth}')
        if curr_depth >= self.max_depth:
            log.info(f'Max depth reached, hydrogenation = {self.curr_hydrog}')
            return self

        # Recurse into neighbours.
        obneighbours = []
        for neigh in ob.OBAtomAtomIter(obatom):
            # Ignore previously visited neighbours.
            if neigh.GetIdx() in self.prev_idxs:
                continue
            # Increment on hydrogen neighbours.
            elif neigh.GetType() == 'H':
                log.info('Found hydrogen.')
                self.curr_hydrog += 1
            # Recurse for any continuations of the chain.
            else:
                log.info(f'Found non-hydrogen neighbour at index {neigh.GetIdx()}')
                obneighbours.append(neigh)

        if len(obneighbours) == 0:
            log.info('Reached end of chain.')
            return self
        else:
            for neigh in obneighbours:
                log.info(f'Recursing into neighbour at index {neigh.GetIdx()}')
                self.__call__(neigh, obatom.GetIdx(), curr_depth)

        return self


class RadicalResolver:
    def __init__(self, obatom: ob.OBAtom, start_direction: Union[int, None]=None):
        '''Recurses through a molecule, resolving dangling radical bonds.
    
        Starting from the atom defined by `obatom`, recurses through a
        molecule's bonds. If a bond is present between two atoms with
        free radicals, increments the bond order of this bond by one,
        implicitly removing these radicals.

        The `bonds_changed` field can be queried after running the function
        to determine if the radical structure has converged, or if another
        pass over the molecule should be done to continue resolving bonds.

        Arguments:
            obatom: OBAtom object to start resolving radicals around.
            start_direction (Default = None): Optional argument, denoting which neighbour to recurse down first.
        '''
        self.obatom = obatom
        self.obmol = obatom.GetParent()
        self.start_idx = obatom.GetIdx()
        self.prev_idx = obatom.GetIdx()
        self.prev_idxs = [self.start_idx]
        self.bonds_changed = False
        self.start_direction = start_direction

        rings = self.obmol.GetSSSR()
        self.n_rings = len(rings)
        if self.n_rings == 0:
            self.rings = None
            self.ringpaths = None
            self.ringsizes = None
        else:
            self.rings = rings
            self.ringpaths = []
            self.ringsizes = []
            for i, ring in enumerate(rings):
                self.ringpaths.append([idx for idx in ring._path])
                self.ringsizes.append(len(self.ringpaths[i]))


    def __call__(self, obatom: Union[ob.OBAtom, None]=None, prev_idx: Union[int, None]=None):
        if obatom is None:
            obatom = self.obatom
        if prev_idx is None:
            prev_idx = self.prev_idx
        else:
            self.prev_idxs.append(prev_idx)

        obneighidx = []
        for neigh in ob.OBAtomAtomIter(obatom):
            if neigh.GetType() == 'H':
                continue
            elif neigh.GetIdx() in self.prev_idxs:
                continue
            else:
                obneighidx.append(neigh.GetIdx())

        if len(obneighidx) == 0:
            log.info('Finished exploring chain.')
            return self

        if self.start_direction is not None:
            obneighidx.pop(np.where(np.array(obneighidx) == self.start_direction)[0][0])
            obneighidx.insert(0, self.start_direction)
            self.start_direction = None
        
        obneighbours = [self.obmol.GetAtom(idx) for idx in obneighidx]
        
        for neigh in obneighbours:
            if (get_radical_state(obatom) > 0) and (get_radical_state(neigh) > 0):
                # Check if atoms are in a ring.
                # This implicitly cannot include atoms in multiple rings due to valence constraints.
                ignore_bond = False
                if self.n_rings > 0 and obatom.IsInRing() and neigh.IsInRing():
                    ignore_bond = self._check_rings(obatom, neigh)
                # Also check if this bond should be ignored in favour of conjugation.
                elif obatom.GetIdx() != prev_idx:
                    prev_bond = self.obmol.GetBond(obatom.GetIdx(), prev_idx)
                    if prev_bond.GetBondOrder() > 1:
                        ignore_bond = self._check_conjugation(obatom, neigh)

                if not ignore_bond:        
                    bond = obatom.GetBond(neigh)
                    bond.SetBondOrder(bond.GetBondOrder()+1)
                    self.bonds_changed = True
            self.__call__(neigh, obatom.GetIdx())

        return self


    def _check_rings(self, obatom: ob.OBAtom, neigh: ob.OBAtom):
        '''Checks neighbouring bonds in rings to determine if a bond update should be ignored.
        
        Returns True if the bond update should be ignored, otherwise
        returns False.
        '''
        # Ensure no triple bonds can exist in rings.
        bond = self.obmol.GetBond(obatom.GetIdx(), neigh.GetIdx())
        if bond.GetBondOrder() >= 2:
            return True

        # Determine which ring the atoms are in (if applicable).
        ringnum = None
        if self.n_rings == 1:
            ringnum = 0
        else:
            for i, path in enumerate(self.ringpaths):
                if obatom.GetIdx() in path and neigh.GetIdx() in path:
                    ringnum = i
        if ringnum is None:
            raise RuntimeError('Atoms could not be found in any ring.')

        # Figure out if the ring path needs to be reversed based on the current direction of exploration.
        ringpath = np.array(self.ringpaths[ringnum])
        ringsize = self.ringsizes[ringnum]
        obatom_ringpos = np.where(ringpath == obatom.GetIdx())[0][0]
        neigh_ringpos = np.where(ringpath == neigh.GetIdx())[0][0]
        if obatom_ringpos == 0 and neigh_ringpos == ringsize-1:
            flipped = True
        elif obatom_ringpos == ringsize-1 and neigh_ringpos == 0:
            flipped = False
        elif obatom_ringpos > neigh_ringpos:
            flipped = True
        else:
            flipped = False
        if flipped:
            ringpath = np.flip(ringpath)
            obatom_ringpos = np.where(ringpath == obatom.GetIdx())[0][0]
            neigh_ringpos = np.where(ringpath == neigh.GetIdx())[0][0]

        # Check bond orders of previous and next bonds in ring.
        next_bond_atom_idx = ringpath.tolist()[(neigh_ringpos + 1) % ringsize]
        prev_bond_atom_idx = ringpath.tolist()[(obatom_ringpos - 1) % ringsize]
        next_bond = self.obmol.GetBond(neigh.GetIdx(), next_bond_atom_idx)
        prev_bond = self.obmol.GetBond(obatom.GetIdx(), prev_bond_atom_idx)
        # If either are > 1, do not apply a bond order change.
        if (next_bond.GetBondOrder() > 1) or (prev_bond.GetBondOrder() > 1):
            return True

        return False


    def _check_conjugation(self, obatom: ob.OBAtom, neigh: ob.OBAtom):
        '''Checks neighbouring atoms/bonds to determine if forming a bond would break conjugation.
        
        Returns True if the bond update should be ignored, otherwise
        returns False.
        '''
        obneighidx = []
        for next_neigh in ob.OBAtomAtomIter(neigh):
            if next_neigh.GetType() == 'H':
                continue
            elif next_neigh.GetIdx() == obatom.GetIdx():
                continue
            else:
                obneighidx.append(next_neigh.GetIdx())

        if len(obneighidx) == 0:
            return False

        for nnidx in obneighidx:
            # Don't form the bond if there is already conjugation in place.
            if self.obmol.GetBond(neigh.GetIdx(), nnidx).GetBondOrder() > 1:
                return True
            # Don't form the bond if skipping it would form a conjugated system.
            else:
                atom = self.obmol.GetAtom(nnidx)
                atom_rads = get_radical_state(atom)
                if atom_rads >= 1:
                    return True

        return False
        

def find_starting_radical(targets: List[int], obmol: ob.OBMol) -> int:
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
            log.info(f'### Exploring around atom with index {targ_idx} ###')
            hresolve = HydrogenationResolver(obmol.GetAtom(targ_idx), curr_depth)
            hyd[i] = hresolve().curr_hydrog
        
        # If there is only one maximum, the start point has been found.
        if hyd.count(max(hyd)) == 1:
            start_found = True
            start_atom = targets[np.argmax(hyd)]
            log.info(f'Start at atom with index {start_atom}\n')
            return start_atom
        # If there are multiple maxima, check against max depth.
        else:
            if curr_depth == max_depth:
                log.info('Max depth reached, starting from any of the targets should be valid.\n')
                return targets[0]
            # Shortlist targets with joint max hydrogenation and iterate again.
            else:
                targets = [targets[i] for i in np.argwhere(hyd == np.amax(hyd)).flatten()]
                log.info(f'Multiple targets found: {targets}')
                log.info('Retrying at greater depth.\n')


def get_best_resolution(pbmols: List[pybel.Molecule]):
    '''Determines the best radical resolution from a given list.
    
    Given multiple radical resolutions of the same molecule (e.g. from
    different starting directions or atoms), determines the best by 
    whichever has the highest overall bond order. If there is a tie,
    determines the best by whichever has the least free radical electrons.

    Arguments:
        pbmols: List of Pybel Molecule objects.
    '''
    obmols = [pbmol.OBMol for pbmol in pbmols]

    total_bond_orders = [0 for _ in range(len(obmols))]
    radical_elecs = [{1: 0, 2: 0, 3: 0} for _ in range(len(obmols))]
    for i, obmol in enumerate(obmols):
        for bond in ob.OBMolBondIter(obmol):
            total_bond_orders[i] += bond.GetBondOrder()
        for atom in ob.OBMolAtomIter(obmol):
            rad_state = get_radical_state(atom)
            if rad_state > 0:
                radical_elecs[i][rad_state] += 1

    if total_bond_orders.count(max(total_bond_orders)) == 1:
        best_pbmol = pbmols[np.argmax(total_bond_orders)]
    else:
        # Form subset of best molecules so far.
        best_molids = np.argwhere(total_bond_orders == np.amax(total_bond_orders)).flatten()
        best_radical_elecs = [radical_elecs[i] for i in best_molids]
        # If one molecule has a lower maximum radical state than the others, choose this molecule.
        state_3_counts = [rdict[3] for rdict in best_radical_elecs]
        state_2_counts = [rdict[2] for rdict in best_radical_elecs]
        state_1_counts = [rdict[1] for rdict in best_radical_elecs]
        if state_3_counts.count(min(state_3_counts)) == 1:
            best_pbmol = pbmols[best_molids[np.argmin(state_3_counts)]]
        elif state_2_counts.count(min(state_2_counts)) == 1:
            best_pbmol = pbmols[best_molids[np.argmin(state_2_counts)]]
        elif state_1_counts.count(min(state_1_counts)) == 1:
            best_pbmol = pbmols[best_molids[np.argmin(state_1_counts)]]
        else:
            best_pbmol = pbmols[best_molids[0]]

    return best_pbmol


def fix_radicals(pbmol: pybel.Molecule, verbose: bool=False):
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
    if verbose:
        log.basicConfig(level=log.INFO, format='%(message)s')

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

    log.info(f'Starting radical resolution from atom with index {start_idx}')
    obneighidx = []
    for neigh in ob.OBAtomAtomIter(obmol.GetAtom(start_idx)):
        if neigh.GetType() == 'H':
            continue
        else:
            obneighidx.append(neigh.GetIdx())
    n_neighbours = len(obneighidx)

    if n_neighbours > 1:
        log.info(f'Starting atom has {n_neighbours} neighbours, resolving best direction.')
        pbmols = [pbmol.clone for _ in range(n_neighbours)]
        obmols = [mol.OBMol for mol in pbmols]
        for i, obmol in enumerate(obmols):
            log.info(f'Going towards neighbour {obneighidx[i]}')
            start_atom = obmol.GetAtom(start_idx)
            bonds_changed = True
            while bonds_changed:
                rresolve = RadicalResolver(start_atom, start_direction=obneighidx[i])
                bonds_changed = rresolve().bonds_changed

        pbmol = get_best_resolution(pbmols)

    else:
        start_atom = obmol.GetAtom(start_idx)
        bonds_changed = True
        while bonds_changed:
            rresolve = RadicalResolver(start_atom)
            bonds_changed = rresolve().bonds_changed

    return pbmol