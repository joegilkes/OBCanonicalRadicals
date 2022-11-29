## NOTE: This package is still under active development. Expect things to change and become significantly more user friendly soon.

# OBCanonicalRadicals

OpenBabel takes a very conservative approach to identifying bonding structure from a 3D geometry (such as an XYZ file) when free radicals are present, often resulting in molecules with very high radical character that, to a trained chemist, should spontaneously collapse to a more stable structure. Finding such a stable structure is simple - just increase the bond order every time two radical atoms are adjacent - but this approach can be inconsistent depending on factors such as the starting point for radical resolution and the direction taken when iterating through the molecule. Without careful consideration of bonding patterns, resulting structures can also end up less stable than their highly radical counterparts.

This package aims to remedy this issue by tidying up OpenBabel's bonding detection in the presence of radicals in a way that is both reproducible (returning a 'canonical' structure for every given radical species) and stable. This allows OpenBabel-detected radical species to be utilised in applications like reaction networks, where obviously unstable radical species can be resolved to their more stable counterparts in a reproducible way in order to cut down the number of species in the network.

## Requirements

* Python >= 3.7.x
* OpenBabel >= 3.1.x

Note that this package WILL NOT work with OpenBabel 2!

## Installation

Due to OpenBabel 3 being quite difficult to compile, most people will have installed OpenBabel through a package distributor like `conda`. If this is the case, installing OBCanonicalRadicals through `pip` will not detect `conda`'s version of OpenBabel, and will then fail to install its own since a compiled version of OpenBabel will not be available on the system. Users should therefore install with `pip install --no-deps ./OBCanonicalRadicals` to utilise `conda`'s version of OpenBabel.

In the event that users do have a locally compiled installation of OpenBabel, running `pip install ./OBCanonicalRadicals` should do the trick.

## Usage

The package is implemeted around running a single function, `obcr.fix_radicals()`. Given a `pybel.Molecule` object, this function will detect the current radical/bonding state of the given molecule, and attempt to rectify the bonding structure by joining adjacent free radicals into bonds wherever it is sensible to do so.

An example of running `fix_radicals()` can be seen in the `run.py` script, which can be used to do canonical conversion of radical SMILES strings or XYZ-based geometries on the command line. Of note is the usage of the `is_radical()` function, which is used to reduce unnecessary computation by only calling `fix_radicals()` on molecules where this may actually be needed.

## Theory

Canonically resolving the radical/bonding structure of an OpenBabel `OBMol` (via Pybel's `Molecule` interface) proceeds through two major steps:

### 1) Determination of canonical starting atom

Resolving the radical structure requires recursively iterating through atoms/bonds within a given molecule, but the resulting structure can change dramatically based on which atom this iteration is chosen to start from. It is therefore important for the same atom to always be chosen as the starting point for a given molecule if we desire a canonical resolution.

In order to determine this starting point, the initial radical state of each atom (according to OpenBabel) is determined, and the atom with the highest radical state is chosen as the starting point. This is chosen based on the theory that atoms with higher radical character will be the least stable atoms in the molecule, therefore they should be resolved preferentially before those with less radical character.

In the event that multiple atoms share the highest radical state (e.g. two diradical atoms exist in the same molecule), the canonical starting point is chosen by finding the atom with the highest neighbouring hydrogenation. This proceeds through recursive searches of neighbour atoms with increasing depth until the correct atom is found. If the maximum depth is reached and the hydrogenation around each potential starting point is still equivalent, it is assumed that there is some form of symmetry in the molecule and that both starting points are equivalent.

#### *Note on hydrogenation metric*

It may be noted that choosing via greatest neighbouring hydrogenation is a mostly arbitrary measure and does not necessarily have much chemical relevance. We believe this to still be a valid metric because there is really no way to know which starting point will lead to the better overall final structure, other than by propagating out the bonding changes from all equivalent starting atoms in all directions, which would come with a greater computational cost. This may be investigated and added as an optional feature in the future. For now, neighbouring hydrogenation is used since it will always select the same starting point for a given molecular structure.

### 2) Propagation of radicals

Once the canonical starting atom has been selected, we iterate down the chain of neighbouring atoms, resolving neighbouring radicals into bonds wherever it is sensible to do so. Specific checks are in place to avoid invalid bonding in rings, and to preferentially bond in patterns that form conjugation rather than forming less stable allenes.

In the event that a starting atom has more than one valid initial direction of radical resolution, we propagate the radicals along all initial directions and determine the best resolution from these candidate final structures.

### 2.5) Determining best resolution

In order to determine which radical resolution should be considered the 'best', multiple checks are performed. Firstly, we check which candidate has the highest overall bond order (as this will implicitly correlate to the remaining number of free radical electrons and thus to stability as a whole) and if there is a winner, this candidate is selected. If overall bond orders are the same, we perform a more extensive check to select the candidate with the least atoms with high radical state (e.g. a candidate with 3 single radical atoms will be chosen over a candidate with a diradical and a radical).

## Examples

in progress