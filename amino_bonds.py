from collections import defaultdict
import numpy as np

from .amino_data import structures, letter_code

#
# The purpose of this file is to determine (from atom position data) where the bonds are in each amino acid.
#

# Covalent Radii from https://en.wikipedia.org/wiki/Van_der_Waals_radius
ELEMENT_TO_VDW_RADII = {
  "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8
}

def _get_atom_rad(atom):
  return ELEMENT_TO_VDW_RADII[atom.element]

def _residue_example_to_bond_graph(multiset, residue_example):
  for i, atom_i in enumerate(residue_example.atoms):
    if atom_i.name == "N_next": continue
    for j in range(i):
      atom_j = residue_example.atoms[j]
      dist = ((atom_i.pos - atom_j.pos)**2).sum()
      if dist < _get_atom_rad(atom_i) + _get_atom_rad(atom_j):
        multiset[j, i] += 1

def _compute_bond_graphs():
  ans = {}
  for letter in letter_code:
    bond_counts = defaultdict(int)
    n_examples = 0
    for residue_example in structures[letter_code[letter]]:
      _residue_example_to_bond_graph(bond_counts, residue_example)
      n_examples += 1
    # consistency check
    for key in bond_counts:
      assert bond_counts[key] == n_examples
    # update the dict
    src, dst = [], []
    for i, j in bond_counts:
      src.append(i)
      dst.append(j)
      src.append(j)
      dst.append(i)
    ans[letter] = (np.array(src, dtype=np.int32), np.array(dst, dtype=np.int32))
  return ans


bond_graphs = _compute_bond_graphs()


if __name__ == "__main__":
  import atoms_display
  print("Visualization test for amino bonds. Pick an Amino:")
  print(" ".join(letter for letter in letter_code))
  letter = input("> ")
  residue_example = structures[letter_code[letter]][0]
  # use [:-1] on next two lines to ditch N_next
  positions = np.stack([atom.pos for atom in residue_example.atoms])[:-1]
  atomic_numbers = np.array([atoms_display.ATOMIC_NUMBERS[atom.element] for atom in residue_example.atoms])[:-1]
  poslist = [positions]
  anlist = [atomic_numbers]
  srcs, dsts = bond_graphs[letter]
  for i, j in zip(srcs, dsts):
    if i < j:
      t = np.linspace(0.3, 0.7, 10)[:, None]
      poslist.append(t*positions[i] + (1 - t)*positions[j])
      anlist.append(np.ones(10, dtype=int))
  positions = np.concatenate(poslist)
  atomic_numbers = np.concatenate(anlist)
  display = atoms_display.launch_atom_display(atomic_numbers, positions, radii_scale=0.7)
