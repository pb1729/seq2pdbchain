import numpy as np
from collections import defaultdict

from pdb_util import parse_pdb, write_xyz


# This file is meant for code generation.
# Given a particular pdb file, we pick a representative of each
# kind of amino acid from the file. We then generate python code
# for a dict that contains these structures.


PDB_PATH = "7sbp.pdb" # can obtain this file here: https://files.rcsb.org/download/7SBP.pdb
DATA_DIR = None
representative_list = [
    "C972", "B509", "A1135", "A53",  "A136",
    "A239", "A191", "A485",  "B49",  "A101",
    "A390", "B129", "A731",  "A168", "A217",
    "A50",  "A376", "B64",   "A170", "A401"
]


def tounit(v):
    return v / np.linalg.norm(v)

def normalize_positions(blk):
    """ standardize positions and orientations of atoms in a residue """
    # translate:
    for atom in blk:
        if atom.name == "N":
            origin = np.array(atom.pos)
    blk = [atom._replace(pos = np.array(atom.pos) - origin) for atom in blk]
    # rotate:
    v_n, v_CA = None, None
    for atom in blk:
        if atom.name == "N_next":
            v_n = atom.pos
        elif atom.name == "CA":
            v_CA = atom.pos
    if v_n is None or v_CA is None: return False
    rot_matrix = np.zeros((3,3))
    rot_matrix[0] = tounit(v_n)
    u = v_CA - v_n*(v_CA @ v_n)/(v_n @ v_n) # component of v_CA orthogonal to V_C
    rot_matrix[1] = tounit(u)
    rot_matrix[2] = np.cross(rot_matrix[0], rot_matrix[1]) # remaining orthonormal vector, right handed
    return [atom._replace(pos = rot_matrix @ atom.pos) for atom in blk]

def get_alpha_rotation(blk):
    """ get the rotation angle between one alpha bond and the next.
        requires that block be an output of normalize_positions. """
    v_N, v_CA, v_Nn, V_CAn = None, None, None, None
    for atom in blk:
        if atom.name == "N":
            v_N = atom.pos
        elif atom.name == "CA":
            v_CA = atom.pos
        elif atom.name == "N_next":
            v_Nn = atom.pos
        elif atom.name == "CA_next":
            v_CAn = atom.pos
    if v_N is None or v_CA is None or v_Nn is None or v_CAn is None:
        return False
    v1 = (v_CA  - v_N )[1:] # project into y-z plane
    v2 = (v_CAn - v_Nn)[1:] # project into y-z plane
    return np.arccos(tounit(v1) @ tounit(v2))


if __name__ == "__main__":
    records = parse_pdb(PDB_PATH)
    atoms = [record for record in records if record != None]
    residue_blks = defaultdict(lambda: [])
    # add regular atoms
    for atom in atoms:
        residue_blks[atom.residue_id].append(atom)
    # append first N atom and CA atom of next residue
    for atom in atoms:
        if atom.name == "N" or atom.name == "CA":
            chain, seq = atom.residue_id[0], int(atom.residue_id[1:])
            prev_residue_id = chain + str(seq - 1)
            if prev_residue_id in residue_blks:
                newatom = atom._replace(residue_id = prev_residue_id, name = atom.name + "_next")
                residue_blks[prev_residue_id].append(newatom)
    for residue_id in representative_list:
        blk = normalize_positions(residue_blks[residue_id])
        if blk == False: print("failed!", residue_id)
        else:
            if DATA_DIR is not None:
              fnm = "%s/%s_%s.xyz" % (DATA_DIR, blk[0].residue, residue_id)
              write_xyz(fnm, blk)
            print("%s: %d\u00b0" % (blk[0].residue, int(get_alpha_rotation(blk) * 180/3.1416)))
    # wacky code generation part:
    print("\n\t--- Residue Construction Info: ---\n")
    print("{")
    for residue_id in representative_list:
        blk = normalize_positions(residue_blks[residue_id])
        print("    \"%s\": [" % blk[0].residue)
        for atom in blk:
            if atom.name != "CA_next":
                print("        Atom(\"%s\", np.%s, \"%s\")," % (atom.name, repr(atom.pos), atom.element))
        print("    ],")
    print("}")
    print()




