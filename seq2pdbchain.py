import sys
import numpy as np

from amino_data import letter_code, structures, Atom
from pdb_util import pdb_line
from structures_codegen import extract_atom_positions


# Maximum distance from any atom in an amino to the peptide nitrogen
MAX_DIST_TO_N = 10.0              # [Å] rounded up from 8.48
ATOM_RADIUS = 0.68                # [Å] atom radius for C, N, O
COLLISION_DIST = 3.*ATOM_RADIUS  # [Å] twice atom radius to form a bond, plus 1.0 safety margin
COMPACTNESS_ITERATIONS = 200      # number of iterations to spend optimizing compactness

def BETA(n_collisions): # [1/collisions] thermodynamic coolness for Monte Carlo method
    return max(0.2, 10/n_collisions)


def blks_to_pdb(residues, blks):
    lines = []
    i_atom = 0
    for i_seq, (residue, blk) in enumerate(zip(residues, blks)):
        for atom in blk:
            lines.append(pdb_line(atom, i_atom, 1 + i_seq, residue))
            i_atom += 1
    lines.append("TER " + pdb_line(Atom(" ", [0., 0., 0.], " "), i_atom, len(residues), residue)[4:27])
    return "\n".join(lines)


def amino_choices_to_blks(structs):
    blks = []
    cursor = np.array([0., 0., 0.])
    transform = np.identity(3)
    for struct in structs:
        blk = []
        for atom in struct.atoms:
            transformed_atom = atom._replace(pos = cursor + (transform @ atom.pos))
            if atom.name != "N_next":
                blk.append(transformed_atom)
            else:
                cursor = transformed_atom.pos
                break # NOTE: N_next should be the last element in the list!
        transform = transform @ struct.transform
        blks.append(blk)
    # do a proper OC1, OC2 termination of the peptide
    (CA_pos, C_pos, O_pos) = extract_atom_positions(["CA", "C", "O"], blks[-1])
    OC2_pos = 3*C_pos - (CA_pos + O_pos)
    idx_O = [i for i in range(len(blks[-1])) if blks[-1][i].name == "O"][0]
    blks[-1][idx_O] = blks[-1][idx_O]._replace(name = "OC1")
    blks[-1].append(Atom("OC2", OC2_pos, "O"))
    return blks


def count_collisions(blks):
    blks = [ # convert each block into an (n, 3) numpy array for faster processing
        np.stack([atom.pos for atom in blk])
        for blk in blks]
    collisions = 0
    for i in range(len(blks)):
        for j in range(i - 1): # directly adjacent aminos are assumed to be okay
            blk_i, blk_j = blks[i][None, :], blks[j][:, None]
            if np.linalg.norm(blk_i[0, 0] - blk_j[0, 0]) > 2*MAX_DIST_TO_N:
                continue # guaranteed no collisions thanks to triangle inequality
            collisions += (np.linalg.norm(blk_i - blk_j, axis=2) < COLLISION_DIST).sum()
    return collisions


def rough_radius(blks):
    pos_sum = np.array([0., 0., 0.])
    n = 0
    for blk in blks:
        for atom in blk:
            pos_sum += atom.pos
            n += 1
    pos_avg = pos_sum / n
    return max([np.linalg.norm(atom.pos - pos_avg)
        for blk in blks for atom in blk])


def pick_random_struct(residue):
    choices = structures[residue]
    return choices[np.random.randint(len(choices))]


def pdb_chain(sequence, showprogress=True):
    structs = []
    residues = []
    for c in sequence:
        res = letter_code[c]
        residues.append(res)
        structs.append(pick_random_struct(res))
    def mutate(struct_list):
        ans = [struct for struct in struct_list]
        j = np.random.randint(len(residues) - 2)
        ans[j] = pick_random_struct(residues[j])
        ans[j + 1] = pick_random_struct(residues[j + 1])
        ans[j + 2] = pick_random_struct(residues[j + 2])
        return ans
    n_collisions = 1000000000000
    while n_collisions > 0:
        new_structs = mutate(structs)
        blks = amino_choices_to_blks(new_structs)
        new_n_collisions = count_collisions(blks)
        # markov chain monte-carlo step:
        if np.exp(BETA(n_collisions)*(n_collisions - new_n_collisions)) > np.random.rand():
            n_collisions = new_n_collisions
            structs = new_structs
            if showprogress:
                sys.stderr.write("current collisions: %d\n" % n_collisions)
    radius = rough_radius(blks)
    for i in range(COMPACTNESS_ITERATIONS):
        new_structs = mutate(structs)
        blks = amino_choices_to_blks(new_structs)
        new_radius = rough_radius(blks)
        if new_radius <= radius and count_collisions(blks) == 0:
            radius = new_radius
            structs = new_structs
            if showprogress:
                sys.stderr.write("current radius: %f\n" % radius)
    return blks_to_pdb(residues, blks)


if __name__ == "__main__":
    print(pdb_chain(input()))
    




