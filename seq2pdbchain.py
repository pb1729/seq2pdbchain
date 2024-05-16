import numpy as np

from amino_data import letter_code, structures


helix_step_angle = 160 * np.pi / 180



def leftpad(s, n):
    return " "*(n-len(s)) + s

def pdb_name(s):
    return " " + s + " "*(3-len(s))

def pdb_float(f):
    return "%8.3f" % f

def pdb_line(atom, i_atom, i_seq, residue):
    """ based on this: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html """
    return ("ATOM  "                # 1-6
        + leftpad(str(i_atom), 5)   # 7-11
        + " "                       # 12
        + pdb_name(atom.name)       # 13-16
        + " "                       # 17
        + residue                   # 18-20
        + " A"                      # 21-22
        + leftpad(str(i_seq), 4)    # 23-26
        + "    "                    # 27-30
        + pdb_float(atom.pos[0])    # 31-38
        + pdb_float(atom.pos[1])    # 39-46
        + pdb_float(atom.pos[2])    # 47-54
        + "  1.00"                  # 55-60
        + "100.00"                  # 61-66
        + "          "              # 67-76
        + leftpad(atom.element, 2)  # 77-78
        + "  "                      # 79-80
    )



def pdb_chain(sequence):
    lines = []
    cursor = np.array([0., 0., 0.])
    transform = np.identity(3)
    i_atom = 0
    for i_seq, c in enumerate(sequence):
        residue = letter_code[c]
        struct = structures[residue]
        rotate = np.array([
            [1.0,           0.0,            0.0],
            [0.0, np.cos(helix_step_angle), -np.sin(helix_step_angle)],
            [0.0, np.sin(helix_step_angle),  np.cos(helix_step_angle)]])
        new_transform = transform @ rotate
        for atom in struct:
            transformed_atom = atom._replace(pos = cursor + (new_transform @ atom.pos))
            if atom.name != "N_next":
                lines.append(pdb_line(transformed_atom, i_atom, i_seq, residue))
                i_atom += 1
            else: # NOTE: This relies on N_next being the last element in the list!
                cursor = transformed_atom.pos
                break
        transform = new_transform
    return "\n".join(lines)


if __name__ == "__main__":
    print(pdb_chain(input("enter a sequence to be converted to a pdb file > ")))
    




