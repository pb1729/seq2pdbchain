import numpy as np

from amino_data import letter_code, structures



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
        struct = structures[residue][np.random.randint(len(structures[residue]))]
        for atom in struct.atoms:
            transformed_atom = atom._replace(pos = cursor + (transform @ atom.pos))
            if atom.name != "N_next":
                lines.append(pdb_line(transformed_atom, i_atom, i_seq, residue))
                i_atom += 1
            else:
                cursor = transformed_atom.pos
                break # NOTE: N_next should be the last element in the list!
        transform = transform @ struct.transform
    return "\n".join(lines)


if __name__ == "__main__":
    print(pdb_chain(input()))
    




