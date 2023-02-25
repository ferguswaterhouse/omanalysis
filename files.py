from dataclasses import dataclass


def read_topology(top_file):
    molecules = {}
    contents = []
    with open(top_file, 'r') as f:
        for line in f.readlines():
            contents.append(line.split())
    contents = contents[contents.index(['[', 'molecules', ']'])+2:] # discards top information
    for content in contents:
        if content[0] not in molecules.keys():
            molecules[content[0]] = int(content[1])
        else:
            molecules[content[0]] += int(content[1])
    return molecules


def read_index(ndx_file):
    group_indices = []
    with open(ndx_file) as f:
        for line in f.readlines():
            line_data = line.split()
            if line_data[0] == '[':
                group_indices.append(line_data[1])
    return group_indices

# READING .GRO FILES
@dataclass
class Frame:
    step: int
    time: float
    dimensions: any
    residues: any # Dict of Dicts: {Resn : {Beadnm : Bead}}


@dataclass
class Bead:
    name: str
    numb: int
    x: float
    y: float
    z: float


def read_gro_file(file):
    residues = {}
    with open(file, 'r') as f:
        lines = f.readlines()
        time, step = float(lines[0].split()[3]), int(lines[0].split()[5])
        no_beads = int(lines[1].replace(' ', ''))
        dimensions = (
            float(lines[-1].split()[0]),
            float(lines[-1].split()[1]),
            float(lines[-1].split()[2])
        )
        for line in lines[2:-1]:
            resn = int(line[0:5].replace(' ', ''))
            # resnm = line[5:9].replace(' ', '')
            beadnm = line[12:15]
            beadn = int(line[15:20].replace(' ', ''))
            beadx = float(line[20:28].replace(' ', ''))
            beady = float(line[28:36].replace(' ', ''))
            beadz = float(line[36:].replace(' ', '').replace('\n', ''))
            bead = Bead(
                name=beadnm,
                numb=beadn,
                x=beadx,
                y=beady,
                z=beadz
            )
            if resn not in residues.keys():
                residues[resn] = {beadnm: bead}
            else:
                residues[resn][beadnm] = bead
    return Frame(
        step=step,
        time=time,
        dimensions=dimensions,
        residues=residues
    )


def read_gro(gro):
    contents = []
    with open(gro, 'r') as f:
        for line in f.readlines()[2:-1]:
            content = line.split()
            content[0] = content[0].replace('1', '')
            if len(content[1]) > 3:
                name = content[1][:3]
                number = content[1][3:]
                content[1:2] = (name, number)
            contents.append(content)
    return contents

# DIVIDING .GRO FILES INTO MOLECULES
def divide_contents_into_molecule(contents):
    molecules = {}
    for content in contents:
        if content[0] not in molecules.keys():
            molecules[content[0]] = [content]
        else:
            molecules[content[0]].append(content)
    return molecules


def split_contents_into_residues(molecule_contents):
    residues = []

    first_atom_name = molecule_contents[0][1]
    molecule_length = 1
    while molecule_contents[molecule_length][1] != first_atom_name:
        molecule_length += 1
    
    residue_number, residue_length = 0, 0
    for content in molecule_contents:
        if residue_length >= molecule_length:
            residue_number += 1
            residue_length = 0
        residues.append(residue_number + 1)
        residue_length += 1

    return residues


def write_gro_file(gro_file, residues, molecule, out):
    with open(gro_file, 'r') as f:
        lines = f.readlines()

    with open(out, 'w') as w:
        number_of_molecules_parsed = 0
        for line in lines:
            if line.split()[0].replace('1', '') == molecule:
                new_line = (5-len(str(residues[number_of_molecules_parsed]))) * ' ' +  str(residues[number_of_molecules_parsed]) + line[5:]
                number_of_molecules_parsed += 1
            else:
                new_line = line
            w.write(new_line)


def seperate_residues(gro_file, molecule, out_file):
    contents = read_gro(gro_file)
    molecules = divide_contents_into_molecule(contents)
    residues = split_contents_into_residues(molecules[molecule])
    write_gro_file(gro_file, residues, molecule, out_file)