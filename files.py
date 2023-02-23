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
