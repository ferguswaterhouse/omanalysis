

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
