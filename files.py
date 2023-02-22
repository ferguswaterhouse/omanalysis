

def read_topology(top_dir):
    molecules = {}
    contents = []
    with open(top_dir, 'r') as f:
        for line in f.readlines():
            contents.append(line.split())
    contents = contents[contents.index(['[', 'molecules', ']'])+2:] # discards top information
    for content in contents:
        if content[0] not in molecules.keys():
            molecules[content[0]] = int(content[1])
        else:
            molecules[content[0]] += int(content[1])
    return molecules
