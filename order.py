import MDAnalysis as mda
import numpy as np


CHAIN_TOPOLOGY = {
    "REMP": {
        1: ['GL2', 'C1A', 'C2A', 'C3A'],
        2: ['GL1', 'C1B', 'C2B', 'C3B'],
        3: ['GL3', 'C1D', 'C2D', 'C3D'],
        4: ['GL4', 'C1C', 'C2C', 'C3C'],
        5: ['GL6', 'C1E', 'C2E'],
        6: ['GL8', 'C1F', 'C2F']
    },
    "RAMP": {
        1: ['GL2', 'C1A', 'C2A', 'C3A'],
        2: ['GL1', 'C1B', 'C2B', 'C3B'],
        3: ['GL3', 'C1D', 'C2D', 'C3D'],
        4: ['GL4', 'C1C', 'C2C', 'C3C'],
        5: ['GL6', 'C1E', 'C2E'],
        6: ['GL8', 'C1F', 'C2F']
    }
}


def get_bonds(chain):
    bonds = []
    for i in range(len(chain)-1):
        bonds.append((chain[i], chain[i+1]))
    return bonds


def vec3_magnitude(vec3):
    return (vec3[0] ** 2 + vec3[1] ** 2 + vec3[2] ** 2) ** (1 / 2)


def calculate_angle_to_normal(vectorij, normal):
    vec_normal = np.array((normal[0], normal[1], normal[2]))
    vectorij_mag = vec3_magnitude(vectorij)
    uniform_vectorij = vectorij / vectorij_mag
    return np.arccos(np.dot(vec_normal, uniform_vectorij))


def order_function(angle):
    return (3 / 2) * np.cos(angle) ** 2 - (1 / 2)


def calculate_bond_order(vectorij, normal):
    angle = calculate_angle_to_normal(vectorij, normal)
    return order_function(angle)


def run(gro_file, xtc_file, molecule, normal):

    print(' > LIPID ORDER PARAMETER ANALYSIS OF {m}...'.format(m=molecule))

    # LOADS INPUT AND SELECTS LIPID TAIL BEADS  
    u = mda.Universe(gro_file, xtc_file) # NEEDS NO JUMP XTC FILE!

    lipid_tail_bead_names = [bead for chain in CHAIN_TOPOLOGY[molecule].values() for bead in chain]
    lipid_tail_beads = u.select_atoms('resname {} and name {}'.format(molecule, ' '.join(lipid_tail_bead_names)))

    # FLATTENS CHAIN TOPOLOGY INTO LIST OF BONDS: [(bead_name_i, bead_name_j), ...]
    list_of_all_bond_types = [(bead_name_i, bead_name_j) for chain in CHAIN_TOPOLOGY[molecule].values() for bead_name_i, bead_name_j in get_bonds(chain)]
    
    # FINDS CORRESPONDING ATOMS IN UNIVERSE FOR ALL BONDS: [[(atom_i, atom_j), ...], ...]
    corresponding_atoms_for_all_bond_types = [[] for i in range(len(list_of_all_bond_types))]
    for residue in lipid_tail_beads.residues:
        for i, (bead_name_i, bead_name_j) in enumerate(list_of_all_bond_types):
            atom_i = residue.atoms.select_atoms('name {}'.format(bead_name_i))[0]
            atom_j = residue.atoms.select_atoms('name {}'.format(bead_name_j))[0]
            corresponding_atoms_for_all_bond_types[i].append((atom_i, atom_j))
    
    # ITERATES OVER ALL FRAMES OF THE TRAJECTORY AND CALCULATES BOND ORDER FOR ALL BONDS
    last_frame = len(u.trajectory)
    all_bond_orders_list = [[] for i in range(len(list_of_all_bond_types))]
    # Iterates over all frames of the trajectory
    for timestep in u.trajectory:
        print(' > {cf}/{lf}  '.format(cf=str(timestep.frame), lf=str(last_frame)), end='\r', flush=True)
        # Iterates over all bond types
        for bond_i, bonds in enumerate(corresponding_atoms_for_all_bond_types):
            bond_orders = []
            # Iterates over all residues
            for atom_i, atom_j in bonds:
                position_i, position_j = atom_i.position, atom_j.position
                vectorij = position_j - position_i
                bond_orders.append(calculate_bond_order(vectorij, normal))
            all_bond_orders_list[bond_i].extend(bond_orders)
    
    all_bond_orders_array = np.array(all_bond_orders_list)
    
    # PRINTS RESULTS
    list_of_corresponding_chains = [chain_i for chain_i, chain in CHAIN_TOPOLOGY[molecule].items() for i in range(len(get_bonds(chain)))]

    avg_bond_orders = np.mean(all_bond_orders_array, axis=1)
    print(' > LIPID ORDER PARAMETERS ANALYSIS COMPLETE')
    print(' > RESULTS:')

    print(' '.join([bead_name_i + '-' + bead_name_j for bead_name_i, bead_name_j in list_of_all_bond_types]))
    print('   ' + '       '.join(str(chain_i) for chain_i in list_of_corresponding_chains))
    print('  ' + '   '.join([str(round(order, 3)) for order in avg_bond_orders]))

    return all_bond_orders_array, list_of_all_bond_types, list_of_corresponding_chains
