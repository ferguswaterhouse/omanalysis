import argparse
import subprocess
from dataclasses import dataclass
import os
import numpy as np
import files
import gromacs as gmx


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

# ====== CALCULATIONS =======
def vec3_magnitude(vec3):
    return (vec3[0] ** 2 + vec3[1] ** 2 + vec3[2] ** 2) ** (1 / 2)


def calculate_angle_to_normal(frame, resn, i, j, normal):
    vecij = np.array((
        frame.residues[resn][j].x - frame.residues[resn][i].x,
        frame.residues[resn][j].y - frame.residues[resn][i].y,
        frame.residues[resn][j].z - frame.residues[resn][i].z
    ))
    vec_normal = np.array((normal[0], normal[1], normal[2]))
    vecij_mag = vec3_magnitude(vecij)
    uvecij = vecij / vecij_mag
    return np.arccos(np.dot(vec_normal, uvecij))


def order_function(angle):
    return (3 / 2) * np.cos(angle) ** 2 - (1 / 2)

# ====== INPUT =======
def read_all_frames(frame_dir, file_name):
    frames = []
    for i in range(len(os.listdir(frame_dir))-1):
        file = '{dir}/{fnm}_{n}.gro'.format(dir=frame_dir, fnm=file_name, n=str(i))
        frames.append(files.read_gro_file(file))
    return frames


def get_all_bonds(chain_topology):
    bonds = []
    for beads in chain_topology.values():
        for i in range(len(beads)-1):
            bonds.append((beads[i], beads[i+1]))
    return bonds

# ====== PROCESSING =======
def calculate_order_for_all_bonds(frame, normal, topology):
    all_bond_orders_for_all_frames = []
    for i, j in get_all_bonds(topology):
        bond_orders_for_frame = []
        for resn in frame.residues.keys():
            angle = calculate_angle_to_normal(frame, resn, i, j, normal)
            bond_orders_for_frame.append(order_function(angle))
        all_bond_orders_for_all_frames.append(np.array(bond_orders_for_frame))
    return np.array(all_bond_orders_for_all_frames)


# def calculate_bond_orders_for_all_frames(frames, normal, topology):
#     all_orders = []
#     for frame in frames:
#         frame_orders = calculate_order_for_all_bonds(frame, normal, topology)
#         all_orders.append(frame_orders)
#     return np.array(all_orders)


# def average_order_data(orders, chain_topologies):
#     bond_traj_averages = {}
#     chain_traj_averages = {}
#     average_per_frame = np.array([np.mean(frame, axis=1) for frame in orders])
#     average_over_trajectory = np.mean(average_per_frame, axis=0)
#     for i, bond in enumerate(get_all_bonds(chain_topologies)):
#         bond_traj_averages[bond] = average_over_trajectory[i]
#     for i, chain in chain_topologies.items():
#         chain_orders = []
#         chain_contents = chain.split()
#         for x in range(len(chain_contents)-1):
#             bond = (chain_contents[x], chain_contents[x+1])
#             chain_orders.append(bond_traj_averages[bond])
#         chain_traj_averages[i] = np.mean(chain_orders)
#     return bond_traj_averages, chain_traj_averages, np.mean(list(chain_traj_averages.values()))
    
# ====== RUNNING =======
def run(xtc_file, tpr_file, lipid, initial_time, final_time, normal, frame_dump_dir):

    print(' > RUNNING LIPID ORDER PARAMETER ANALYSIS...')

    topology = CHAIN_TOPOLOGY[lipid]  
    frame_name = 'frame_dump'

    gmx.frame_dump(xtc_file, tpr_file, initial_time, final_time, 20, lipid, frame_dump_dir, frame_name)
    frames = read_all_frames(frame_dump_dir, frame_name)

    last_frame = len(frames)

    all_orders = []
    for i, frame in enumerate(frames):
        print(' > {}'.format(str(i)) + '/' + str(last_frame), end='\r', flush=True) # Show progress as current frame
        frame_orders = calculate_order_for_all_bonds(frame, normal, topology)
        all_orders.append(frame_orders)
    orders = np.array(all_orders)

    print(' > LIPID ORDER PARAMETER ANALYSIS COMPLETE')
    return orders

    # bond_averages, chain_averages, average = average_order_data(orders, topology)
    # return bond_averages, chain_averages, average
