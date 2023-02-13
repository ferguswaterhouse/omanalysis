# Takes molecule and bead name to treat as molecule centre and calculates the nearest neighbour index for each molecule in the trajectory.

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import numpy as np
import argparse


def nearest_neighbour_calculation(unique_neighbour_set, number_of_neighbours, total_molecules):
    number_of_unique_neighbours_over_traj = len(unique_neighbour_set) - number_of_neighbours
    all_possible_neighbours = total_molecules - number_of_neighbours - 1
    return number_of_unique_neighbours_over_traj / all_possible_neighbours


def run(gro_file, xtc_file, molecule, frames, number_of_neighbours, bead_name='GM1'):

    u = mda.Universe(gro_file, xtc_file)

    print(' > RUNNING NEAREST NEIGHBOUR ANALYSIS...')

    molecules = u.select_atoms('resname {} and name {}'.format(molecule, bead_name))
    total_molecules = len(molecules)

    closest_neighbours_over_traj = [set() for i in range(total_molecules)]

    for timestep in u.trajectory[:frames]:
        print(' > {}'.format(str(timestep.frame)) + '/' + str(frames), end='\r', flush=True) # Show progress as current frame
        molecule_positions = molecules.positions
        for i in range(total_molecules):
            distances_between_all_molecules = dist.distance_array(reference=molecule_positions[i], configuration=molecule_positions, box=timestep.dimensions)
            closest_neighbours_including_self = np.argpartition(distances_between_all_molecules[0], number_of_neighbours + 1)[:number_of_neighbours + 1]
            closest_neighbours = np.delete(closest_neighbours_including_self, np.where(closest_neighbours_including_self == i))
            closest_neighbours_over_traj[i].update(closest_neighbours)

    nearest_neighbour_indices = [nearest_neighbour_calculation(unique_neighbour_set, number_of_neighbours, total_molecules) for unique_neighbour_set in closest_neighbours_over_traj]

    avg_nn_index, std_nn_index = np.mean(nearest_neighbour_indices), np.std(nearest_neighbour_indices)
    
    print(' > NN Index:', avg_nn_index, '+/-', std_nn_index)
    print(' > COMPLETE')
    
    return np.array(nearest_neighbour_indices)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gro')
    parser.add_argument('-xtc')
    parser.add_argument('-mol')
    parser.add_argument('-nn')
    parser.add_argument('-frames')
    args = parser.parse_args()
    run(args.gro, args.xtc, args.mol, int(args.frames), int(args.nn))
