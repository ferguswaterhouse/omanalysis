import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import numpy as np


def nearest_neighbour_calculation(unique_neighbour_set, number_of_neighbours, total_molecules):
    # NNI = (number of unique neighbours) / (total number of molecules - number of neighbours per frame - 1)
    number_of_unique_neighbours_over_traj = len(unique_neighbour_set) - number_of_neighbours
    all_possible_neighbours = total_molecules - number_of_neighbours - 1
    return number_of_unique_neighbours_over_traj / all_possible_neighbours


def run(gro_file, xtc_file, molecule, frames, number_of_neighbours, bead_name='GM1'):

    print(' > RUNNING NEAREST NEIGHBOUR ANALYSIS...')

    u = mda.Universe(gro_file, xtc_file)

    molecules = u.select_atoms('resname {} and name {}'.format(molecule, bead_name))
    total_molecules = len(molecules)

    closest_neighbours_over_traj = [set() for i in range(total_molecules)]

    for timestep in u.trajectory[:frames]:
        print(' > {}'.format(str(timestep.frame)) + '/' + str(frames), end='\r', flush=True) # Show progress as current frame
        molecule_positions = molecules.positions
        for molecule_nr in range(total_molecules):
            distances_between_all_molecules = dist.distance_array(reference=molecule_positions[molecule_nr], configuration=molecule_positions, box=timestep.dimensions)

            closest_neighbours_including_self = np.argpartition(distances_between_all_molecules[0], number_of_neighbours + 1)[:number_of_neighbours + 1]
            closest_neighbours = np.delete(closest_neighbours_including_self, np.where(closest_neighbours_including_self == molecule_nr)) # Remove self from neighbours

            closest_neighbours_over_traj[molecule_nr].update(closest_neighbours)

    nearest_neighbour_indices = [nearest_neighbour_calculation(unique_neighbour_set, number_of_neighbours, total_molecules) for unique_neighbour_set in closest_neighbours_over_traj]

    avg_nn_index, std_nn_index = np.mean(nearest_neighbour_indices), np.std(nearest_neighbour_indices)
    
    print(' > NN Index =', avg_nn_index, '+/-', std_nn_index)
    print(' > COMPLETE\n')
    
    return np.array(nearest_neighbour_indices) # Array; shape=(number of molecules); values=Nearest Neighbour Index
