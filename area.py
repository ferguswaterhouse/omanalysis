import MDAnalysis as mda
import numpy as np
import files


def run(gro_file, xtc_file, top_file, molecule, frames):

    u = mda.Universe(gro_file, xtc_file)

    print(' > CALCULATING AREA PER LIPID...')

    number_of_molecules_in_leaflet = files.read_topology(top_file)[molecule]

    time = []
    x_dimensions, y_dimensions = [], []
    for timestep in u.trajectory[:frames]:
        print(' > {}'.format(str(timestep.frame)) + '/' + str(frames), end='\r', flush=True) # Show progress as current frame
        time.append(timestep.time * 10 ** -3) # from ps to ns
        x_dimensions.append(timestep.dimensions[0])
        y_dimensions.append(timestep.dimensions[1])

    area = np.array(x_dimensions) * np.array(y_dimensions)
    area_per_lipid = area * (10 ** -2) / number_of_molecules_in_leaflet

    print(' > AREA PER LIPID =', np.mean(area_per_lipid), '+/-', np.std(area_per_lipid), 'nm^2')
    print(' > COMPLETE\n')
    
    return np.array([time, area_per_lipid]) # Array; shape=(timesteps); values=Area per lipid
