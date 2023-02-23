import MDAnalysis as mda
import numpy as np

def run(gro_file, xtc_file, innerleaflet_lipid, outerleaflet_lipid, innerleaflet_bead, outerleaflet_bead, frames):

    u = mda.Universe(gro_file, xtc_file)

    print(' > RUNNING NEAREST NEIGHBOUR ANALYSIS...')

    innerleaflet_beads = u.select_atoms('resname {} and name {}'.format(innerleaflet_lipid, innerleaflet_bead))
    outerleaflet_beads = u.select_atoms('resname {} and name {}'.format(outerleaflet_lipid, outerleaflet_bead))

    membrane_thicknesses = []
    for timestep in u.trajectory[:frames]:
        print(' > {}'.format(str(timestep.frame)) + '/' + str(frames), end='\r', flush=True) # Show progress as current frame
        innerleaflet_boundary = innerleaflet_beads.centroid()[2]
        outerleaflet_boundary = outerleaflet_beads.centroid()[2]
        membrane_thickness = outerleaflet_boundary - innerleaflet_boundary
        membrane_thicknesses.append(membrane_thickness)
    
    print(' > COMPLETE')
    print(' > MEMBRANE THICKNESS =', np.mean(membrane_thicknesses), '+/-', np.std(membrane_thicknesses))

    return np.array(membrane_thicknesses)
