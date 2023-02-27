import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import argparse


def get_unique_bead_types(bead_group):
    unique_bead_types = []
    for bead_type in bead_group.names:
        if bead_type not in unique_bead_types:
            unique_bead_types.append(bead_type)
    return unique_bead_types


def convert_to_array_in_terms_of_bead_types(bead_types, all_bead_contacts_per_frame):

    all_ion_contacts_in_traj = [] # Array: Frame x Bead Type x Residue

    for ion_contacts in all_bead_contacts_per_frame:
        all_lps_ion_contacts_in_frame = [] # Array: Bead Type x Residue

        for bead_type in bead_types:
            all_lps_ion_contacts_for_bead_type_in_frame = ion_contacts[np.where(bead_types.names == bead_type)] # Array: Residue

            all_lps_ion_contacts_in_frame.append(all_lps_ion_contacts_for_bead_type_in_frame)
        all_ion_contacts_in_traj.append(all_lps_ion_contacts_in_frame)

    return np.array(all_ion_contacts_in_traj)


def run(gro_dir, xtc_dir, ion, lps, frames, radius):

    print(f' > RUNNING CONTACT ANALYSIS BETWEEN {lps} AND {ion}...')

    # INPUT

    bead_selections = {
        'REMP': "resname REMP and name GM* P* S*",
        'RAMP': "resname RAMP and name GM* P* S*"
    }

    u = mda.Universe(gro_dir, xtc_dir)

    lps_bead_group = u.select_atoms(bead_selections[lps])
    ion_bead_group = u.select_atoms(f'name {ion}')

    number_of_lps_beads = len(lps_bead_group)

    lps_bead_types = get_unique_bead_types(lps_bead_group) # Gets a list of all the unique bead types in the lps

    lps_ion_contacts_per_frame = np.zeros((frames, number_of_lps_beads), dtype=int) # Defines a matrix to store the number of ion contacts for each bead at each frame
    # Array: Frame x Bead

    # CALCULATIONS
    for ts in u.trajectory[:frames]:

        print(' > {}'.format(str(ts.frame)) + '/' + str(frames), end='\r', flush=True) # Show progress as current frame

        lps_ion_distance = contacts.distance_array(lps_bead_group.positions, ion_bead_group.positions) # Calculates distance between each lps bead and each ion
        lps_ion_contacts = contacts.contact_matrix(lps_ion_distance, radius) # Returns a boolean matrix representing: True = ion less than radius away from bead

        # Sums all ions contacts for each bead
        for bead_i in range(number_of_lps_beads):
            number_of_ion_contacts_for_bead_i = lps_ion_contacts[bead_i].sum()
            lps_ion_contacts_per_frame[ts.frame, bead_i] = number_of_ion_contacts_for_bead_i

    bead_type_contacts = convert_to_array_in_terms_of_bead_types(lps_bead_types, lps_ion_contacts_per_frame) # Array: Frame x Bead Type x Residue = # Ion Contacts

    # PRINT RESULTS
    avg_bead_type_contacts_per_residue = np.mean(bead_type_contacts, axis=2)
    avg_bead_type_contacts_per_residue_per_frame = np.mean(avg_bead_type_contacts_per_residue, axis=0)
    
    print(' > COMPLETE.')
    print(' > RESULTS:')
    for i in range(len(lps_bead_types)):
        print(' {} = {}'.format(lps_bead_types[i], round(avg_bead_type_contacts_per_residue_per_frame[i], 3)))

    return bead_type_contacts # Array: Frame x Bead Type x Residue = # Ion Contacts


if __name__ == '__main__':
    default_contact_radius = 5.0
    parser = argparse.ArgumentParser()
    parser.add_argument('-gro', type=str, help='Path to the .gro file')
    parser.add_argument('-xtc', type=str, help='Path to the .xtc file')
    parser.add_argument('-ion', type=str, help='Ion name')
    parser.add_argument('-lps', type=str, help='LPS type')
    parser.add_argument('-frames', type=int, help='Number of frames to analyse')
    parser.add_argument('-rad', type=float, help='Max contact distance', default=default_contact_radius)
    args = parser.parse_args()
    run(args.gro, args.xtc, args.ion, args.lps, args.frames, args.rad)
