import argparse
import subprocess
import gromacs as gmx
import neighbour as nn
import area
import index
import thickness
import torsion
import order
import contacts
import files

# TO DO:
# - Saving all data
# - Add option to run all scripts directly + give help
# - Analysis using units straight
# - Add option to limit frames
# - Processing of data to average repeats
# - Allow analysis of single repeats

def complete(repeats, gro_files, xtc_files, tpr_files, top_files, lps, ions, frames, out_dir):

    print(' > RUNNING COMPLETE ANALYSIS...')

    # ====== CREATE OUTPUT DIRECTORIES ======
    subprocess.run(f'mkdir {out_dir}', shell=True)
    for directory in ['density', 'msd', 'nearest_neighbour', 'area', 'order', 'torsion', 'thickness', 'setup']:
        subprocess.run(f'mkdir {out_dir}/{directory}', shell=True)
    print(' > CREATED OUTPUT DIRECTORIES')

    # ====== SETUP ======
    # File directories
    out_sep_gro_files = [f'sep_r{str(repeat+1)}.gro' for repeat in range(repeats)]
    out_sep_gro_dirs = ['/'.join([out_dir, 'setup', out_sep_gro_files[repeat]]) for repeat in range(repeats)]

    out_ndx_files = [f'complete_index_r{str(repeat+1)}.ndx' for repeat in range(repeats)]
    default_ndx_files = [f'default_index_r{str(repeat+1)}.ndx' for repeat in range(repeats)]
    out_ndx_dirs = ['/'.join([out_dir, 'setup', out_ndx_files[repeat]]) for repeat in range(repeats)]
    default_ndx_dirs = ['/'.join([out_dir, 'setup', default_ndx_files[repeat]]) for repeat in range(repeats)]

    out_no_jump_xtc_files = [f'nojump_r{str(repeat+1)}.xtc' for repeat in range(repeats)]
    out_no_jump_xtc_dirs = ['/'.join([out_dir, 'setup', out_no_jump_xtc_files[repeat]]) for repeat in range(repeats)]

    # Define each LPS residue - this was done as for an unknown reason all lps in output .gro files were given the same residue number
    for repeat in range(repeats):
        files.seperate_residues(gro_files[repeat], lps, out_sep_gro_dirs[repeat])
    print(' > DEFINED LPS RESIDUES FOR EACH .GRO FILE')

    # Create complete index files
    for repeat in range(repeats):
        index.run(out_sep_gro_dirs[repeat], lps, ions, default_ndx_dirs[repeat], out_ndx_dirs[repeat])
    print(' > CREATED COMPLETE INDEX FILES')

    # Create no-jump trajectories
    for repeat in range(repeats):
        gmx.no_jump(xtc_files[repeat], tpr_files[repeat], out_ndx_dirs[repeat], out_no_jump_xtc_dirs[repeat])
    print(' > NO JUMP TRAJECTORIES CREATED')

    # ====== DENSITY ======
    # Group names to calculate density for - corresponding to the index file
    density_group_names = ['W']
    universal_lps_group_names = [
        '{}_HEAD_PO41'.format(lps),
        '{}_HEAD_PO42'.format(lps),
        '{}_HEAD_PO4'.format(lps),
        '{}_HEAD_COO1'.format(lps),
        '{}_HEAD_COO2'.format(lps),
        '{}_HEAD_COO'.format(lps),
        '{}_PO4'.format(lps),
        '{}_CORE'.format(lps),
        '{}_HEAD'.format(lps),
        '{}_GLYC'.format(lps),
        '{}_CARB'.format(lps),
        '{}_TAIL'.format(lps),
        '{}_TERM'.format(lps),
    ]
    ramp_lps_group_names = [
        '{}_CORE_PO41'.format(lps),
        '{}_CORE_PO42'.format(lps),
        '{}_CORE_PO4'.format(lps)
    ]
    inner_leaflet_group_names = [
        'POPE_HEAD',
        'POPE_PO4',
        'POPE_GLYC',
        'POPE_CARB',
        'POPE_TAIL',
        'POPE_TERM',
        'POPG_HEAD',
        'POPG_PO4',
        'POPG_GLYC',
        'POPG_CARB',
        'POPG_TAIL',
        'POPG_TERM',
        'CDL2_HEAD',
        'CDL2_PO4',
        'CDL2_GLYC',
        'CDL2_CARB',
        'CDL2_TAIL',
        'CDL2_TERM'
    ]
    grouped_group_names = [
        'ALL_HEAD',
        'ALL_PO4',
        'ALL_GLYC',
        'ALL_CARB',
        'ALL_TAIL',
        'ALL_TERM'
    ]
    density_group_names.extend(universal_lps_group_names)
    # if lps == 'RAMP': density_group_names.extend(ramp_lps_group_names)
    # density_group_names.extend(inner_leaflet_group_names)
    # density_group_names.extend(grouped_group_names)
    # density_group_names.extend(ions)

    # Run density
    for repeat in range(repeats):
        for group_name in density_group_names:
            out_density_file = f'{group_name.lower()}_density_r{str(repeat+1)}.xvg'
            out_density_dir = '/'.join([out_dir, 'density', out_density_file])
            gmx.density(xtc_files[repeat], tpr_files[repeat], out_ndx_dirs[repeat], out_density_dir, group_name, 'ALL_TERM')
    print(' > DENSITY ANALYSIS COMPLETE')

    # ====== MEAN SQUARE DISPLACEMENT ======
    msd_group_names = [lps]#, 'POPE', 'POPG', 'CDL2'] # Corresponding to the ndx file
    for repeat in range(repeats):
        for group_name in msd_group_names:
            out_msd_rmcomm_file = f'{group_name.lower()}_msd_rmcomm_r{str(repeat+1)}.xvg'
            out_msd_rmcomm_dir = '/'.join([out_dir, 'msd', out_msd_rmcomm_file])
            gmx.msd_rmcomm(out_no_jump_xtc_dirs[repeat], tpr_files[repeat], out_ndx_dirs[repeat], out_msd_rmcomm_dir, group_name)
    print(' > MEAN SQUARE DISPLACEMENT ANALYSIS COMPLETE')

    # ====== NEAREST NEIGHBOUR ======
    out_nn_files = [f'{lps}_nn_r{str(repeat+1)}.nni' for repeat in range(repeats)]
    out_nn_dirs = ['/'.join([out_dir, 'nearest_neighbour', out_nn_files[repeat]]) for repeat in range(repeats)]
    for repeat in range(repeats):
        nn_indices_per_residue = nn.run(out_sep_gro_dirs[repeat], xtc_files[repeat], lps, frames, 6, 'GM1')
    
    # ====== AREA ======
    out_area_files = [f'{lps}_area_r{str(repeat+1)}.apl' for repeat in range(repeats)]
    out_area_dirs = ['/'.join([out_dir, 'area', out_area_files[repeat]]) for repeat in range(repeats)]
    for repeat in range(repeats):
        area_per_lipid_per_frame = area.run(out_sep_gro_dirs[repeat], xtc_files[repeat], top_files[repeat], lps, frames)

    # ====== ORDER ======
    out_order_files = [f'{lps}_order_r{str(repeat+1)}.ord' for repeat in range(repeats)]
    out_order_dirs = ['/'.join([out_dir, 'order', out_order_files[repeat]]) for repeat in range(repeats)]
    for repeat in range(repeats):
        bond_orders_per_frame, list_of_all_bond_types, corresponding_chains = order.run(out_sep_gro_dirs[repeat], out_no_jump_xtc_dirs[repeat], lps, (0, 0, 1))

    # ====== THICKNESS ======
    out_membrane_thickness_files = [f'thickness_r{str(repeat+1)}.thk' for repeat in range(repeats)]
    out_membrane_thickness_dirs = ['/'.join([out_dir, 'thickness', out_membrane_thickness_files[repeat]]) for repeat in range(repeats)]
    for repeat in range(repeats):
        membrane_thickness_per_frame = thickness.run(out_sep_gro_dirs[repeat], out_no_jump_xtc_dirs[repeat], 'POPE', lps, 'PO4', 'PO1', frames)

    # ====== TORSION ======
    for repeat in range(repeats):
        glcn_glcn_angles, glcn_kdo_angles, kdo_kdo_angles = torsion.run(out_sep_gro_dirs[repeat], out_no_jump_xtc_dirs[repeat])

    # ====== ION CONTACTS ======
    for repeat in range(repeats):
        for ion in ions:
            ion_contacts = contacts.run(out_sep_gro_dirs[repeat], xtc_files[repeat], ion, lps, frames, radius=5.0) # Array: Frame x Bead Type x Residue = # Ion Contacts


if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('command', help='Command to run. Options: [complete, density, no-jump, msd, msd-rmcomm, nearest-neighbour, area, index-complete, order, sep-res, torsion, thickness]')

    gmx_files = ['gro', 'xtc', 'tpr', 'ndx', 'edr', 'top']
    for file_type in gmx_files:
        parser.add_argument(f'-{file_type}', type=str)

    special_arguments = ['out', 'defndx', 'gnm', 'mol', 'frames', 'time', 'nn', 'lps', 'bead', 'xvg', 'il', 'ol', 'ib', 'ob', 'r']
    for arg in special_arguments:
        parser.add_argument(f'-{arg}', type=str)

    repeated_arguments = ['parents', 'gros', 'xtcs', 'tprs', 'tops']
    for arg in repeated_arguments:
        parser.add_argument(f'-{arg}', nargs='+', type=str)

    list_arguments = ['ions']
    for arg in list_arguments:
        parser.add_argument(f'-{arg}', nargs='+', type=str)

    args = parser.parse_args()

    # Run the command
    if args.command == 'density':
        gmx.density(args.xtc, args.tpr, args.ndx, args.out, args.gnm, 'ALL_TERM')
    elif args.command == 'no-jump':
        gmx.no_jump(args.xtc, args.tpr, args.ndx, args.out)
    elif args.command == 'msd':
        gmx.msd(args.xtc, args.tpr, args.ndx, args.out, args.gnm)
    elif args.command == 'msd-rmcomm':
        gmx.msd_rmcomm(args.xtc, args.tpr, args.ndx, args.out, args.gnm)
    elif args.command == 'nearest-neighbour':
        nn.run(args.gro, args.xtc, args.mol, int(args.frames), int(args.nn), bead_name=args.bead)
    elif args.command == 'area':
        area.run(args.gro, args.xtc, args.top, args.mol, int(args.frames))
    elif args.command == 'index-complete':
        index.run(args.gro, args.lps, args.ions, args.defndx, args.out)
    elif args.command == 'thickness':
        thickness.run(args.gro, args.xtc, args.il, args.ol, args.ib, args.ob, int(args.frames))
    elif args.command == 'torsion':
        torsion.run(args.gro, args.xtc, args.out)
    elif args.command == 'order':
        order.run(args.gro, args.xtc, args.mol, (0, 0, 1))
    elif args.command == 'sep-res':
        files.seperate_residues(args.gro, args.mol, args.out)
        print(' SUCCESS: SEPERATED {m} RESIDUES AND WRITTEN TO {o}'.format(m=args.mol, o=args.out))
    elif args.command == 'complete':
        if args.parents is not None:
            gro_files = [parent + '.gro' for parent in args.parents]
            xtc_files = [parent + '.xtc' for parent in args.parents]
            tpr_files = [parent + '.tpr' for parent in args.parents]
        else:
            gro_files = args.gros
            xtc_files = args.xtcs
            tpr_files = args.tprs
        complete(int(args.r), gro_files, xtc_files, tpr_files, args.tops, args.lps, args.ions, int(args.frames), args.out)
    else:
        print('Command not recognized')
