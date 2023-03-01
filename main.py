import argparse
import subprocess
import os
from multiprocessing import Process
import numpy as np
import gromacs as gmx
import neighbour
import area
import index
import thickness
import torsion
import order
import contacts
import files
import process

# TO DO:
# - Sort column format in output files
# - Angle of glucosamine backbone to membrane plane?
# - COG option for Nearest Neighbour?
# - Head group conformation?
# - Complete analysis with averaging
# - Possibly add error processing...
# - Add option to run all scripts directly + give help + jupyter?
    

def single(repeat, gro_file, xtc_file, tpr_file, top_file, lps, ions, frames, out_dir):

    print(f' ========== RUNNING ANALYSIS OF SYSTEM {repeat} ==========')

    # ====== CREATE OUTPUT DIRECTORIES ======
    sub_directories = ['density', 'msd', 'nearest_neighbour', 'area', 'order', 'torsion', 'thickness', 'ion_contacts', 'setup']
    if os.path.exists(out_dir):
        print(' > OUTPUT DIRECTORY ALREADY EXISTS')
        for directory in sub_directories:
            if os.path.exists(f'{out_dir}/{directory}'):
                print(f' > {directory.upper()} DIRECTORY ALREADY EXISTS')
            else:
                subprocess.run(f'mkdir {out_dir}/{directory}', shell=True)
    else:
        subprocess.run(f'mkdir {out_dir}', shell=True)
        for directory in sub_directories:
            subprocess.run(f'mkdir {out_dir}/{directory}', shell=True)
        print(' > CREATED OUTPUT DIRECTORIES')

    # ====== SETUP ======
    out_sep_gro_file = f'sep_r{str(repeat)}.gro'
    out_sep_gro_dir = '/'.join([out_dir, 'setup', out_sep_gro_file])

    out_ndx_file = f'complete_index_r{str(repeat)}.ndx'
    default_ndx_file = f'default_index_r{str(repeat)}.ndx'
    out_ndx_dir = '/'.join([out_dir, 'setup', out_ndx_file])
    default_ndx_dir = '/'.join([out_dir, 'setup', default_ndx_file])

    out_no_jump_xtc_file = f'no_jump_r{str(repeat)}.xtc'
    out_no_jump_xtc_dir = '/'.join([out_dir, 'setup', out_no_jump_xtc_file])

    files.seperate_residues(gro_file, lps, out_sep_gro_dir)
    print(' > DEFINED LPS RESIDUES FOR EACH .GRO FILE')

    index.run(out_sep_gro_dir, lps, ions, default_ndx_dir, out_ndx_dir)
    print(' > CREATED COMPLETE INDEX FILE')

    gmx.no_jump(xtc_file, tpr_file, out_ndx_dir, out_no_jump_xtc_dir)
    print(' > REMOVED JUMPS FROM .XTC FILE')

    # ====== DENSITY ======
    density_group_names = ['W'] # Includes default group names

    density_group_names.extend(index.UNIVERSAL_LPS_GROUP_NAMES(lps))
    if lps == 'RAMP': density_group_names.extend(index.UNIQUE_RAMP_LPS_GROUP_NAMES(lps))
    density_group_names.extend(index.INNER_LEAFLET_GROUP_NAMES)
    density_group_names.extend(index.GROUPED_GROUP_NAMES)
    density_group_names.extend(ions)

    for group_name in density_group_names:
        out_density_file = f'{group_name.lower()}_density_r{str(repeat)}.xvg'
        out_density_dir = '/'.join([out_dir, 'density', out_density_file])
        gmx.density(xtc_file, tpr_file, out_ndx_dir, out_density_dir, group_name, 'ALL_TERM')

    # ====== MEAN SQUARE DISPLACEMENT ======
    msd_group_names = [lps, 'POPE', 'POPG', 'CDL2']

    for group_name in msd_group_names:
        out_msd_rmcomm_file = f'{group_name.lower()}_msd_rmcomm_r{str(repeat)}.xvg'
        out_msd_rmcomm_dir = '/'.join([out_dir, 'msd', out_msd_rmcomm_file])
        gmx.msd_rmcomm(out_no_jump_xtc_dir, tpr_file, out_ndx_dir, out_msd_rmcomm_dir, group_name)
    print(' > MSD ANALYSSIS COMPLETE\n')

    # ====== NEAREST NEIGHBOUR ======
    out_nearest_neighbour_file = f'{lps.lower()}_nn_r{str(repeat)}.nni'
    out_nearest_neighbour_dir = '/'.join([out_dir, 'nearest_neighbour', out_nearest_neighbour_file])
    nearest_neighbour_indices = neighbour.run(out_sep_gro_dir, xtc_file, lps, frames, 6, 'GM1') # Array; shape=(number of molecules); values=Nearest Neighbour Index
    np.savetxt(fname=out_nearest_neighbour_dir, X=nearest_neighbour_indices, header='NEAREST NEIGHBOUR INDICES; shape(#molecules); value=NNI')

    # ====== AREA ======
    out_area_file = f'{lps.lower()}_area_r{str(repeat)}.apl'
    out_area_dir = '/'.join([out_dir, 'area', out_area_file])
    area_per_lipid = area.run(gro_file, xtc_file, top_file, lps, frames) # Array; shape=(number of frames); values=Area per lipid
    np.savetxt(fname=out_area_dir, X=area_per_lipid, header='AREA PER LIPID; shape(2(time, apl), #frames); value=APL')

    # ====== ORDER ======
    out_order_file = f'{lps.lower()}_order_r{str(repeat)}.ord'
    out_order_dir = '/'.join([out_dir, 'order', out_order_file])
    bond_orders_per_frame, list_of_all_bond_types, corresponding_chains = order.run(out_sep_gro_dir, out_no_jump_xtc_dir, lps, (0, 0, 1), frames) # Array; (Bond types, #Residues x #Frames); Value=order parameter)
    np.savetxt(fname=out_order_dir, X=bond_orders_per_frame.T, header='BOND ORDER PARAMETERS; shape(#bond types, #residues x #frames); value=order parameter')

    # ====== THICKNESS ======
    out_membrane_thickness_file = f'thickness_r{str(repeat)}.thk'
    out_membrane_thickness_dir = '/'.join([out_dir, 'thickness', out_membrane_thickness_file])
    membrane_thickness_per_frame = thickness.run(out_sep_gro_dir, out_no_jump_xtc_dir, 'POPE', lps, 'PO4', 'PO1', frames) # Array; shape=(number of frames); values=Membrane thickness
    np.savetxt(fname=out_membrane_thickness_dir, X=membrane_thickness_per_frame.T, header='MEMBRANE THICKNESS; shape(#frames); value=Membrane Thickness (Angstroms)')

    # ====== TORSION ======
    glcn_glcn_file, glcn_kdo_file, kdo_kdo_file = f'{lps.lower()}_glcn_glcn_r{str(repeat)}.agl', f'{lps.lower()}_glcn_kdo_r{str(repeat)}.agl', f'{lps.lower()}_kdo_kdo_r{str(repeat)}.agl'
    glcn_glcn_dir, glcn_kdo_dir, kdo_kdo_dir = '/'.join([out_dir, 'torsion', glcn_glcn_file]), '/'.join([out_dir, 'torsion', glcn_kdo_file]), '/'.join([out_dir, 'torsion', kdo_kdo_file])
    glcn_glcn_angles, glcn_kdo_angles, kdo_kdo_angles = torsion.run(out_sep_gro_dir, out_no_jump_xtc_dir, lps) # Array; shape=(#frames, #residues); values=Torsion Angles
    np.savetxt(fname=glcn_glcn_dir, X=glcn_glcn_angles, header='GLCN-GLCN TORSION ANGLES; shape(#frames, #residues); value=Angle')
    np.savetxt(fname=glcn_kdo_dir, X=glcn_kdo_angles, header='GLCN-KDO TORSION ANGLES; shape(#frames, #residues); value=Angle')
    np.savetxt(fname=kdo_kdo_dir, X=kdo_kdo_angles, header='KDO-KDO TORSION ANGLES; shape(#frames, #residues); value=Angle')

    # ====== ION CONTACTS ======
    for ion in ions:
        out_ion_contacts_file = f'{ion.lower()}_contacts_r{str(repeat)}'
        out_ion_contacts_dir = '/'.join([out_dir, 'ion_contacts', out_ion_contacts_file])
        ion_contacts = contacts.run(out_sep_gro_dir, xtc_file, ion, lps, frames, radius=5.0) # Array; Shape=(#frames, bead types, #residues); values=#Contacts/Bead/Frame
        np.save(file=out_ion_contacts_dir, arr=ion_contacts) # Saved as binary to allow 3D data to be saved
    
    print(f' ========== ANALYSIS OF SYSTEM {repeat} COMPLETE ==========')


def complete(repeats, gro_files, xtc_files, tpr_files, top_files, lps, ions, frames, out_dir):

    print(f' ========== RUNNING ANALYSIS OF ALL SYSTEMS ==========')

    # ====== CREATE OUTPUT DIRECTORIES ======
    sub_directories = ['density', 'msd', 'nearest_neighbour', 'area', 'order', 'torsion', 'thickness', 'ion_contacts', 'setup']
    if os.path.exists(out_dir):
        print(' > OUTPUT DIRECTORY ALREADY EXISTS')
        for directory in sub_directories:
            if os.path.exists(f'{out_dir}/{directory}'):
                print(f' > {directory.upper()} DIRECTORY ALREADY EXISTS')
            else:
                subprocess.run(f'mkdir {out_dir}/{directory}', shell=True)
    else:
        subprocess.run(f'mkdir {out_dir}', shell=True)
        for directory in sub_directories:
            subprocess.run(f'mkdir {out_dir}/{directory}', shell=True)
        print(' > CREATED OUTPUT DIRECTORIES')
    print('\n > ALL SYSTEMS ANALYSED\n')

    print(' > AVERAGING DATA ACROSS SYSTEMS...')

    # ====== RUN PARALLEL ANALYSES OF EACH SYSTEM ======
    system_processes = []
    for repeat in range(repeats):
        single_system_run = Process(target=single, args=(repeat+1, gro_files[repeat], xtc_files[repeat], tpr_files[repeat], top_files[repeat], lps, ions, frames, out_dir))
        system_processes.append(single_system_run)
        single_system_run.start()
    
    for single_system_run in system_processes:
        single_system_run.join()

    # ====== PROCESS DATA FROM ALL SYSTEMS ======

    # === DENSITY ===
    density_files = [f for f in os.listdir('/'.join([out_dir, 'density'])) if os.path.isfile('/'.join([out_dir, 'density', f])) and f.endswith('.xvg')]
    density_file_names, all_density_file_repeats = process.find_common_files(density_files, 'xvg', repeats)

    for i, density_file_repeats in enumerate(all_density_file_repeats):
        densities = [np.loadtxt('/'.join([out_dir, 'density', density_file]), comments=('@', '#'), unpack=True) for density_file in density_file_repeats]
        avg_density = process.average_density_repeats(densities, 100)
        np.savetxt(
            fname='/'.join([out_dir, 'density', f'avg_{density_file_names[i]}.dns']), 
            X=avg_density.T, 
            header='DENSITY: Z, Avg Density, Std')

    # === MSD ===
    msd_files = [f for f in os.listdir('/'.join([out_dir, 'msd'])) if os.path.isfile('/'.join([out_dir, 'msd', f])) and f.endswith('.xvg')]
    msd_file_names, all_msd_file_repeats = process.find_common_files(msd_files, 'xvg', repeats)

    for i, msd_file_repeats in enumerate(all_msd_file_repeats):

        all_msd_data = np.array([np.loadtxt('/'.join([out_dir, 'msd', msd_file]), comments=('@', '#'), unpack=True) for msd_file in msd_file_repeats])
        msd_values = np.array([data[1] for data in all_msd_data])

        diffusion_coefficients = [process.diffusion(msd_data) for msd_data in all_msd_data]
        avg_diffusion_coefficient = np.mean(diffusion_coefficients)
        std_diffusion_coefficient = np.std(diffusion_coefficients)

        time_axis = all_msd_data[0][0]
        avg_msd_over_time = np.mean(msd_values, axis=0)
        std_msd_over_time = np.std(msd_values, axis=0)

        avg_msd_data_over_time = np.array([time_axis, avg_msd_over_time, std_msd_over_time])

        np.savetxt(
            fname='/'.join([out_dir, 'msd', f'avg_{msd_file_names[i]}.msd']), 
            X=avg_msd_data_over_time.T, 
            header='MSD: Time, Avg MSD, Std')
        
        with open('/'.join([out_dir, 'msd', f'avg_{msd_file_names[i]}.txt']), 'w') as f:
            f.write('DIFFUSION COEFFICIENT: Avg, Std\n')
            f.write('\t'.join(['avg', str(avg_diffusion_coefficient)])+'\n')
            f.write('\t'.join(['std', str(std_diffusion_coefficient)])+'\n')

    # === NEAREST NEIGHBOUR ===
    nn_files = [f for f in os.listdir('/'.join([out_dir, 'nearest_neighbour'])) if os.path.isfile('/'.join([out_dir, 'nearest_neighbour', f])) and f.endswith('.nni')]
    nn_file_names, all_nn_file_repeats = process.find_common_files(nn_files, 'nni', repeats)

    for i, nn_file_repeats in enumerate(all_nn_file_repeats):
        all_nn_data = np.array([np.loadtxt('/'.join([out_dir, 'nearest_neighbour', nn_file]), comments=('@', '#'), unpack=True) for nn_file in nn_file_repeats])

        flattened_nn_data = all_nn_data.flatten()
        avg_nni = np.mean(flattened_nn_data)
        std_nni = np.std(flattened_nn_data)

        np.savetxt(
            fname='/'.join([out_dir, 'nearest_neighbour', f'all_{nn_file_names[i]}.nni']), 
            X=flattened_nn_data.T, 
            header='NEAREST NEIGHBOUR: Residue x Nearest Neighbour Index')
        with open('/'.join([out_dir, 'nearest_neighbour', f'avg_{nn_file_names[i]}.txt']), 'w') as f:
            f.write('NEAREST NEIGHBOUR INDEX\n')
            f.write('\t'.join(['avg', str(avg_nni)]) + '\n')
            f.write('\t'.join(['std', str(std_nni)]))

    # === AREA ===
    area_files = [f for f in os.listdir('/'.join([out_dir, 'area'])) if os.path.isfile('/'.join([out_dir, 'area', f])) and f.endswith('.apl')]
    area_file_names, all_area_file_repeats = process.find_common_files(area_files, 'apl', repeats)

    for i, area_file_repeats in enumerate(all_area_file_repeats):
        all_area_data = np.array([np.loadtxt('/'.join([out_dir, 'area', area_file]), comments=('@', '#'), unpack=True).T for area_file in area_file_repeats])
        all_area_values = np.array([data[1] for data in all_area_data])

        time_axis = all_area_data[0][0]
        avg_area_over_time = np.mean(all_area_values, axis=0)
        std_area_over_time = np.std(all_area_values, axis=0)

        avg_area_data_over_time = np.array([time_axis, avg_area_over_time, std_area_over_time])

        flattened_area_data = all_area_values.flatten()
        avg_area = np.mean(flattened_area_data)
        std_area = np.std(flattened_area_data)

        np.savetxt(
            fname='/'.join([out_dir, 'area', f'all_{area_file_names[i]}.apl']), 
            X=avg_area_data_over_time.T, 
            header='AREA PER LIPID: Time, Avg Area, Std')
        with open('/'.join([out_dir, 'area', f'avg_{area_file_names[i]}.txt']), 'w') as f:
            f.write('AREA PER LIPID\n')
            f.write('\t'.join(['avg', str(avg_area)]) + '\n')
            f.write('\t'.join(['std', str(std_area)]))

    # === ORDER PARAMETER ===
    order_files = [f for f in os.listdir('/'.join([out_dir, 'order'])) if os.path.isfile('/'.join([out_dir, 'order', f])) and f.endswith('.ord')]
    order_file_names, all_order_file_repeats = process.find_common_files(order_files, 'ord', repeats)
    for i, order_file_repeats in enumerate(all_order_file_repeats):
        all_order_data = np.array([np.loadtxt('/'.join([out_dir, 'order', order_file]), comments=('@', '#'), unpack=True) for order_file in order_file_repeats])
        flattened_order_data = np.reshape(all_order_data, (all_order_data.shape[1], all_order_data.shape[0] * all_order_data.shape[2]))
        
        avg_bond_orders = np.mean(flattened_order_data, axis=1)
        std_bond_orders = np.std(flattened_order_data, axis=1)

        list_of_corresponding_chains = [chain_i for chain_i, chain in order.CHAIN_TOPOLOGY[lps].items() for i in range(len(order.get_bonds(chain)))]
        chain_bond_orders = {}
        for chain in order.CHAIN_TOPOLOGY[lps].keys():
            chain_bond_orders[chain] = []

        for bond_i in range(len(avg_bond_orders)):
            chain_bond_orders[list_of_corresponding_chains[bond_i]].append(avg_bond_orders[bond_i])

        avg_chain_orders = {}
        for chain, bond_orders in chain_bond_orders.items():
            avg_chain_orders[chain] = np.mean(bond_orders)
            
        with open('/'.join([out_dir, 'order', f'avg_{order_file_names[i]}.txt']), 'w') as f:
            f.write('AVERAGE LIPID TAIL ORDERS: Chain, Order\n')
            for chain, avg_order in avg_chain_orders.items():
                f.write('\t'.join([str(chain), str(avg_order)]) + '\n')

    # === TORSION ===
    torsion_files = [f for f in os.listdir('/'.join([out_dir, 'torsion'])) if os.path.isfile('/'.join([out_dir, 'torsion', f])) and f.endswith('.agl')]
    torsion_file_names, all_torsion_file_repeats = process.find_common_files(torsion_files, 'agl', repeats)
    for i, torsion_file_repeats in enumerate(all_torsion_file_repeats):
        all_torsion_data = np.array([np.loadtxt('/'.join([out_dir, 'torsion', torsion_file]), comments=('@', '#'), unpack=True) for torsion_file in torsion_file_repeats])
        flattened_torsion_data = all_torsion_data.flatten()
        np.savetxt(
            fname='/'.join([out_dir, 'torsion', f'all_{torsion_file_names[i]}.agl']),
            X=flattened_torsion_data.T,
            header='TORSION ANGLES: Angle'
        )

    # === THICKNESS ===
    thickness_files = [f for f in os.listdir('/'.join([out_dir, 'thickness'])) if os.path.isfile('/'.join([out_dir, 'thickness', f])) and f.endswith('.thk')]
    thickness_file_names, all_thickness_file_repeats = process.find_common_files(thickness_files, 'thk', repeats)

    for i, thickness_file_repeats in enumerate(all_thickness_file_repeats):
        all_thickness_data = np.array([np.loadtxt('/'.join([out_dir, 'thickness', thickness_file]), comments=('@', '#'), unpack=True).T for thickness_file in thickness_file_repeats])
        flattened_thickness_data = all_thickness_data.flatten()

        avg_membrane_thickness = np.mean(flattened_thickness_data)
        std_membrane_thickness = np.std(flattened_thickness_data)

        np.savetxt(
            fname='/'.join([out_dir, 'thickness', f'all_{thickness_file_names[i]}.thk']),
            X=flattened_thickness_data.T,
            header='THICKNESS: Thickness'
        )

        with open('/'.join([out_dir, 'thickness', f'avg_{thickness_file_names[i]}.txt']), 'w') as f:
            f.write('MEMBRANE THICKNESS\n')
            f.write('\t'.join(['avg', str(avg_membrane_thickness), '(A)']) + '\n')
            f.write('\t'.join(['std', str(std_membrane_thickness), '(A)']))



if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('command', help='Command to run. Options: [complete, single-system, density, no-jump, msd, msd-rmcomm, nearest-neighbour, area, index-complete, order, sep-res, torsion, thickness]')

    gmx_files = ['gro', 'xtc', 'tpr', 'ndx', 'edr', 'top']
    for file_type in gmx_files:
        parser.add_argument(f'-{file_type}', type=str)

    special_arguments = ['out', 'defndx', 'gnm', 'mol', 'frames', 'time', 'nn', 'lps', 'bead', 'xvg', 'il', 'ol', 'ib', 'ob', 'r', 'parent']
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
        neighbour.run(args.gro, args.xtc, args.mol, int(args.frames), int(args.nn), bead_name=args.bead)
    elif args.command == 'area':
        area.run(args.gro, args.xtc, args.top, args.mol, int(args.frames))
    elif args.command == 'index-complete':
        index.run(args.gro, args.lps, args.ions, args.defndx, args.out)
    elif args.command == 'thickness':
        thickness.run(args.gro, args.xtc, args.il, args.ol, args.ib, args.ob, int(args.frames))
    elif args.command == 'torsion':
        torsion.run(args.gro, args.xtc, args.lps)
    elif args.command == 'order':
        order.run(args.gro, args.xtc, args.mol, (0, 0, 1))
    elif args.command == 'sep-res':
        files.seperate_residues(args.gro, args.mol, args.out)
        print(' SUCCESS: SEPERATED {m} RESIDUES AND WRITTEN TO {o}'.format(m=args.mol, o=args.out))
    elif args.command == 'single-system':
        if args.parent is not None:
            gro_file = args.parent + '.gro'
            xtc_file = args.parent + '.xtc'
            tpr_file = args.parent + '.tpr'
        else:
            gro_file = args.gro
            xtc_file = args.xtc
            tpr_file = args.tpr
        single(int(args.r), gro_file, xtc_file, tpr_file, args.top, args.lps, args.ions, int(args.frames), args.out)
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
