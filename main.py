import argparse
import gromacs as gmx
import neighbour as nn
import area
import index
import order
import thickness

# TO DO:
# - Add sugar length analysis
# - Order parameters - Change to MDAnalysis?
# - Add membrane thickness analysis - How to do this? - Average phosphate position for each frame
# - Torsional angle analysis
# - Contact data ??? - Not necessary yet
# - Add them all into one pipeline x 3

if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('command', help='Command to run. Options: [complete, density, np-jump, msd, msd-rmcomm, nearest-neighbour, area, index-complete, order]')

    gmx_files = ['gro', 'xtc', 'tpr', 'ndx', 'edr']
    for file_type in gmx_files:
        parser.add_argument(f'-{file_type}', type=str)

    special_arguments = ['out', 'deffnm', 'gnm', 'mol', 'frames', 'time', 'nn', 'lps', 'bead', 'xvg', 'il', 'ol', 'ib', 'ob']
    for arg in special_arguments:
        parser.add_argument(f'-{arg}', type=str)

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
        area.run(args.edr, args.xvg, 252)
    elif args.command == 'index-complete':
        index.run(args.gro, args.lps, args.ions, args.ndx)
    elif args.command == 'order':
        order.run(args.xtc, args.tpr, args.lps, 0, str(args.time), (0, 0, 1), args.out)
    elif args.command == 'thickness':
        thickness.run(args.gro, args.xtc, args.il, args.ol, args.ib, args.ob, int(args.frames))
    else:
        print('Command not recognized')
