import argparse
import gromacs as gmx
import neighbour as nn
import area
import index

# TO DO:
# - Order parameters
# - Add membrane thickness analysis
# - Contact data ???
# - Add sugar length analysis
# - Add them all into one pipeline x 3


if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('command', help='Command to run. Options: [complete, msd, density]')

    gmx_files = ['gro', 'xtc', 'tpr', 'ndx', 'edr']
    for file_type in gmx_files:
        parser.add_argument(f'-{file_type}', type=str)

    special_arguments = ['out', 'deffnm', 'gnm', 'mol', 'frames', 'nn', 'lps', 'bead', 'xvg']
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
    else:
        print('Command not recognized')
