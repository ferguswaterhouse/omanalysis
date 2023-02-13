import argparse
import gromacs as gmx
import neighbour as nn


if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('command', help='Command to run. Options: [complete, msd, density]')

    gmx_files = ['gro', 'xtc', 'tpr', 'ndx', 'edr']
    for file_type in gmx_files:
        parser.add_argument(f'-{file_type}', type=str)

    special_arguments = ['out', 'deffnm', 'gnm', 'mol', 'frames', 'nn', 'bead']
    for arg in special_arguments:
        parser.add_argument(f'-{arg}', type=str)

    args = parser.parse_args()

    # Run the command
    if args.command == 'density':
        gmx.density(args.xtc, args.tpr, args.ndx, args.out, args.gnm, density_input='ALL_TERM\n{gnm}\n'.format(gnm=args.gnm))
    elif args.command == 'nearest-neighbour':
        nn.run(args.gro, args.xtc, args.mol, int(args.frames), int(args.nn), bead_name=args.bead)
    else:
        print('Command not recognized')
