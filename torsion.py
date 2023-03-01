import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals


def run(gro_dir, xtc_dir, lps):

    print(f' > RUNNING TORSION ANGLE ANALYSIS OF {lps}')

    u = mda.Universe(gro_dir, xtc_dir)

    # Defines bonds
    glcn_glcn = "name GM3 GM1 GM5 GM4"
    glcn_kdo = "name GM1 GM3 S02 S03"
    kdo_kdo = "name S07 S06 S03 S02"

    all_residues_all_atoms = u.select_atoms(f'resname {lps}').residues

    # Selects all bonds over all residues
    all_lps_glcn_glcn = [residue.atoms.select_atoms(glcn_glcn) for residue in all_residues_all_atoms]
    all_lps_glcn_kdo = [residue.atoms.select_atoms(glcn_kdo) for residue in all_residues_all_atoms]
    all_lps_kdo_kdo = [residue.atoms.select_atoms(kdo_kdo) for residue in all_residues_all_atoms]

    # Calculates torsion angles
    glcn_glcn_angles = dihedrals.Dihedral(all_lps_glcn_glcn)
    glcn_kdo_angles = dihedrals.Dihedral(all_lps_glcn_kdo)
    kdo_kdo_angles = dihedrals.Dihedral(all_lps_kdo_kdo)

    glcn_glcn_angles.run()
    glcn_kdo_angles.run()
    kdo_kdo_angles.run()

    print(f' > TORSION ANGLE ANALYSIS OF {lps} COMPLETE\n')

    return glcn_glcn_angles.results.angles, glcn_kdo_angles.results.angles, kdo_kdo_angles.results.angles # Array; shape=(#frames, #residues); values=Torsion Angles
