import numpy as np
import gromacs


def calculate_area(boxxy_file):
    boxxy = np.loadtxt(boxxy_file, comments=('#', '@'), unpack=True)
    x_values, y_values = boxxy[1], boxxy[2]
    area = x_values * y_values
    return area


def run(edr_file, boxxy_file, number_of_lipids):
    gromacs.boxxy(edr_file, boxxy_file)
    area_per_lipid_over_time = calculate_area(boxxy_file) / number_of_lipids
    return area_per_lipid_over_time
