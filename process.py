import numpy as np
from sklearn.linear_model import LinearRegression


def get_constant_axis(density_axes, bins_n):
    min_z, max_z = -1000, 1000
    for axis in density_axes:
        if axis[-1] < max_z: max_z = axis[-1]
        if axis[0] > min_z: min_z = axis[0]
    return np.linspace(start=min_z, stop=max_z, num=bins_n)


def average_density_repeats(densities, bins):
    avg_density = []
    std_density = []
    z_axis = get_constant_axis([density[0] for density in densities], bins)
    for z in z_axis:
        density_values = []
        for density in densities:
            z_axis = density[0]
            higher_bin = 0
            while z_axis[higher_bin] < z:
                higher_bin += 1
            lower_bin = higher_bin - 1
            if z == z_axis[0]:
                lower_bin = 0
                higher_bin = 1
            interpolation = (z - z_axis[lower_bin]) / (z_axis[higher_bin] - z_axis[lower_bin])
            density_values.append(density[1][lower_bin] + interpolation * (density[1][higher_bin] - density[1][lower_bin]))
        avg_density.append(np.mean(density_values))
        std_density.append(np.std(density_values))
    return np.array([z_axis, avg_density, std_density])


def find_common_files(files, ftype_sfx, repeats):
    file_pfxs = []
    for file in files:
        for sfx in [f'_r{r+1}.{ftype_sfx}' for r in range(repeats)]:
            if file.endswith(sfx): 
                file_pfxs.append(file.replace(sfx, ''))
    unique_file_pfxs = list(set(file_pfxs))
    file_groups = []
    for pfx in unique_file_pfxs:
        file_groups.append([f'{pfx}_r{r+1}.{ftype_sfx}' for r in range(repeats)])
    return unique_file_pfxs, file_groups


def diffusion(msd_data):
    lower_quartile = int(len(msd_data[0]) / 4)
    upper_quartile = 3 * lower_quartile
    x = msd_data[0][lower_quartile: upper_quartile].reshape(-1, 1)
    y = msd_data[1][lower_quartile: upper_quartile]
    model = LinearRegression().fit(x, y)
    diff = model.coef_[0] / 4 # nm2/ns
    diff_cm = diff * 10**9 * (10**-7)**2
    return diff_cm
