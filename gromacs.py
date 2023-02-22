import subprocess

# Runs gmx density for the given index group along the Z axis taking 100 slices
def density(xtc_file, tpr_file, ndx_file, out_file, ndx_group_name, centre_group_name):
    density_command = "gmx density -f {xtc} -s {tpr} -n {ndx} -d Z -o {out} -center -sl 100 -dens number".format(
        xtc=xtc_file,
        tpr=tpr_file,
        ndx=ndx_file,
        out=out_file)
    density_input = '{i}\n{c}\n'.format(i=ndx_group_name, c=centre_group_name)
    print(' > RUNNING GMX DENSITY FOR {}...'.format(ndx_group_name))
    subprocess.run(density_command, shell=True, input=density_input, encoding='ascii')
    print(' > GMX DENSITY FOR {} COMPLETE'.format(ndx_group_name))

# Creates a no jump trajectory for whole system - important for calculating MSD
def no_jump(xtc_file, tpr_file, ndx_file, out_file):
    trjconv_cmd = 'gmx trjconv -f {xtc} -s {tpr} -n {ndx} -pbc nojump -o {out}'.format(
        xtc=xtc_file,
        tpr=tpr_file,
        ndx=ndx_file,
        out=out_file)
    trjconv_input = 'System\n'
    print(' > CREATING NO JUMP TRAJECTORY...')
    subprocess.run([trjconv_cmd], shell=True, input=trjconv_input, encoding='ascii')
    print(' > NO JUMP TRAJECTORY COMPLETE')

# Runs gmx msd for the given index group in the membrane plane (normal to Z axis)
def msd(xtc_file, tpr_file, ndx_file, out_file, ndx_group_name):
    msd_cmd = 'gmx msd -f {xtc} -s {tpr} -n {ndx} -lateral z -o {out}'.format(
        xtc=xtc_file,
        tpr=tpr_file,
        ndx=ndx_file,
        out=out_file)
    msd_input = '{i}\n'.format(i=ndx_group_name)
    print(' > RUNNING GMX MSD FOR {}...'.format(ndx_group_name))
    subprocess.run([msd_cmd], shell=True, input=msd_input, encoding='ascii')
    print(' > GMX MSD FOR {} COMPLETE'.format(ndx_group_name))

# Runs gmx msd for the given index group in the membrane plane (normal to Z axis) removing the center of mass motion
def msd_rmcomm(xtc_file, tpr_file, ndx_file, out_file, ndx_group_name):
    msd_cmd = 'gmx msd -f {xtc} -s {tpr} -n {ndx} -rmcomm -lateral z -o {out}'.format(
        xtc=xtc_file,
        tpr=tpr_file,
        ndx=ndx_file,
        out=out_file)
    msd_input = '{i}\n{i}\n'.format(i=ndx_group_name)
    print(' > RUNNING GMX MSD RMCOMM FOR {}...'.format(ndx_group_name))
    subprocess.run([msd_cmd], shell=True, input=msd_input, encoding='ascii')
    print(' > GMX MSD RMCOMM FOR {} COMPLETE'.format(ndx_group_name))
