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


def boxxy(edr_file, out_file):
    energy_cmd = "gmx energy -f {edr} -o {xvg}".format(
        edr=edr_file, 
        xvg=out_file)
    energy_input = '12\n13\n\n'
    print(' > GETTING BOX X AND Y DIMENSIONS...')
    subprocess.run([energy_cmd], shell=True, input=energy_input, encoding='ascii')
    print(' > BOX X AND Y DIMENSIONS RETRIEVED...')


def index(gro_file, out_file, index_input):
    index_cmd = 'gmx make_ndx -f {gro} -o {ndx}'.format(
        gro=gro_file,
        ndx=out_file)
    print(' > CREATING INDEX FILE...')
    subprocess.run([index_cmd], shell=True, input=index_input, encoding='ascii')
    print(' > INDEX FILE CREATED...')


def frame_dump(xtc_file, tpr_file, initial_time, final_time, frame_skip, lipid, outdir, file_name):

    cmd = 'gmx trjconv -f {xtc} -s {tpr} -b {int} -e {fit} -sep -skip {skp} -pbc whole -o {odr}/{fnm}_.gro > /dev/null'.format(
        xtc=xtc_file,
        tpr=tpr_file,
        int=initial_time,
        fit=final_time,
        skp=frame_skip,
        odr=outdir,
        fnm=file_name
    )
    subprocess.run([cmd], shell=True, input=lipid, encoding='ascii')
