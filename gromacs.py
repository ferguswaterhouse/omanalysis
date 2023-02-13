import subprocess

def density(xtc_file, tpr_file, ndx_file, out_file, ndx_group_name, density_input):
    # Runs gmx density for the given index group along the Z axis taking 100 slices
    density_command = "gmx density -f {xtc} -s {tpr} -n {ndx} -d Z -o {out} -center -sl 100 -dens number".format(
        xtc=xtc_file,
        tpr=tpr_file,
        ndx=ndx_file,
        out=out_file)
    subprocess.run(density_command, shell=True, input=density_input, encoding='ascii')
