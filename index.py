import gromacs
import files


def learn_default_index_groups(gro_file, default_ndx_file):
    gromacs.index(gro_file, default_ndx_file, 'q\n')
    return files.read_index(default_ndx_file)


def create_all_groups(start_index, lps, ions):

    # Gives the command for each group
    ndx_group_commands = {
        '{}_HEAD_PO41'.format(lps): '"{}" & a PO1'.format(lps),
        '{}_HEAD_PO42'.format(lps): '"{}" & a PO2'.format(lps),
        '{}_HEAD_PO4'.format(lps):  '"{}" & a PO1 PO2'.format(lps),
        '{}_HEAD_COO1'.format(lps): '"{}" & a S01'.format(lps),
        '{}_HEAD_COO2'.format(lps): '"{}" & a S10'.format(lps),
        '{}_HEAD_COO'.format(lps):  '"{}" & a S01 S10'.format(lps),
        '{}_PO4'.format(lps):       '"{}" & a P*'.format(lps),
        '{}_CORE'.format(lps):      '"{}" & a S*'.format(lps),
        '{}_HEAD'.format(lps):      '"{}" & a GM*'.format(lps),
        '{}_GLYC'.format(lps):      '"{}" & a GL*'.format(lps),
        '{}_CARB'.format(lps):      '"{}" & a C*'.format(lps),
        '{}_TAIL'.format(lps):      '"{}_GLYC" | "{}_CARB"'.format(lps, lps),
        '{}_TERM'.format(lps):      '"{}" & a C3A C3B C3D C3C C2E C2F'.format(lps),
        'POPE_HEAD':                '"POPE" & a NH3',
        'POPG_HEAD':                '"POPG" & a GL0',
        'CDL2_HEAD':                '"CDL2" & a GL0',
        'POPE_PO4':                 '"POPE" & a PO4',
        'POPG_PO4':                 '"POPG" & a PO4',
        'CDL2_PO4':                 '"CDL2" & a P*',
        'POPE_GLYC':                '"POPE" & a GL*',
        'POPG_GLYC':                '"POPG" & a GL1 GL2',
        'CDL2_GLYC':                '"CDL2" & a GL11 GL12 GL21 GL22',
        'POPE_CARB':                '"POPE" & a C* D*',
        'POPG_CARB':                '"POPG" & a C* D*',
        'CDL2_CARB':                '"CDL2" & a C* D*',
        'POPE_TAIL':                '"POPE_GLYC" | "POPE_CARB"',
        'POPG_TAIL':                '"POPG_GLYC" | "POPG_CARB"',
        'CDL2_TAIL':                '"CDL2_GLYC" | "CDL2_CARB"',
        'POPE_TERM':                '"POPE" & a C4A C4B',
        'POPG_TERM':                '"POPG" & a C4A C4B',
        'CDL2_TERM':                '"CDL2" & a C5A1 C5B1 C5A2 C5B2',
        'ALL_HEAD':                 '"{}_HEAD" | "POPE_HEAD" | "POPG_HEAD" | "CDL2_HEAD"'.format(lps),
        'ALL_PO4':                  '"{}_PO4" | "POPE_PO4" | "POPG_PO4" | "CDL2_PO4"'.format(lps),
        'ALL_GLYC':                 '"{}_GLYC" | "POPE_GLYC" | "POPG_GLYC" | "CDL2_GLYC"'.format(lps),
        'ALL_CARB':                 '"{}_CARB" | "POPE_CARB" | "POPG_CARB" | "CDL2_CARB"'.format(lps),
        'ALL_TAIL':                 '"{}_TAIL" | "POPE_TAIL" | "POPG_TAIL" | "CDL2_TAIL"'.format(lps),
        'ALL_TERM':                 '"{}_TERM" | "POPE_TERM" | "POPG_TERM" | "CDL2_TERM"'.format(lps)
    }

    ndx_groups_to_define = [] # List of groups to define in command

    # List of possible index groups
    universal_lps_group_names = [
        '{}_HEAD_PO41'.format(lps),
        '{}_HEAD_PO42'.format(lps),
        '{}_HEAD_PO4'.format(lps),
        '{}_HEAD_COO1'.format(lps),
        '{}_HEAD_COO2'.format(lps),
        '{}_HEAD_COO'.format(lps),
        '{}_PO4'.format(lps),
        '{}_CORE'.format(lps),
        '{}_HEAD'.format(lps),
        '{}_GLYC'.format(lps),
        '{}_CARB'.format(lps),
        '{}_TAIL'.format(lps),
        '{}_TERM'.format(lps),
    ]

    inner_leaflet_group_names = [
        'POPE_HEAD',
        'POPE_PO4',
        'POPE_GLYC',
        'POPE_CARB',
        'POPE_TAIL',
        'POPE_TERM',
        'POPG_HEAD',
        'POPG_PO4',
        'POPG_GLYC',
        'POPG_CARB',
        'POPG_TAIL',
        'POPG_TERM',
        'CDL2_HEAD',
        'CDL2_PO4',
        'CDL2_GLYC',
        'CDL2_CARB',
        'CDL2_TAIL',
        'CDL2_TERM'
    ]

    grouped_group_names = [
        'ALL_HEAD',
        'ALL_PO4',
        'ALL_GLYC',
        'ALL_CARB',
        'ALL_TAIL',
        'ALL_TERM'
    ]

    # Adds the group names to a list of groups to define
    ndx_groups_to_define.extend(universal_lps_group_names)
    ndx_groups_to_define.extend(inner_leaflet_group_names)
    ndx_groups_to_define.extend(grouped_group_names)

    ndx_commands = ['\n'.join([ndx_group_commands[ndx_group_name], 'name {} {}'.format(str(start_index+i), ndx_group_name)]) for i, ndx_group_name in enumerate(ndx_groups_to_define)]
    ndx_commands.append('q\n')
    ndx_input = '\n'.join(ndx_commands)

    return ndx_input


def run(gro_file, default_ndx_file):
    default_index_groups = learn_default_index_groups(gro_file, default_ndx_file)
    ndx_input = create_all_groups(len(default_index_groups), 'REMP', ('NA', 'CA'))
    gromacs.index(gro_file, 'test.ndx', ndx_input)
