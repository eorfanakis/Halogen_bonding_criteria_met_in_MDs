import os
import re


threshold = 0.3270


groups = ['O21', 'O22', 'O31', 'O32', 'O11', 'O12', 'O13', 'O14']
ndx_files = {
    'O21': 'XB_O21.ndx',
    'O22': 'XB_O22.ndx',
    'O31': 'XB_O31.ndx',
    'O32': 'XB_O32.ndx',
    'O11': 'XB_O11.ndx',
    'O12': 'XB_O12.ndx',
    'O13': 'XB_O13.ndx',
    'O14': 'XB_O14.ndx',
}
output_file = 'filtered_distances_summary.txt'
def get_atom_numbers(index_file, group):
    atoms = []
    with open(index_file, 'r') as idx_file:
        in_group = False
        for line in idx_file:
            if f'[ {group} ]' in line:  
                in_group = True
                continue
            if in_group:
                if '[' in line: 
                    break
                atoms.extend(map(int, line.split()))
    return atoms
with open(output_file, 'w') as outfile:
    outfile.write("Frame\tUmbrella\tGroup\tAtom_Pair\tDistance (nm)\n")
    for group in groups:
        o_atoms = get_atom_numbers(ndx_files[group], group)
        for i in range(51):  # Assuming umbrella0 to umbrella50
            xvg_file = f'dist_halogen_{group}_{i}.xvg'
            if not os.path.exists(xvg_file):
                print(f"File {xvg_file} not found, skipping.")
                continue
            with open(xvg_file, 'r') as infile:
                for line in infile:
                    if line.startswith('@') or line.startswith('#'):
                        continue 
                    values = re.split(r'\s+', line.strip())
                    frame = float(values[0])  
                    distances = values[1:] 

                    for j, dist in enumerate(distances):
                        distance = float(dist)
                        if distance <= threshold:
                            atom_pair = f"CL-{group}({o_atoms[j]})" 
                            outfile.write(f"{frame:.3f}\tumbrella{i}\t{group}\t{atom_pair}\t{distance:.3f}\n")

print(f"Filtered distances saved to {output_file}")

