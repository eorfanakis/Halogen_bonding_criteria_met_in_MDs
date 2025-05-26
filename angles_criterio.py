import re
import os
import numpy as np
import MDAnalysis as mda

carbon_atom = 19
cl_atom = 20

def load_filtered_triplets(distances_file, carbon_atom, cl_atom):
    triplets = {}
    with open(distances_file, 'r') as f:
        next(f)  
        for line in f:
            parts = line.strip().split('\t')
            frame = float(parts[0])
            umbrella = int(parts[1].replace("umbrella", ""))
            atom_pair = parts[3]

            match = re.search(r'CL-O\d+\((\d+)\)', atom_pair)
            if match:
                o_atom = int(match.group(1))
                key = (frame, umbrella)
                if key not in triplets:
                    triplets[key] = []
                triplets[key].append([carbon_atom, cl_atom, o_atom])
    return triplets

def calculate_angle(atom1, atom2, atom3):
    v1 = atom1 - atom2
    v2 = atom3 - atom2
    angle = np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))
    return np.degrees(angle)

distances_file = 'filtered_distances_summary.txt'
output_file = 'filtered_angles_output.txt'

filtered_triplets = load_filtered_triplets(distances_file, carbon_atom, cl_atom)

with open(output_file, "w") as f_out:
    for gro_file in sorted([f for f in os.listdir('.') if f.startswith("frame_") and f.endswith(".gro")]):
        match = re.search(r'frame_(\d+\.\d+)_umbrella(\d+).gro', gro_file)
        if match:
            frame = float(match.group(1))
            umbrella = int(match.group(2))
            key = (frame, umbrella)

            if key in filtered_triplets:
                u = mda.Universe(gro_file)
                for triplet in filtered_triplets[key]:
                    try:
                        atom1 = u.atoms.select_atoms(f'bynum {triplet[0]}').positions[0]
                        atom2 = u.atoms.select_atoms(f'bynum {triplet[1]}').positions[0]
                        atom3 = u.atoms.select_atoms(f'bynum {triplet[2]}').positions[0]

                        angle = calculate_angle(atom1, atom2, atom3)
                        f_out.write(f"Angle between atoms {triplet} in {gro_file}: {angle:.2f} degrees\n")
                    except IndexError:
                        print(f"Atoms {triplet} not found in {gro_file}. Skipping.")

print(f"Angles have been successfully calculated and written to {output_file}")

