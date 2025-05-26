import os
import subprocess
import numpy as np

def extract_frame(tpr_file, xtc_file, output_gro, time, umbrella):

    command = [
        'gmx', 'trjconv', '-f', xtc_file, '-s', tpr_file, '-o', output_gro,
        '-dump', str(time)
    ]
    process = subprocess.run(command, input='0\n', text=True, capture_output=True)
    if process.returncode != 0:
        print(f"Error extracting frame at time {time}: {process.stderr}")
        return None
    print(f"Frame extracted at time {time} ps: {output_gro}")
    return output_gro

def parse_gro_file(gro_file, atom_indices):

    coordinates = {}
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        for line in lines[2:-1]: 
            parts = line.split()
            try:
                atom_index = int(parts[2])  
                if atom_index in atom_indices:
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    coordinates[atom_index] = np.array([x, y, z])
            except (IndexError, ValueError):
                continue
    return coordinates

def process_filtered_distances(filtered_distances_file):
    with open(filtered_distances_file, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("Frame"):  
                continue

            parts = line.split('\t')
            if len(parts) < 4:
                print(f"Skipping malformed line: {line}")
                continue

            try:
                time = float(parts[0]) 
                umbrella_window = parts[1]  
                tpr_file = f"{umbrella_window}.tpr"
                xtc_file = f"{umbrella_window}.xtc"
                output_filename = f"frame_{time}_{umbrella_window}.gro"

                if os.path.exists(output_filename):
                    print(f"Frame already exists: {output_filename}. Skipping...")
                    continue

                if not (os.path.exists(tpr_file) and os.path.exists(xtc_file)):
                    print(f"Missing files for {umbrella_window}: {tpr_file} or {xtc_file}")
                    continue

                extract_frame(tpr_file, xtc_file, output_filename, time, umbrella_window)

            except Exception as e:
                print(f"Error processing line: {line}, Error: {e}")
                continue

if __name__ == "__main__":
    filtered_distances_file = "filtered_distances_summary.txt"
    process_filtered_distances(filtered_distances_file)

