import MDAnalysis as mda
import numpy as np
import os

gro_files = [f'umbrella{i}.gro' for i in range(51)]  

def calculate_com_z(universe):
    alp_atoms = universe.select_atoms('resname ALP')  
    if len(alp_atoms) == 0:
        raise ValueError("No atoms found for ALP molecule in the given file.")
    
    com_z = alp_atoms.center_of_mass()[2] 
    com_z_nm = com_z / 10 
    return com_z_nm


output_file = 'alp_com_z_corrected.txt'
with open(output_file, 'w') as f:
    f.write("File\tCOM_z/2 (nm)\n")  

    for gro_file in gro_files:
        if os.path.exists(gro_file):
            try:
                u = mda.Universe(gro_file)
                com_z_nm = calculate_com_z(u)
                com_z_nm_div2 = com_z_nm / 2 
                print(f"COM of ALP molecule in {gro_file} (z-axis / 2, nm): {com_z_nm_div2:.4f} nm")
                
                f.write(f"{gro_file}\t{com_z_nm_div2:.4f} nm\n")
            except Exception as e:
                print(f"Error processing {gro_file}: {e}")
        else:
            print(f"{gro_file} not found.")

print(f"Results saved to {output_file}")

