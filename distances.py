import os


base_command = "gmx distance -f {xtc_file} -s umbrella{num}.tpr -n XB_{group}.ndx -select 2 -oall dist_halogen_{group}_{num}.xvg"


groups = ['O21', 'O22', 'O31', 'O32', 'O11', 'O12', 'O13', 'O14']


num_files = 51  # umbrella0 to umbrella50 so 51 Umbrella windows


for i in range(num_files):
    xtc_file = f"umbrella{i}.xtc"
    for group in groups:
        command = base_command.format(xtc_file=xtc_file, group=group, num=i)
        os.system(command)

print("All commands have been executed.")

