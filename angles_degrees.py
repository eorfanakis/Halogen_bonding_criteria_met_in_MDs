import re


input_file = 'filtered_angles_output.txt'
output_file = 'filtered_angles_0_20_160_180.txt'

def is_angle_in_range(angle):
    return (0 <= angle <= 20) or (160 <= angle <= 180)

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        match = re.search(r'(\d+\.\d+) degrees', line)
        if match:
            angle = float(match.group(1))  
            if is_angle_in_range(angle):
                outfile.write(line)

print(f"Filtered results have been written to {output_file}")

