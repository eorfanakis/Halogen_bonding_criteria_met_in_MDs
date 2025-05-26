import re
import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

# File paths
com_file = 'alp_com_z_corrected.txt'
angles_file = 'filtered_angles_0_20_160_180.txt'
distances_file = 'filtered_distances_summary.txt'

# Load data functions
def load_com_data(com_file):
    com_data = {}
    with open(com_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.split()
            umbrella = int(re.search(r'umbrella(\d+)', parts[0]).group(1))
            com_z = float(parts[1].replace('nm', '')) * 2  # Multiply by 2 as requested
            com_data[umbrella] = com_z
    return com_data

def load_angle_data(angles_file):
    angle_data = []
    with open(angles_file, 'r') as f:
        for line in f:
            match = re.search(r'Angle between atoms \[(\d+), (\d+), (\d+)\] in frame_(\d+\.\d+)_umbrella(\d+)\.gro: ([\d\.]+)', line)
            if match:
                frame_time, umbrella = float(match.group(4)), int(match.group(5))
                angle = float(match.group(6))
                # Adjust angles from 0-20 degrees to 180-angle
                if angle <= 20:
                    angle = 180 - angle
                angle_data.append((frame_time, umbrella, angle))
    return angle_data

def load_distance_data(distances_file, angle_data):
    distance_data = {}
    with open(distances_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            frame_time = float(parts[0])
            umbrella = int(parts[1].replace('umbrella', ''))
            distance = float(parts[4])
            for angle_frame_time, angle_umbrella, _ in angle_data:
                if frame_time == angle_frame_time and umbrella == angle_umbrella:
                    distance_data[(frame_time, umbrella)] = distance
    return distance_data

# Load all data
com_data = load_com_data(com_file)
angle_data = load_angle_data(angles_file)
distance_data = load_distance_data(distances_file, angle_data)

# Prepare data for plotting
x_vals, y_vals, z_vals, colors = [], [], [], []

for frame_time, umbrella, angle in angle_data:
    key = (frame_time, umbrella)
    if key in distance_data and umbrella in com_data:
        x = com_data[umbrella]
        y = distance_data[key]
        z = angle
        x_vals.append(x)
        y_vals.append(y)
        z_vals.append(z)
        # All dots blue after 3.6 nm
        colors.append('blue' if x > 3.6 else 'purple')

mean_distance = np.mean(y_vals)
mean_angle = np.mean(z_vals)

# === 3D Plot ===
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x_vals, y_vals, z_vals, c=colors, marker='o', depthshade=True, s=120)

ax.set_xlabel('Distance from center of bilayer (nm)')
ax.set_ylabel('Halogen bond distance (nm)')
ax.set_zlabel('C - Cl - O Angle (degrees)')

# ✅ Show integers only on the Z (angle) axis, leave other axes untouched
ax.set_zticklabels([f"{int(t)}" for t in ax.get_zticks()])

# Text annotations
textstr_left = f'Geometrical criteria met: {len(x_vals)}\nUmbrella Sampling windows: 51\nMean Distance: {mean_distance:.2f} nm\nMean Angle: {mean_angle:.2f}°'
textstr_right = 'Purple dots: DOPC Backbone Oxygens\nBlue dots: DOPC Phosphate Oxygens'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text2D(0.05, 0.95, textstr_left, transform=ax.transAxes, fontsize=10, bbox=props)
ax.text2D(0.70, 0.95, textstr_right, transform=ax.transAxes, fontsize=10, bbox=props)

# Save static image
plt.savefig('3D-plot.png', dpi=1200)

# === 2D Plot ===
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Normalize color based on angle proximity to 180
angle_proximity = [180 - abs(180 - angle) for angle in z_vals]
scatter2 = ax2.scatter(x_vals, y_vals, c=angle_proximity, cmap='Purples', edgecolors='k')

ax2.set_xlabel('Distance from center of bilayer (nm)')
ax2.set_ylabel('Halogen bond distance (nm)')
ax2.set_title('2D Plot of Alprazolam')

cbar = plt.colorbar(scatter2, ax=ax2)
cbar.set_label('Angle proximity to 180° (degrees)')

plt.savefig('modified_2d_plot.png', dpi=1200)
plt.show()

