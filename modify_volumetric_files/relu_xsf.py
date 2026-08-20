from tqdm import tqdm
import numpy as np
import sys

z1 = 6.0
z2 = 34.0

filename = sys.argv[1]
output = f"{filename[:-4]}_cut_{int(z1)}_{int(z2)}.xsf"

# ----------------------------
# Read file
# ----------------------------
with open(filename) as f:
    lines = f.readlines()

# Locate data block
for i, line in enumerate(lines):
    if line.strip().startswith("BEGIN_DATAGRID_3D"):
        header_start = i
        break

dims_line = header_start + 1
origin_line = header_start + 2
data_start = header_start + 6

nx, ny, nz = map(int, lines[dims_line].split())
nvalues = nx * ny * nz

# Read volumetric data
values = []
i = data_start
with tqdm(total=nvalues, desc="Reading density", unit="values") as pbar:
    while len(values) < nvalues:
        nums = list(map(float, lines[i].split()))
        values.extend(nums)
        pbar.update(len(nums))
        i += 1
values = np.array(values)

# Reshape
data = values.reshape((nx, ny, nz))

# ----------------------------
# Choose region to remove (Å)
# ----------------------------

origin_z = float(lines[origin_line].split()[2])

# Third lattice vector
cell_z = float(lines[origin_line + 3].split()[2])

# z coordinate corresponding to each z-slice
z = origin_z + np.linspace(0, cell_z, nz)

# Mask slices between z1 and z2
mask = (z >= z1) & (z <= z2)

# Set those slices to zero
data[mask, :, :] = 0.0

# ----------------------------
# Save new file
# ----------------------------

values = data.ravel()

with open(output, "w") as f:
    # Write header
    f.writelines(lines[:data_start])

    # Write volumetric data with progress bar
    for v in tqdm(values, total=len(values), desc="Writing density", unit="values"):
        f.write(f"{v:.5e}\n")

    # Write footer (END_DATAGRID_3D, etc.)
    f.writelines(lines[data_start + nvalues:])

print(f"Saved {output}")
