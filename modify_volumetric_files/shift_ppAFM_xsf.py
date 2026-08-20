import numpy as np
import sys

filename = sys.argv[1]
shift_amount = sys.argv[2]
output = f"{filename[:-4]}_shifted_{shift_amount}AA.xsf"


# Read file
with open(filename) as f:
    lines = f.readlines()

# Find beginning of data
for i, line in enumerate(lines):
    if line.strip().startswith("BEGIN_DATAGRID_3D"):
        header_start = i
        break

dims_line = header_start + 1
origin_line = header_start + 2
vectors_end = header_start + 5
data_start = header_start + 6

nx, ny, nz = map(int, lines[dims_line].split())

nvalues = nx * ny * nz

# Read all values
values = np.array(
    [float(x) for x in lines[data_start:data_start+nvalues]]
)

# XSF DATAGRID_3D is FORTRAN/column-major order
data = values.reshape((nx, ny, nz), order='F')

# Z grid spacing
cell_z = float(lines[origin_line+3].split()[2])
dz = cell_z / nz

shift = round(float(shift_amount) / dz)

# Shift along Z
data = np.roll(data, shift=shift, axis=2)

# Flatten back in XSF's required FORTRAN order
values = data.ravel(order='F')

# Replace data
new_lines = lines[:data_start]
new_lines.extend(f"{v:.5e}\n" for v in values)
new_lines.extend(lines[data_start+nvalues:])

with open(output, "w") as f:
    f.writelines(new_lines)
