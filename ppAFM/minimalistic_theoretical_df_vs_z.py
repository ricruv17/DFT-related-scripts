"""
Interactive XSF Volumetric Data Viewer

This script loads and visualizes 3D volumetric data from an XCrySDen df.xsf file
containing frequency-shift values (e.g., from NC-AFM simulations or experiments).
It provides an interactive Matplotlib interface to explore the data in both
XY-plane slices and Z-line spectra.

Main features:
- Parses df.xsf volumetric data into a structured pandas DataFrame containing
  spatial coordinates (x, y, z) and frequency shift values.
- Caches parsed data to CSV for faster reloads on subsequent runs.
- Scroll through Z slices (mouse wheel) to view XY-plane images at different heights.
- Click on any point in the XY-plane to plot the corresponding frequency-vs-Z spectrum
  in an adjacent subplot.
- Supports multiple clicked points with unique colors and a legend for easy comparison.

Usage:
- Place this script in the same directory as your target df.xsf file.
- Run it directly; the first df.xsf file found will be loaded.
- Use the mouse wheel to change Z slices, and left-click to select points for Z-line plots.

Requirements:
- numpy, pandas, matplotlib
- df.xsf file output from ppAFM (https://github.com/Probe-Particle/ppafm) program 

Author: Ricardo Ruvalcaba & Faisal Almuhaisen
Date: 2025-08-14
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pandas as pd
from tkinter import filedialog
from tqdm import tqdm

def export_to_csv(ax, x_label, y_label):
    lines = ax.get_lines()
    if not lines:
        print("No lines to export.")
        return

    export_data = {}

    for line in lines:
        label = line.get_label()
        x = line.get_xdata()
        y = line.get_ydata()
        export_data[f"{label}-{x_label}"] = pd.Series(x)
        export_data[f"{label}-{y_label}"] = pd.Series(y)

    df = pd.concat(export_data, axis=1)

    # Ask the user where to save the file
    filepath = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if filepath:
        df.to_csv(filepath, index=False)
        print(f"Exported to {filepath}")


def parse_xsf_minimal(filename):
    print('Opening file...')
    with open(filename, "r") as f:
        lines = f.readlines()

    # Locate volumetric block
    for i, line in enumerate(lines):
        if "BEGIN_DATAGRID_3D" in line:
            start = i + 1
            break
    else:
        raise RuntimeError("BEGIN_DATAGRID_3D not found")

    # Grid size
    nx, ny, nz = map(int, lines[start].split())

    # Origin and lattice vectors
    origin = np.array(list(map(float, lines[start + 1].split())))
    a1 = np.array(list(map(float, lines[start + 2].split())))
    a2 = np.array(list(map(float, lines[start + 3].split())))
    a3 = np.array(list(map(float, lines[start + 4].split())))
    
    # Read scalar values
    scalar_values = []
    data_lines = lines[start + 5:]

    for line in tqdm(
        data_lines,
        desc="Reading voxels",
        unit="lines"
    ):
        if "END_DATAGRID_3D" in line:
            break
        scalar_values.extend(map(float, line.split()))

    data = np.array(scalar_values, dtype=float).reshape((nz, ny, nx))

    print('Building coordinate axes...')
    # Build coordinate axes
    x = np.linspace(origin[0], origin[0] + a1[0], nx)
    y = np.linspace(origin[1], origin[1] + a2[1], ny)
    z = np.linspace(origin[2], origin[2] + a3[2], nz)

    return data, x, y, z


def df_z_at_xy(data, x, y, z, x_target, y_target):
    ix = np.argmin(np.abs(x - x_target))
    iy = np.argmin(np.abs(y - y_target))

    df_z = data[:, iy, ix]

    return z, df_z, x[ix], y[iy]


def on_key(event):
    if event.key.lower() == 'd':
        print("Exporting Z-line spectra to CSV...")
        export_to_csv(ax, x_label="Z", y_label="df")


# === Main Execution ===
data, x, y, z = parse_xsf_minimal("df.xsf")

z_1, df_1, x_used, y_used = df_z_at_xy(
    data, x, y, z,
    x_target=Au111_x,
    y_target=Au111_y
)

z_2, df_2, x_used, y_used = df_z_at_xy(
    data, x, y, z,
    x_target=organoAu_x,
    y_target=organoAu_y
)

print('Plotting...')
fig, ax = plt.subplots()
ax.plot(z_1, df_1, label='Au111')
ax.plot(z_2, df_2, label='organoAu')
plt.legend()
fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()

