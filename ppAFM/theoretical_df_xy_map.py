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

def parse_xsf_volumetric(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "BEGIN_DATAGRID_3D_whatever" in line:
            start = i + 1
            break

    nx, ny, nz = map(int, lines[start].split())
    origin = np.array(list(map(float, lines[start + 1].split())))
    a1 = np.array(list(map(float, lines[start + 2].split())))
    a2 = np.array(list(map(float, lines[start + 3].split())))
    a3 = np.array(list(map(float, lines[start + 4].split())))
    voxel_vector = (nx, ny, nz)
    lattice_vectors = (a1, a2, a3)

    data_lines = lines[start + 5:]
    scalar_values = []
    for line in tqdm(data_lines, desc="Parsing scalar values", unit="line"):
        if "END_DATAGRID_3D" in line:
            break
        scalar_values.extend(map(float, line.split()))

    data = np.array(scalar_values).reshape((nz, ny, nx))

    x_coords = np.linspace(origin[0], origin[0] + a1[0], nx)
    y_coords = np.linspace(origin[1], origin[1] + a2[1], ny)
    dz = round(a3[2]/nz, -int(np.floor(np.log10(abs(a3[2]/nz))))) # round to significant digits (numerical issue)
    z_coords = origin[2] + 4.5 + dz * np.arange(nz)

    proportionality_constant = 1 #80/0.567 #empirically extracted
    coords = []
    for iz in tqdm(range(nz), desc="Building coordinates", unit="slice"):
        for iy in range(ny):
            for ix in range(nx):
                coords.append([z_coords[iz], y_coords[iy], x_coords[ix], data[iz, iy, ix]*proportionality_constant])

    df_table = pd.DataFrame(coords, columns=["z", "y", "x", "frequency"])
    df_table = df_table.sort_values(by=["x", "y", "z"])

    return voxel_vector, lattice_vectors, df_table

def load_or_parse_dataframe(xsf_path):
    # Determine CSV filename from XSF filename
    base_name = os.path.splitext(os.path.basename(xsf_path))[0]
    csv_filename = f"{base_name}.csv"

    if os.path.exists(csv_filename):
        print(f"Loading cached DataFrame from '{csv_filename}'...")
        df_table = pd.read_csv(csv_filename)
        voxel_vector = lattice_vectors = None  # Skip voxel/lattice unless needed
    else:
        print(f"No cached DataFrame found. Parsing '{xsf_path}'...")
        voxel_vector, lattice_vectors, df_table = parse_xsf_volumetric(xsf_path)
        df_table.to_csv(csv_filename, index=False)
        print(f"DataFrame saved to '{csv_filename}'")

    return voxel_vector, lattice_vectors, df_table

def get_closest_z_line(df_table, x_target, y_target):
    df_table['distance'] = np.sqrt((df_table['x'] - x_target)**2 + (df_table['y'] - y_target)**2)
    closest_row = df_table.loc[df_table['distance'].idxmin()]
    x_closest = closest_row['x']
    y_closest = closest_row['y']
    z_line = df_table[(df_table['x'] == x_closest) & (df_table['y'] == y_closest)].copy()
    z_line = z_line.sort_values(by='z')
    z_line = z_line.iloc[:-1]  # optional: remove last point
    df_table.drop(columns='distance', inplace=True)
    return z_line, (x_closest, y_closest)

def plot_z_line(z_line, x_found, y_found, color):
    if z_line.empty:
        print("Z-line is empty. Nothing to plot.")
        return

    fig_z, ax = plt.subplots()
    ax.plot(z_line['z'], z_line['frequency'], color=color, label=f"(x={x_found:.4f}, y={y_found:.4f})")
    ax.set_xlabel("Z (Å)")
    ax.set_ylabel("Frequency shift (Hz)")
    ax.set_title(f"df vs. Z at (x={x_found:.4f}, y={y_found:.4f})")
    ax.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

def interactive_xy_plane_scroll(df_table, initial_z=None, n_colors=10):
    z_values = np.sort(df_table['z'].unique())
    z_index = 0
    if initial_z is not None:
        z_index = (np.abs(z_values - initial_z)).argmin()
    current_z = z_values[z_index]

    df_zslice = df_table[np.isclose(df_table['z'], current_z)].copy()
    pivot_table = df_zslice.pivot(index='y', columns='x', values='frequency')
    pivot_table = pivot_table.sort_index(ascending=True).sort_index(axis=1, ascending=True)

    # --- Create a single figure with two subplots ---
    fig, (ax_z, ax_xy) = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

    # XY-plane image
    img = ax_xy.imshow(
        pivot_table.values,
        extent=[
            pivot_table.columns.min(),
            pivot_table.columns.max(),
            pivot_table.index.min(),
            pivot_table.index.max()
        ],
        origin='lower',
        aspect='auto',
        cmap='gray'
    )
    ax_xy.set_xlabel("X (Å)")
    ax_xy.set_ylabel("Y (Å)")
    title_xy = ax_xy.set_title(f"Frequency shift at Z ≈ {current_z:.4f} Å")
    cbar = fig.colorbar(img, ax=ax_xy, label="Frequency shift (Hz)")
    
    # Z-line subplot
    ax_z.set_xlabel("Z (Å)")
    ax_z.set_ylabel("Frequency shift (Hz)")
    ax_z.set_title("Z-line spectra for clicked points\nPress \"d\" to export")

    # Colormap for multiple points
    colormap = plt.get_cmap('nipy_spectral', n_colors)
    color_list = [colormap(i) for i in range(n_colors)]
    color_index = 0

    # Markers on XY-plane
    markers_data = []  # (x, y, color)
    marker = ax_xy.scatter([], [], s=20, c=[], marker='o')

    def update_xy_plot(z_index):
        nonlocal current_z
        current_z = z_values[z_index]
        df_zslice = df_table[np.isclose(df_table['z'], current_z)].copy()
        if not df_zslice.empty:
            pivot_table = df_zslice.pivot(index='y', columns='x', values='frequency')
            pivot_table = pivot_table.sort_index(ascending=True).sort_index(axis=1, ascending=True)
            img.set_data(pivot_table.values)
            img.set_extent([
                pivot_table.columns.min(),
                pivot_table.columns.max(),
                pivot_table.index.min(),
                pivot_table.index.max()
            ])
            img.set_clim(vmin=pivot_table.values.min(), vmax=pivot_table.values.max())
            title_xy.set_text(f"Frequency shift at Z ≈ {current_z:.4f} Å")
            fig.canvas.draw_idle()

    def on_scroll(event):
        nonlocal z_index
        if event.button == 'up':
            z_index = min(z_index + 1, len(z_values) - 1)
        elif event.button == 'down':
            z_index = max(z_index - 1, 0)
        update_xy_plot(z_index)

    def on_click(event):
        nonlocal color_index
        if event.inaxes != ax_xy:
            return
        x_click = event.xdata
        y_click = event.ydata
        if x_click is None or y_click is None:
            return

        # Get closest z-line
        z_line, (x_found, y_found) = get_closest_z_line(df_table, x_click, y_click)

        # Assign color
        color = color_list[color_index % n_colors]
        color_index += 1

        # Update XY-plane markers
        markers_data.append((x_found, y_found, color))
        xs, ys, colors = zip(*markers_data)
        marker.set_offsets(np.column_stack([xs, ys]))
        marker.set_color(colors)
        fig.canvas.draw_idle()

        # Plot Z-line on the right subplot
        ax_z.plot(z_line['z'], z_line['frequency'], color=color, label=f"(x={x_found:.4f}, y={y_found:.4f})")
        ax_z.legend()
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key.lower() == 'd':
            print("Exporting Z-line spectra to CSV...")
            export_to_csv(ax_z, x_label="Z", y_label="df")

    fig.canvas.mpl_connect('scroll_event', on_scroll)
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    plt.tight_layout()
    plt.show()

# === Main Execution ===
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
xsf_files = [f for f in os.listdir(script_dir) if f=="LCPD_volts.xsf"]
if not xsf_files:
    raise FileNotFoundError("File LCPD_volts.xsf not found in the current folder.")
xsf_filename = os.path.join(script_dir, xsf_files[0])
print(f"Using file: {xsf_filename}")

voxel_vector, lattice_vectors, df_table = load_or_parse_dataframe(xsf_filename)
interactive_xy_plane_scroll(df_table, initial_z=11.16, n_colors=10)

