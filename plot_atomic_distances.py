"""
Script to graph the geometric structure of a molecule from a xyz file, while representing the bond length with different
 colors.
Written by Ricardo Ruvalcaba at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
This is version 1 (28/03/23)
    Issues to fix:
    - When plotting the bonds between all atoms except H, bonds further away from the real bonding distances are
    detected in metallic atoms. Those bonds are being hidden, but the colorbar is scaled incorrectly.

Atomic radii and colors taken from xyz2graph program (https://github.com/zotko/xyz2graph).

When executed, the script will ask for all the necessary information.
To execute, just run on your bash terminal:
python3 FILENAME.xyz
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

bond_constant = 2.2
font_size = 16

atomic_radii = dict(
    Ac=1.88,
    Ag=1.59,
    Al=1.35,
    Am=1.51,
    As=1.21,
    Au=1.50,
    B=0.83,
    Ba=1.34,
    Be=0.35,
    Bi=1.54,
    Br=1.21,
    C=0.68,
    Ca=0.99,
    Cd=1.69,
    Ce=1.83,
    Cl=0.99,
    Co=1.33,
    Cr=1.35,
    Cs=1.67,
    Cu=1.52,
    D=0.23,
    Dy=1.75,
    Er=1.73,
    Eu=1.99,
    F=0.64,
    Fe=1.34,
    Ga=1.22,
    Gd=1.79,
    Ge=1.17,
    H=0.23,
    Hf=1.57,
    Hg=1.70,
    Ho=1.74,
    I=1.40,
    In=1.63,
    Ir=1.32,
    K=1.33,
    La=1.87,
    Li=0.68,
    Lu=1.72,
    Mg=1.10,
    Mn=1.35,
    Mo=1.47,
    N=0.68,
    Na=0.97,
    Nb=1.48,
    Nd=1.81,
    Ni=1.50,
    Np=1.55,
    O=0.68,
    Os=1.37,
    P=1.05,
    Pa=1.61,
    Pb=1.54,
    Pd=1.50,
    Pm=1.80,
    Po=1.68,
    Pr=1.82,
    Pt=1.50,
    Pu=1.53,
    Ra=1.90,
    Rb=1.47,
    Re=1.35,
    Rh=1.45,
    Ru=1.40,
    S=1.02,
    Sb=1.46,
    Sc=1.44,
    Se=1.22,
    Si=1.20,
    Sm=1.80,
    Sn=1.46,
    Sr=1.12,
    Ta=1.43,
    Tb=1.76,
    Tc=1.35,
    Te=1.47,
    Th=1.79,
    Ti=1.47,
    Tl=1.55,
    Tm=1.72,
    U=1.58,
    V=1.33,
    W=1.37,
    Y=1.78,
    Yb=1.94,
    Zn=1.45,
    Zr=1.56
)


atomic_colors = dict(
    Ar="cyan",
    B="salmon",
    Ba="darkgreen",
    Be="darkgreen",
    Br="darkred",
    C="black",
    Ca="darkgreen",
    Cl="green",
    Cs="violet",
    F="green",
    Fe="darkorange",
    Fr="violet",
    H="silver",
    He="cyan",
    I="darkviolet",
    K="violet",
    Kr="cyan",
    Li="violet",
    Mg="darkgreen",
    N="blue",
    Na="violet",
    Ne="cyan",
    O="red",
    P="orange",
    Ra="darkgreen",
    Rb="violet",
    S="yellow",
    Sr="darkgreen",
    Ti="gray",
    Xe="cyan",
    other="palevioletred"
)


class MoleculeGraph:
    # Represents a graph of a molecule with bonds colored as function of length.
    def __init__(self):
        self.elements = []
        self.x = []
        self.y = []
        self.z = []
        self.atomic_radii = []
        self.interatomic_distances = []
        self.min_bond_distance = 100
        self.max_bond_distance = 0
        self.no_atoms = 0
        self.draw_only_carbon_bonds = False

    def read_xyz_file(self, filename):
        # Reads an XYZ file, searches for elements and their cartesian coordinates and adds them to corresponding arrays
        error_message = f'File {filename} was not found in the current folder.'
        try:
            with open(filename) as file_in:
                lines = []
                for line in file_in:
                    lines.append(line)
            lines = [line for line in lines if line != '\n']
            if filename == 'final_positions.xyz':
                print(f'\nFile {filename} found in the current folder.\nOpening now...')
            elif filename[-4:] != '.xyz':
                print(f'\nFile {filename} is not an xyz file.')
                return self.read_xyz_file(input('Please enter the name of a valid xyz file:\n'))
            self.no_atoms = int(lines[0])
            for line in lines[-self.no_atoms:]:
                line = line.split()
                self.elements.append(line[0])
                self.x.append(float(line[1]))
                self.y.append(float(line[2]))
                self.z.append(float(line[3]))
            self.atomic_radii = [atomic_radii[element] for element in self.elements]
            self.ask_type_of_bonds_drawn()
            self.calculate_interatomic_distances()
        except FileNotFoundError:
            if filename == 'final_positions.xyz':
                print('\n' + error_message)
            else:
                print('\nERROR: ' + error_message)
            return self.read_xyz_file(input('Please enter the name of a valid xyz file:\n'))

    def ask_type_of_bonds_drawn(self):
        # Asks the user if only C-C or all types of bonds (except H) will be plotted.
        selection = input('\nSelect the type of bonds you want to plot: (1) C-C only or (2) all except H: ')
        while selection != '1' and selection != '2':
            print('\nERROR: Only values 1 or 2 are allowed.')
            selection = input('Please enter a valid number: ')
        if selection == '1':
            self.draw_only_carbon_bonds = True
        elif selection == '2':
            self.draw_only_carbon_bonds = False

    def calculate_interatomic_distances(self):
        # Generates a matrix with distances between each atom as its elements.
        for atom1_index in range(self.no_atoms):
            atom1_element = self.elements[atom1_index]
            atom1_radius = atomic_radii[self.elements[atom1_index]]
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            z1 = self.z[atom1_index]
            atom1_distances = []
            for atom2_index in range(self.no_atoms):
                atom2_element = self.elements[atom2_index]
                atom2_radius = atomic_radii[self.elements[atom2_index]]
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                z2 = self.z[atom2_index]
                distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                atom1_distances.append(distance)
                max_radius = max(atom1_radius, atom2_radius)
                if self.draw_only_carbon_bonds:
                    both_atoms_are_carbons = atom1_element == 'C' and atom2_element == 'C'
                    if self.max_bond_distance < distance <= max_radius*bond_constant and both_atoms_are_carbons:
                        self.max_bond_distance = distance
                    elif 0 != distance < self.min_bond_distance and both_atoms_are_carbons:
                        self.min_bond_distance = distance
                else:
                    if self.max_bond_distance < distance <= max_radius*bond_constant:
                        self.max_bond_distance = distance
                    elif 0 != distance < self.min_bond_distance and atom1_element != 'H' and atom2_element != 'H':
                        self.min_bond_distance = distance
            self.interatomic_distances.append(atom1_distances)

    def translate_molecules_geometric_center_to_origin(self):
        # Name is self-explanatory.
        geometric_center_x = sum(self.x)/self.no_atoms
        geometric_center_y = sum(self.y)/self.no_atoms
        geometric_center_z = sum(self.z)/self.no_atoms
        for index in range(self.no_atoms):
            self.x[index] -= geometric_center_x
            self.y[index] -= geometric_center_y
            self.z[index] -= geometric_center_z

    def rotate_molecule_towards_z_axis(self):
        # Name is self-explanatory. Sections of the function will be explained.
        # 0. Find plane that fits best all the atoms in the molecule.
        min_squares_matrix = np.array([self.x, self.y]).T
        n = np.dot(np.matmul(np.linalg.inv(np.matmul(min_squares_matrix.T, min_squares_matrix)), min_squares_matrix.T),
                   np.array(self.z))
        nz = (n[0]**2 + n[1]**2 + 1)**-0.5
        nx = -n[0] * nz
        ny = -n[1] * nz
        # 1. Find rotation axis and angle between that plane and the z-axis.
        ux, uy, uz = np.cross([nx, ny, nz], [0, 0, 1])
        theta = np.arccos(nz)
        row1 = [np.cos(theta) + ux ** 2 * (1 - np.cos(theta)),
                ux * uy * (1 - np.cos(theta)) - uz * np.sin(theta),
                ux * uz * (1 - np.cos(theta)) + uy * np.sin(theta)]
        row2 = [uy * ux * (1 - np.cos(theta)) + uz * np.sin(theta),
                np.cos(theta) + uy ** 2 * (1 - np.cos(theta)),
                uy * uz * (1 - np.cos(theta)) - ux * np.sin(theta)]
        row3 = [uz * ux * (1 - np.cos(theta)) - uy * np.sin(theta),
                uz * uy * (1 - np.cos(theta)) + ux * np.sin(theta),
                np.cos(theta) + uz ** 2 * (1 - np.cos(theta))]
        rot_matrix = np.array([row1, row2, row3])
        # 2. Rotate and save new positions
        for index in range(self.no_atoms):
            r0 = np.array([self.x[index], self.y[index], self.z[index]])
            rf = np.matmul(rot_matrix, r0)
            self.x[index], self.y[index], self.z[index] = rf

    def align_molecule_with_x_axis(self):
        # Name is self-explanatory. Sections of the function will be explained.
        # 0. Find the straight line that fits best the atoms in the xy-plane and its angle-difference with the x-axis.
        min_squares_matrix = np.array([[1] * self.no_atoms, self.x]).T
        n = np.dot(np.matmul(np.linalg.inv(np.matmul(min_squares_matrix.T, min_squares_matrix)), min_squares_matrix.T),
                   np.array(self.y))
        u = np.array([1, n[1]]) - np.array([0, n[0]])
        u = u/np.linalg.norm(u)
        optimal_theta = np.arccos(u[0])
        # 1. Rotate and save new positions.
        row1 = [np.cos(optimal_theta), -np.sin(optimal_theta)]
        row2 = [np.sin(optimal_theta), np.cos(optimal_theta)]
        rot_matrix = np.array([row1, row2])
        for index in range(self.no_atoms):
            r0 = np.array([self.x[index], self.y[index]])
            rf = np.matmul(rot_matrix, r0)
            self.x[index], self.y[index] = rf

    def plot_molecule(self):
        # Name is self-explanatory. Sections of the function will be explained.
        # 0. Create the figure and colormap objects and position the molecule for correct visualization.
        fig = plt.figure(figsize=[int(max(self.x) * 1.5), int(max(self.y) * 2)])
        ax = fig.add_subplot(1, 1, 1)
        cmap = plt.get_cmap("coolwarm")
        norm = mpl.colors.Normalize(self.min_bond_distance, self.max_bond_distance)
        self.translate_molecules_geometric_center_to_origin()
        self.rotate_molecule_towards_z_axis()
        self.align_molecule_with_x_axis()
        # 1. Plot the atoms.
        for atom in range(len(self.x)):
            x = self.x[atom]
            y = self.y[atom]
            size = atomic_radii[self.elements[atom]]
            if self.elements[atom] in atomic_colors:
                color = atomic_colors[self.elements[atom]]
            else:
                color = atomic_colors['other']
            ax.scatter(x, y, s=100*size, c=color)
        # 2. Plot the bonds between atoms with colors depending on length according to user's choice.
        for atom1_index in range(self.no_atoms):
            atom1_radius = atomic_radii[self.elements[atom1_index]]
            atom1_element = self.elements[atom1_index]
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            for atom2_index in range(atom1_index, self.no_atoms):
                atom2_radius = atomic_radii[self.elements[atom2_index]]
                atom2_element = self.elements[atom2_index]
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                distance = self.interatomic_distances[atom1_index][atom2_index]
                max_radius = max(atom1_radius, atom2_radius)
                if self.draw_only_carbon_bonds:
                    both_atoms_are_carbons = atom1_element == 'C' and atom2_element == 'C'
                    if distance <= max_radius * bond_constant and both_atoms_are_carbons:
                        ax.plot([x1, x2], [y1, y2], linewidth=10 * atomic_radii['C'], color=cmap(norm(distance)),
                                zorder=-distance)
                else:
                    if distance <= max_radius*bond_constant and atom1_element != 'H' and atom2_element != 'H':
                        ax.plot([x1, x2], [y1, y2], linewidth=10*atomic_radii['C'], color=cmap(norm(distance)),
                                zorder=-distance)
        # 3. Add the colorbar and rest of details to the plot.
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        ticks = [round(a, 3) for a in np.linspace(1.001 * self.min_bond_distance, 0.999 * self.max_bond_distance, 5)]
        cbar = fig.colorbar(sm, ax=ax, shrink=0.7, orientation='horizontal', extend='both', ticks=ticks)
        cbar.set_label(r'Bond lengths ($\AA$)', fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim([min(self.x) - 0.5, max(self.x) + 0.5])
        ax.set_ylim([min(self.y) - 0.5, max(self.y) + 0.5])
        plt.tight_layout()
        plt.show()


# Generate the molecule object and call it when the script is run.
# Default file is "final_positions.xyz". Script will ask for another file if not found.
molecule = MoleculeGraph()
molecule.read_xyz_file('final_positions.xyz')
molecule.plot_molecule()
