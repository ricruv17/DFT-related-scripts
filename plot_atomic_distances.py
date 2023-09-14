"""
Script to graph the geometric structure of a molecule from a xyz file, while representing the bond length with different
 colors.
Written by Ricardo Ruvalcaba at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
Version 1 (28/03/23) known bugs:
    - When plotting the bonds between all atoms except H, bonds further away from the real bonding distances are
    detected in metallic atoms. Those bonds are being hidden, but the colorbar is scaled incorrectly.
Version 2 (06/04/23). No known bugs.
Version 3 (04/05/23).
    - Added capability to read file name from command line input.
    Known bugs:
    - Fails to plot structures with more than 162 atoms. More extensive testing and fixing is needed.
Version 4 (14/09/23).
    - Added capability to alignt the symmetry axis of a molecule with the x-axis.
    Known bugs:
    - Fails to plot structures with more than 162 atoms. More extensive testing and fixing is needed.

Atomic radii and colors taken from xyz2graph program (https://github.com/zotko/xyz2graph).

When executed, the script will ask for all the necessary information.
To execute, just run on your bash terminal:
python3 plot_atomic_distances.py FILENAME.xyz

Or the default:
python3 plot_atomic_distances.py  # ->  will seek for the file "geometry.out.xyz",
                                         otherwise will ask for a correct file name.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

bond_constant = 2.2
font_size = 16
angle_criteria = 10 * np.pi/180  # in radians

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
        self.interatomic_angles = []
        self.is_bond_plotted_matrix = []
        self.min_bond_distance = 100
        self.max_bond_distance = 0
        self.no_atoms = 0
        self.draw_only_carbon_bonds = False

    def find_file(self):
        # Searches for a file with a given name in the current folder.
        arguments = sys.argv
        if len(arguments) > 1:
            try:
                filename = sys.argv[1]
                open(filename)
                print(f'\nFile {filename} found in the current folder.\nOpening now...')
                if filename[-4:] != '.xyz':
                    print(f'\nFile {filename} does not have an xyz extension.\nProgram may fail...\n')
                return filename
            except FileNotFoundError:
                print(f'\nERROR: File {filename} was not found in the current folder.')
                sys.argv[1] = input('Please enter the name of a valid xyz file: ')
                return self.find_file()
        else:
            try:
                open('geometry.out.xyz')
                print(f'\nFile geometry.out.xyz found in the current folder.\nOpening now...')
                return 'geometry.out.xyz'
            except FileNotFoundError:
                print(f'\nFile geometry.out.xyz was not found in the current folder.')
                sys.argv.append(input('Please enter the name of your xyz file: '))
                return self.find_file()

    def read_xyz_file(self):
        # Reads an XYZ file, searches for elements and their cartesian coordinates and adds them to corresponding arrays
        filename = self.find_file()
        with open(filename) as file_in:
            lines = []
            for line in file_in:
                lines.append(line)
        lines = [line for line in lines if line != '\n']
        self.no_atoms = int(lines[0])
        for line in lines[-self.no_atoms:]:
            line = line.split()
            self.elements.append(line[0])
            self.x.append(float(line[1]))
            self.y.append(float(line[2]))
            self.z.append(float(line[3]))
        self.atomic_radii = [atomic_radii[element] for element in self.elements]

    def translate_molecules_centroid_to_origin(self):
        # Name is self-explanatory.
        centroid_x = sum(self.x)/self.no_atoms
        centroid_y = sum(self.y)/self.no_atoms
        centroid_z = sum(self.z)/self.no_atoms
        for index in range(self.no_atoms):
            self.x[index] -= centroid_x
            self.y[index] -= centroid_y
            self.z[index] -= centroid_z

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
        # 0. Convert the points to a numpy array for easier calculations
        points = []
        for index in range(self.no_atoms):
            points.append((self.x[index], self.y[index]))
        points = np.array(points)

        # 1. Perform Singular Value Decomposition (SVD) on the points
        _, _, vh = np.linalg.svd(points)

        # 2. Extract the right singular vector corresponding to the smallest singular value
        symmetry_axis_vector = vh[0]

        # 3. Calculate the angle of the symmetry axis with the x-axis
        angle = -np.arctan2(symmetry_axis_vector[1], symmetry_axis_vector[0])

        # 4. Rotate and save new positions.
        row1 = [np.cos(angle), -np.sin(angle)]
        row2 = [np.sin(angle), np.cos(angle)]
        rot_matrix = np.array([row1, row2])
        for index in range(self.no_atoms):
            r0 = np.array([self.x[index], self.y[index]])
            rf = np.matmul(rot_matrix, r0)
            self.x[index], self.y[index] = rf

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

    def calculate_interatomic_distances_and_angles(self):
        # Generates matrices with distances, angles and bonding availability between each atom as its elements.
        for atom1_index in range(self.no_atoms):
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            z1 = self.z[atom1_index]
            atom1_distances = []
            atom1_angles = []
            for atom2_index in range(self.no_atoms):
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                z2 = self.z[atom2_index]
                distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                atom1_distances.append(distance)
                bond_angle = self.determine_bond_angle(x1, y1, x2, y2)
                atom1_angles.append(bond_angle)
            self.interatomic_distances.append(atom1_distances)
            self.interatomic_angles.append(atom1_angles)
            self.is_bond_plotted_matrix.append([1] * len(atom1_distances))

    @staticmethod
    def determine_bond_angle(x1, y1, x2, y2):
        # Name is self-explanatory.
        delta_x = x2 - x1
        delta_y = y2 - y1
        if delta_y == 0 and delta_x == 0:
            return 361
        if delta_x == 0 and delta_y > 0:
            return np.pi / 2
        elif delta_x == 0 and delta_y < 0:
            return -np.pi / 2
        angle = np.arctan(delta_y / delta_x)
        if delta_y >= 0:
            if delta_x >= 0:
                angle += 0
            elif delta_x <= 0:
                angle += 2 * np.pi / 2
        elif delta_y <= 0:
            if delta_x <= 0:
                angle += 2 * np.pi / 2
            elif delta_x >= 0:
                angle += 4 * np.pi / 2
        return angle

    def exclude_nonbonding_atoms(self):
        # Name is self-explanatory. Works based on several criteria, further commented in the function.
        bonding_distances = set()
        for atom1_index in range(self.no_atoms):
            atom1_element = self.elements[atom1_index]
            atom1_radius = atomic_radii[self.elements[atom1_index]]
            atom1_distances = self.interatomic_distances[atom1_index]
            forbidden_angles = []
            for distance in sorted(atom1_distances):
                atom2_index = atom1_distances.index(distance)
                atom2_element = self.elements[atom2_index]
                atom2_radius = atomic_radii[self.elements[atom2_index]]
                # 1. Exclude double counting a bond and self-bond.
                if atom1_index >= atom2_index:
                    self.is_bond_plotted_matrix[atom1_index][atom2_index] = 0
                else:
                    self.is_bond_plotted_matrix[atom1_index][atom2_index] = 1
                # 2. Exclude bonds to all atoms but carbon or only hydrogen depending on the user's choice.
                both_atoms_are_carbons = atom1_element == 'C' and atom2_element == 'C'
                one_atom_is_hydrogen = atom1_element == 'H' or atom2_element == 'H'
                if self.draw_only_carbon_bonds:
                    self.is_bond_plotted_matrix[atom1_index][atom2_index] *= both_atoms_are_carbons
                else:
                    self.is_bond_plotted_matrix[atom1_index][atom2_index] *= not one_atom_is_hydrogen
                # 3. Exclude bonds that are too long.
                max_radius = max(atom1_radius, atom2_radius)
                atom_is_within_allowed_radius = distance <= max_radius*bond_constant
                self.is_bond_plotted_matrix[atom1_index][atom2_index] *= atom_is_within_allowed_radius
                # 4. Exclude bonds that have the allowed distance but are in the same direction as a known bond.
                angle = self.interatomic_angles[atom1_index][atom2_index]
                is_allowed_angle = self.determine_if_angle_is_allowed(angle, forbidden_angles)
                self.is_bond_plotted_matrix[atom1_index][atom2_index] *= is_allowed_angle
                if is_allowed_angle:
                    forbidden_angles.append([angle - angle_criteria / 2, angle + angle_criteria / 2])
                bonding_distances.add(
                    self.is_bond_plotted_matrix[atom1_index][atom2_index] * self.interatomic_distances[atom1_index][
                        atom2_index])
        bonding_distances = sorted(list(bonding_distances))
        self.min_bond_distance = bonding_distances[1]
        self.max_bond_distance = bonding_distances[-1]

    @staticmethod
    def determine_if_angle_is_allowed(angle, forbidden_angles):
        # Name is self-explanatory.
        for nth_range in forbidden_angles:
            if nth_range[0] < angle < nth_range[1]:
                return False
        return True


    def plot_molecule(self):
        # Name is self-explanatory. Sections of the function will be explained.
        # 0. Create the figure and colormap objects and position the molecule for correct visualization.
        fig = plt.figure(figsize=[6, 6])
        ax = fig.add_subplot(1, 1, 1)
        cmap = plt.get_cmap("coolwarm")
        norm = mpl.colors.Normalize(self.min_bond_distance, self.max_bond_distance)
        # 1. Plot the atoms.
        for atom in range(self.no_atoms):
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
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            for atom2_index in range(atom1_index + 1, self.no_atoms):
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                distance = self.interatomic_distances[atom1_index][atom2_index]
                if self.is_bond_plotted_matrix[atom1_index][atom2_index]:
                    ax.plot([x1, x2], [y1, y2], linewidth=10 * atomic_radii['C'], color=cmap(norm(distance)), zorder=-1)
        # 3. Add the color-bar and rest of details to the plot.
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        ticks = [round(a, 3) for a in np.linspace(1.001 * self.min_bond_distance, 0.999 * self.max_bond_distance, 5)]
        cbar = fig.colorbar(sm, ax=ax, shrink=0.7, orientation='horizontal', extend='both', ticks=ticks)
        cbar.set_label(r'Bond length ($\AA$)', fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim([min(self.x) - 0.5, max(self.x) + 0.5])
        ax.set_ylim([min(self.y) - 0.5, max(self.y) + 0.5])
        plt.tight_layout()
        plt.show()


# Generate the molecule object and call it when the script is run.
# Default file is "geometry.out.xyz". Script will ask for another file if not found.
molecule = MoleculeGraph()
molecule.read_xyz_file()
molecule.translate_molecules_centroid_to_origin()
molecule.rotate_molecule_towards_z_axis()
molecule.align_molecule_with_x_axis()
molecule.ask_type_of_bonds_drawn()
molecule.calculate_interatomic_distances_and_angles()
molecule.exclude_nonbonding_atoms()
molecule.plot_molecule()
