'''
Script to graph the geometric structure of a molecule from a xyz file, while representing the bond length with different colors.
Written by Ricardo Ruvalcaba at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
This is version 1 (20/03/23)

Class backbone taken from xyz2graph program (https://github.com/zotko/xyz2graph).

When executed, the script will ask for all the information.
To execute, just run on your bash terminal:
python3 FILENAME.out

Things to implement:
- Rotate the molecule perpendicular to the z-direction towards the xyz file
- make the molecule larger and colorbar closer to each other
- make the metal atoms bond only with closest neighbors
- ask user if they want to plot distances only with C-C or all-all

'''

# 0. Import useful packages and define globar variables and functions
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import re

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
    other = "pink"
)



class MoleculeGraph:
    "Represents a graph of a molecule with bonds colored as function of length."

    def __init__(self):
        self.elements = []
        self.x = []
        self.y = []
        self.z = []
        self.atomic_radii = []
        self.interatomic_distances = []
        self.min_bond_distance = 100
        self.max_bond_distance = 0


    def read_xyz_file(self, filename):
        "Reads an XYZ file, searches for elements and their cartesian coordinates and adds them to corresponding arrays."

        error_message = f'File {filename} was not found in the current folder.'
        try:
            pattern = re.compile(
                r"([A-Za-z]{1,3})\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)"
            )
            with open(filename) as file_in:
                for element, x, y, z in pattern.findall(file_in.read()):
                    self.elements.append(element)
                    self.x.append(float(x))
                    self.y.append(float(y))
                    self.z.append(float(z))

            self.atomic_radii = [atomic_radii[element] for element in self.elements]
            self.calculate_interatomic_distances()
            if filename == 'final_positions.xyz':
                print(f'\nFile {filename} found in the current folder.\nOpening now...')
            elif filename[-4:] != '.xyz':
                print(f'\nFile {filename} is not an xyz file.')
                return self.read_xyz_file(input('Please enter the name of a valid xyz file:\n'))

        except FileNotFoundError:
            if filename == 'final_positions.xyz':
                print('\n' + error_message)
            else:
                print('\nERROR: ' + error_message)
            return self.read_xyz_file(input('Please enter the name of a valid xyz file:\n'))


    def calculate_interatomic_distances(self):
        no_atoms = len(self.x)
        for atom1_index in range(no_atoms):
            atom1_element = self.elements[atom1_index]
            atom1_radius = atomic_radii[self.elements[atom1_index]]
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            z1 = self.z[atom1_index]
            atom1_distances = []
            for atom2_index in range(no_atoms):
                atom2_element = self.elements[atom2_index]
                atom2_radius = atomic_radii[self.elements[atom2_index]]
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                z2 = self.z[atom2_index]
                distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                atom1_distances.append(distance)
                max_radius = max(atom1_radius, atom2_radius)
                if self.max_bond_distance < distance <= max_radius*bond_constant:
                    self.max_bond_distance = distance
                elif 0 != distance < self.min_bond_distance and atom1_element != 'H' and atom2_element != 'H':
                    self.min_bond_distance = distance
            self.interatomic_distances.append(atom1_distances)
        

    def plot_molecule(self):
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.04])
        fig = plt.figure(figsize=(6, 7))
        ax = fig.add_subplot(gs[0, 0], projection='3d')
        cax = fig.add_subplot(gs[1, 0])
        cmap = plt.get_cmap("coolwarm")
        norm = mpl.colors.Normalize(self.min_bond_distance, self.max_bond_distance)
        # scatter the atoms
        for atom in range(len(self.x)):
            x = self.x[atom]
            y = self.y[atom]
            z = self.z[atom]
            size = atomic_radii[self.elements[atom]]
            if self.elements[atom] in atomic_colors:
                color = atomic_colors[self.elements[atom]]
            else:
                color = atomic_colors['other']
            ax.scatter(x, y, z, s=100*size, c=color)
        
        # plot the bonds between atoms with colors depending on length and excluding Hydrogen
        no_atoms = len(self.x)
        for atom1_index in range(no_atoms):
            atom1_radius = atomic_radii[self.elements[atom1_index]]
            atom1_element = self.elements[atom1_index]
            x1 = self.x[atom1_index]
            y1 = self.y[atom1_index]
            z1 = self.z[atom1_index]
            for atom2_index in range(atom1_index, no_atoms):
                atom2_radius = atomic_radii[self.elements[atom2_index]]
                atom2_element = self.elements[atom2_index]
                x2 = self.x[atom2_index]
                y2 = self.y[atom2_index]
                z2 = self.z[atom2_index]
                distance = self.interatomic_distances[atom1_index][atom2_index]
                min_radius = min(atom1_radius, atom2_radius)
                max_radius = max(atom1_radius, atom2_radius)
                if distance <= max_radius*bond_constant and atom1_element != 'H' and atom2_element != 'H':
                    ax.plot([x1, x2], [y1, y2], [z1, z2], linewidth=8*min_radius, color=cmap(norm(distance)), zorder=-1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        ticks = [round(a, 3) for a in np.linspace(self.min_bond_distance, self.max_bond_distance, 5)]
        cbar = fig.colorbar(sm, cax=cax, shrink=0.7, orientation='horizontal', extend='both', ticks=ticks)
#        cbar.minorticks_on()
        cbar.set_label(r'Bond lengths ($\AA$)', fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)
        ax.set_aspect('equal')
        ax.view_init(90, 270)
        ax.axis('off')
        ax.set_xlim([min(self.x), max(self.x)])
        ax.set_ylim([min(self.y), max(self.y)])
        plt.tight_layout()
        plt.show()
        


molecule = MoleculeGraph()
#molecule.read_xyz_file('gas_phase_b86pbe25_pdcomplex.xyz')
molecule.read_xyz_file('anthr.xyz')
molecule.plot_molecule()


