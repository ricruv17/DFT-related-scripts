"""
Script that applies a spherical Gaussian filter to molecular orbitals in CUBE file format
and writes the result into a new CUBE file.

Written by Ricardo Ruvalcaba at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
This is version 1 (01/06/23). No known bugs. Features yet to implement:
    - It only supports CUBE files with cubic axis vectors. Add capability to read
    other types of coordinate systems.
    - Plot and format final images of the orbitals (most likely on VMD).
    - Let user choose surface color.
    - Print out energies, occupation, etc.
    - Rotate and align the molecule.

CUBE file reading and writing code adapted from:
https://gist.github.com/aditya95sriram/8d1fccbb91dae93c4edf31cd6a22510f

When executed, the script will ask for all the necessary information.
To execute, just run on your bash terminal:
python3 plot_molecular_orbital_surfaces.py FILENAME.cube

Or the default:
python3 plot_molecular_orbital_surfaces.py  # ->  the program  will ask for a correct file name.
"""
import numpy as np
from scipy import signal
import sys


class OrbitalSurfaces:
    """Represents a graph of a molecule's particular orbital given in a CUBE file."""
    def __init__(self):
        self.filename = None
        self.no_atoms = 0
        self.origin_coordinates = []
        self.no_voxels = []
        self.axis_vectors = []
        self.elements = []
        self.charges = []
        self.atoms_x_coordinate = []
        self.atoms_y_coordinate = []
        self.atoms_z_coordinate = []
        self.volume_data = None
        self.x_axis_mesh = None
        self.y_axis_mesh = None
        self.z_axis_mesh = None

    def find_file(self):
        """Searches for a file with a given name in the current folder, initially taken
        from the bash prompt.
        If the file is not found, the function will ask the user to manually type the location
        of the file.
        """
        arguments = sys.argv
        # If the user enters a file name in the prompt.
        if len(arguments) > 1:
            try:
                filename = sys.argv[1]
                open(filename)
                print(f'\nFile {filename} found in the current folder.\nOpening now...')
                # Warn the user in case their file does not have the correct extension.
                if not filename.endswith('cub') or not filename.endswith('cube'):
                    print(f'\nFile {filename} does not have a CUBE extension.\nProgram may fail...\n')
                return filename
            except FileNotFoundError: # In case the file is not found.
                print(f'\nERROR: File {filename} was not found in the current folder.')
                sys.argv[1] = input('Please enter the name of a valid CUBE file: ')
                return self.find_file()
        # If the user does not enter a file name in the prompt.
        else:
            sys.argv.append(input('Please enter the name of your CUBE file: '))
            return self.find_file()

    def read_cube_file(self):
        """Reads a CUBE file and stores its info in the object."""
        self.filename = self.find_file() # Gets the name of the desired file.
        # Begin reading of the CUBE file.
        with open(self.filename) as file_in:
            # Ignore comments.
            file_in.readline()
            file_in.readline()
            # Extract unit cell data.
            line = file_in.readline()
            self.no_atoms = int(line.split()[0])
            self.origin_coordinates = list(map(float, line.split()[1:]))
            for line_number in [4, 5, 6]:
                line = file_in.readline()
                self.no_voxels.append(int(line.split()[0]))
                self.axis_vectors.append(list(map(float, line.split()[1:])))
            # Extract voxels' length and convert to Angstroms if given in Bohrs.
            voxel_length_x_axis = self.axis_vectors[0][0] * 0.529177
            voxel_length_y_axis = self.axis_vectors[1][1] * 0.529177
            voxel_length_z_axis = self.axis_vectors[2][2] * 0.529177
            if voxel_length_x_axis < 0:
                voxel_length_x_axis = -voxel_length_x_axis / 0.529177
                voxel_length_y_axis = -voxel_length_y_axis / 0.529177
                voxel_length_z_axis = -voxel_length_z_axis / 0.529177
            # Extract atoms' data.
            for line_number in range(6, 6 + self.no_atoms):
                line = file_in.readline()
                line = line.split()
                self.elements.append(int(line[0]))
                self.charges.append(float(line[1]))
                self.atoms_x_coordinate.append(float(line[2]))
                self.atoms_y_coordinate.append(float(line[3]))
                self.atoms_z_coordinate.append(float(line[4]))
            # Extract volumetric data.
            self.volume_data = np.zeros((np.prod(self.no_voxels)))
            line_no = 0
            for line in file_in:
                for val in line.strip().split():
                    self.volume_data[line_no] = float(val)
                    line_no += 1
        self.volume_data = np.reshape(self.volume_data, self.no_voxels)
        r0 = -np.array([self.no_voxels[0] * voxel_length_x_axis, self.no_voxels[1] * voxel_length_y_axis, self.no_voxels[2] * voxel_length_z_axis])/2
        rf = -r0
        mesh = np.mgrid[r0[0]:rf[0]:complex(0, self.no_voxels[0]), r0[1]:rf[1]:complex(0, self.no_voxels[1]), r0[2]:rf[2]:complex(0, self.no_voxels[2])]
        self.x_axis_mesh, self.y_axis_mesh, self.z_axis_mesh = mesh

    def write_cube_file(self, filename):
        """Writes the information of the object into a CUBE file with the correct formatting.
        
        Parameters:
            filename: name of the newly written CUBE file. 
        """
        with open(filename, 'w') as file_out:
            file_out.write(f' Gaussian-smoothed version of {self.filename}.\n *********************************************\n') # Write the comment lines.
            # Write unit cell data.
            new_no_voxels = np.shape(self.volume_data)
            file_out.write("{:5} {:11.6f} {:11.6f} {:11.6f}\n".format(self.no_atoms, self.origin_coordinates[0], self.origin_coordinates[1], self.origin_coordinates[2]))
            for axis in [0, 1, 2]:
                file_out.write("{:5} {:11.6f} {:11.6f} {:11.6f}\n".format(new_no_voxels[axis], self.axis_vectors[axis][0], self.axis_vectors[axis][1], self.axis_vectors[axis][2]))
            # Write atoms' data.
            for atom in range(self.no_atoms):
                file_out.write("{:5} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(self.elements[atom], self.charges[atom], self.atoms_x_coordinate[atom], self.atoms_y_coordinate[atom], self.atoms_z_coordinate[atom]))
            # Format and write volumetric data.
            for i in range(new_no_voxels[0]):
                for j in range(new_no_voxels[1]):
                    for k in range(new_no_voxels[2]):
                        if (i or j or k) and k%6==0:
                            file_out.write("\n")
                        file_out.write(" {0: .5E}".format(self.volume_data[i,j,k]))


    def apply_gaussian_filter(self, sigma=1):
        """Generates a Gaussian kernel in real space and convolutes it with the volumetric data.
        
        Parameters:
            sigma: standard deviation of the normal distribution.
        """
        kernel = 1/sigma * np.exp(-(self.x_axis_mesh**2 + self.y_axis_mesh**2 + self.z_axis_mesh**2) / (2*sigma**2))
        self.volume_data = signal.convolve(self.volume_data, kernel, mode='same')


molecule = OrbitalSurfaces()
molecule.read_cube_file()
# This value for sigma reproduced correctly the results from the article
# https://doi.org/10.1002/ange.202009200
molecule.apply_gaussian_filter(sigma=2)
molecule.write_cube_file(f'smoothened_{molecule.filename}')
