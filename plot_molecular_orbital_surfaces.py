"""
Script that applies a spherical Gaussian filter to all the CUBE files in the current folder that come from 
FHI-aims eigenstate calculations, writes the result into new CUBE files, plots them on VMD, and stitches the
images together.

It only works with files that have standard names output by FHI-aims. Examples:
aims.out                                    Main output file of the calculation
cube_001_eigenstate_00117_spin_1.cube       CUBE file containing one Kohn-Sham eigenstate

Written by Ricardo Ruvalcaba and Shadi Fatayer at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
Version 1 (01/06/23). No known bugs. Features yet to implement:
    - Plot and format final images of the orbitals (most likely on VMD).
    - Ask user if he wants to plot both the original and the filtered orbitals.
    - Let user choose surface color.
    - Rotate and align the molecule.
    - Print filenames of the images that are being generated one by one.
    - It only supports CUBE files with cubic axis vectors. Add capability to read
    other types of coordinate systems.
This is version 2 (30/06/23). Features yet to implement:
    - Ask user if he wants to plot both the original and the filtered orbitals.
    - Let user choose surface color.
    - Print filenames of the images that are being generated one by one.
    - Make the program faster.
    - Orient the molecule along the x-axis
    - Fixing bug that squishes the images when the amout of cube files is smaller than 3
    - Fix bug with Euler's angles

CUBE file reading and writing code adapted from:
https://gist.github.com/aditya95sriram/8d1fccbb91dae93c4edf31cd6a22510f

When executed, the script will automatically detect all CUBE files with the correct nomenclature.
To execute, just run on your bash terminal:
python3 plot_molecular_orbital_surfaces.py
"""
from __future__ import print_function
import argparse
import sys
import re
import subprocess
import os
import datetime
from os import listdir, environ
from os.path import isfile, join
from future.utils import iteritems

import numpy as np
import os
import scipy

import matplotlib.pyplot as plt
import numpy as np
import os
import PIL

##############################################################################
"""This part is for the gaussian filtering."""

RESOLUTION = 2048 # 2048 is the maximum
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
        self.aimsout_list = []
        self.arguments = sys.argv.copy()


    def get_sigma_and_isovalues(self):
        """ Forces the user to give correct values for gaussian sigma and surface isovalue. """
        try:
            sigma = float(self.arguments[1])
            isovalue = float(self.arguments[2])
            print(f'\nThe orbitals will be processed with sigma = {sigma}, isovalue = {isovalue}.')
            return sigma, isovalue
        except IndexError:
            print(f'\nYou did not specify the sigma nor the isovalue.')
            self.arguments.append(input('Please enter the value for the gaussian sigma: '))
            self.arguments.append(input('Please enter the isovalue: '))
            return self.get_sigma_and_isovalues()
        except ValueError or TypeError:
            print(f'\nERROR: sigma and the isovalue must be real numbers.')
            self.arguments[1] = input('Please enter a correct value for sigma: ')
            self.arguments[2] = input('Please enter a correct isovalue: ')
            return self.get_sigma_and_isovalues()


    def read_cube_file(self):
        """Reads a CUBE file and stores its info in the object."""
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


    def write_cube_file(self, filename, sigma):
        """Writes the information of the object into a CUBE file with the correct formatting.
        
        Parameters:
            filename: name of the newly written CUBE file. 
        """
        with open(filename, 'w') as file_out:
            text = f'Gaussian-smoothed (sigma = {sigma}) version of {self.filename}'
            file_out.write(text + '\n') # Write the comment lines.
            file_out.write('*' * len(text) + '\n')
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
                        if (i or j or k) and (k % 6 == 0):
                            file_out.write("\n")
                        file_out.write(" {0: .8E}".format(self.volume_data[i, j, k]))


    def apply_gaussian_filter(self, sigma=1):
        """Generates a Gaussian kernel in real space and convolutes it with the volumetric data.
        
        Parameters:
            sigma: standard deviation of the normal distribution.
        """
        # Define the gaussian kernel in reciprocal space.
        shape = (self.no_voxels[0], self.no_voxels[1], self.no_voxels[2])
        center = (shape[0] // 2, shape[1] // 2, shape[2] // 2)
        r = np.ones(shape)
        for i in range(-center[0], center[0]):
            for j in range(-center[1], center[1]):
                for k in range(-center[2], center[2]):
                    u = i / shape[0]
                    v = j / shape[1]
                    w = k / shape[2]
                    r[i][j][k] = np.sqrt( u**2 + v**2 + w**2)
        kernel = np.exp(- 2 * (np.pi * sigma * r)**2)
        # Multiply by the kernel in reciprocal space and go back to real space.
        volume_data_fft = scipy.fftpack.fftn(self.volume_data)
        self.volume_data = scipy.fftpack.ifftn(volume_data_fft * kernel).real


    @staticmethod
    def find_eigenstate_cube_files():
        """ Returns the names of all the CUBE files in the current folder in a list. """
        cube_files = []
        for filename in os.listdir():
            if filename.startswith('cube')\
                and filename.find('eigenstate') != -1\
                    and filename.endswith('cube') or filename.endswith('cub'):
                cube_files.append(filename)
        return cube_files


    def find_which_orbitals_were_calculated(self):
        """ Returns a list with the numbers of the orbitals that were plotted in 
         CUBE files. """
        files = self.find_eigenstate_cube_files()
        calculated_orbitals = {1: [],
                               2: []}
        for filename in files:
            spin = int(filename[31])
            orbital_number = int(filename[20:25])
            calculated_orbitals[spin].append(orbital_number)
        calculated_orbitals[1].sort()
        calculated_orbitals[2].sort()
        return calculated_orbitals


    def print_orbital_data(self):
        """ Prints the state number, occupation, and energy from the aims.out file. """
        # 1. Store the aim.out file in a list for further processing.
        filename = 'aims.out'
        with open(filename) as file_in:
            for line in file_in:
                self.aimsout_list.append(line)
                # 2. Find the numner of KS states in the calculation.
                if 'Number of Kohn-Sham states (occupied + empty):' in line:
                    no_of_KS_states = line.split()
                    no_of_KS_states = int(no_of_KS_states[-1])
        # 3. Determine the position of the final KS states.
        self.aimsout_list.reverse()
        eigenvalues_index = len(self.aimsout_list) - self.aimsout_list.index('  Writing Kohn-Sham eigenvalues.\n') - 1
        self.aimsout_list.reverse()
        # 4. Print to terminal.
        calculated_orbitals = self.find_which_orbitals_were_calculated()
        print('\n  The eigenvalues that were calculated are:')
        calculation_is_spin_polarized = bool(calculated_orbitals[2])
        if not calculation_is_spin_polarized: # print for closed-shell calculations
            print('  State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]')
            initial_orbital = eigenvalues_index + 3 + calculated_orbitals[1][0] - 1
            final_orbital = eigenvalues_index + 3 + calculated_orbitals[1][-1]
            for line in self.aimsout_list[initial_orbital: final_orbital]:
                print(line[:-1])
        else: # print for open-shell calculations
            print('  State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]')
            print('  Spin-up eigenvalues')
            line_lower_limit = eigenvalues_index + 5 + calculated_orbitals[1][0] - 1
            line_upper_limit = eigenvalues_index + 5 + calculated_orbitals[1][-1]
            for line in self.aimsout_list[line_lower_limit: line_upper_limit]:
                print(line[:-1])
            print('  Spin-down eigenvalues')
            line_lower_limit = eigenvalues_index + 5 + no_of_KS_states + 4 + calculated_orbitals[2][0] - 1
            line_upper_limit = eigenvalues_index + 5 + no_of_KS_states + 4 + calculated_orbitals[2][-1]
            for line in self.aimsout_list[line_lower_limit: line_upper_limit]:
                print(line[:-1])

    def calculate_centroid_coordinates(self):
        # Name is self-explanatory.
        centroid_x = sum(self.atoms_x_coordinate)/self.no_atoms
        centroid_y = sum(self.atoms_y_coordinate)/self.no_atoms
        centroid_z = sum(self.atoms_z_coordinate)/self.no_atoms
        return centroid_x, centroid_y, centroid_z

    def get_angles_to_rotate_molecule_towards_z_axis(self):
        """Calculates the angles that will rotate the molecule towards the z axis."""
        # 0. Find plane that fits best all the atoms in the molecule.
        centroid_x, centroid_y, centroid_z = self.calculate_centroid_coordinates()
        x_for_fitting = np.array(self.atoms_x_coordinate) - centroid_x
        y_for_fitting = np.array(self.atoms_y_coordinate) - centroid_y
        z_for_fitting = np.array(self.atoms_z_coordinate) - centroid_z
        min_squares_matrix = np.array([x_for_fitting, y_for_fitting]).T
        n = np.dot(np.matmul(np.linalg.inv(np.matmul(min_squares_matrix.T, min_squares_matrix)), min_squares_matrix.T),
                   np.array(z_for_fitting))
        nz = (n[0]**2 + n[1]**2 + 1)**-0.5
        nx = -n[0] * nz
        ny = -n[1] * nz
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
        matrix = np.array([row1, row2, row3])

        sy = np.sqrt(matrix[0, 0]**2 + matrix[1, 0]**2)
        if sy < 1e-6:
            # Singular case: sy is close to zero
            # This corresponds to a rotation about the y-axis by -90 degrees
            # We can set the other Euler angles to zero in this case
            theta_x = np.arctan2(-matrix[2, 1], matrix[2, 2])
            theta_y = -np.pi/2
            theta_z = 0
        else:
            theta_x = np.arctan2(matrix[2, 1], matrix[2, 2])
            theta_y = np.arctan2(-matrix[2, 0], sy)
            theta_z = np.arctan2(matrix[1, 0], matrix[0, 0])
        return np.array((theta_x, theta_y, theta_z))*180/np.pi


# Generate the new filtered CUBE file
orbital = OrbitalSurfaces()
raw_CUBE_files = OrbitalSurfaces.find_eigenstate_cube_files()
sigma, isovalue = orbital.get_sigma_and_isovalues()
print(f'\n{len(raw_CUBE_files)} raw CUBE files found in the current folder.\nApplying gaussian filter...\n')
file_no = 1
for file in raw_CUBE_files:
    orbital = OrbitalSurfaces()
    orbital.filename = file
    print(f'File {orbital.filename} ({file_no} out of {len(raw_CUBE_files)}) being filtered. Please stand by...')
    orbital.read_cube_file()
    orbital.apply_gaussian_filter(sigma=sigma)
    orbital.write_cube_file(f'filtered_{orbital.filename}', sigma)
    file_no += 1
euler_angles = orbital.get_angles_to_rotate_molecule_towards_z_axis()

##############################################################################
"""This part is for the plotting of the isosurfaces using VMD."""

sys.argv = [sys.argv[0]]

vmd_cube_help = """vmd_cube is a script to render cube files with vmd.
To generate cube files with Psi4 add the command cubeprop() at the end of your input file."""

vmd_exe = ""

vmd_script_name = ".vmd_mo_script.vmd"

vmd_template = """#
# VMD script to plot MOs from cube files
#

# Load the molecule and change the atom style
mol load cube PARAM_CUBEFILE.cube
mol modcolor 0 PARAM_CUBENUM Element
mol modstyle 0 PARAM_CUBENUM Ribbons
#mol modstyle 0 PARAM_CUBENUM Licorice 0.110000 10.000000 10.000000
#mol modstyle 0 PARAM_CUBENUM CPK 0.400000 0.40000 30.000000 16.000000

# Rotate and translate the molecule
rotate x by PARAM_RX
rotate y by PARAM_RY
rotate z by PARAM_RZ
translate by PARAM_TX PARAM_TY PARAM_TZ
scale by PARAM_SCALE

# Eliminate the axis and perfect the view
axes location Off
display projection orthographic
display depthcue off
display resize PARAM_IMAGEW PARAM_IMAGEH
color Display Background white"""


vmd_template_surface = """#
# Add a surface
mol color ColorID PARAM_ISOCOLOR
mol representation Isosurface PARAM_ISOVALUE 0 0 0 1 1
mol selection all
mol material Diffuse
mol addrep PARAM_CUBENUM
"""

vmd_template_interactive = """#
# Disable rendering
mol off PARAM_CUBENUM
"""

vmd_template_render = """
# Render
render TachyonInternal PARAM_CUBEFILE.tga
mol delete PARAM_CUBENUM
"""

default_path = os.getcwd()

# Default parameters
options = {"ISOVALUE"     : [None,"Isosurface Value(s)"],
           "ISOCOLOR"     : [None,"Isosurface Color(s)"],
           "ISOCUT"       : [None,"Isosurface Value Cutoff"],
           "RX"           : [None,"X-axis Rotation"],
           "RY"           : [None,"Y-axis Rotation"],
           "RZ"           : [None,"Z-axis Rotation"],
           "TX"           : [None,"X-axis Translation"],
           "TY"           : [None,"Y-axis Translation"],
           "TZ"           : [None,"Z-axis Translation"],
           "OPACITY"      : [None,"Opacity"],
           "CUBEDIR"      : [None,"Cubefile Directory"],
           "SCALE"        : [None,"Scaling Factor"],
           "MONTAGE"      : [None,"Montage"],
           "LABEL_MOS"    : [None,"Label MOs"],
           "FONTSIZE"     : [None,"Font size"],
           "IMAGEW"       : [None,"Image width"],
           "IMAGEH"       : [None,"Image height"],
           "VMDPATH"      : [None,"VMD Path"],
           "INTERACTIVE"  : [None,"Interactive Mode"],
           "GZIP"         : [None,"Gzip Cube Files"]}


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def multigsub(subs,str):
    for k,v in subs.items():
        str = re.sub(k,v,str)
    return str

environ['VMDPATH'] = '/mnt/c/Program Files \(x86\)/University of Illinois/VMD/vmd.exe'
# environ['VMDPATH'] = '/usr/local/bin/vmd' # for my previous laptop
def find_vmd(options):
    if environ['VMDPATH']:
        vmdpath = environ['VMDPATH']
        vmdpath = multigsub({" " : r"\ "},vmdpath)
        options["VMDPATH"][0] = vmdpath
    else:
        print("Please set the VMDPATH environmental variable to the path of VMD.")
        exit(1)


def save_setup_command(argv):
    file_name = join(default_path, '.vmd_cube_command')
    f = open(file_name, 'w')
    f.write('# setup command was executed '+datetime.datetime.now().strftime("%d-%B-%Y %H:%M:%S"+"\n"))
    f.write(" ".join(argv[:])+"\n")
    f.close()


def read_options(options):
    parser = argparse.ArgumentParser(description=vmd_cube_help)
    parser.add_argument('data', metavar='<cubefile dir>', type=str, nargs='?',default=".",
                   help='The directory containing the cube files.')
                   
    parser.add_argument('--isovalue', metavar='<isovalue>', type=float, nargs='*',default=[isovalue,-isovalue],
                   help='a list of isosurface values (a list of floats, default = [0.05,-0.05])')
    parser.add_argument('--isocolor', metavar='<integer>', type=int, nargs='*',default=[3,23],
                   help='a list of isosurface color IDs (a list of integers, default = [3,23])')
    parser.add_argument('--isocut', metavar='<isovalue cutoff>', type=float, nargs='?',default=1e-8,
                   help='cutoff value for rendering an isosurface (float, default = 1.0e-8)')
                   
    parser.add_argument('--rx', metavar='<angle>', type=float, nargs='?', default=euler_angles[0],
                   help='the x-axis rotation angle (float, default = 30.0)')
    parser.add_argument('--ry', metavar='<angle>', type=float, nargs='?', default=euler_angles[1],
                   help='the y-axis rotation angle (float, default = 40.0)')
    parser.add_argument('--rz', metavar='<angle>', type=float, nargs='?', default=euler_angles[2],
                   help='the z-axis rotation angle (float, default = 15.0)')


    parser.add_argument('--tx', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the x-axis translation (float, default = 0.0)')
    parser.add_argument('--ty', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the y-axis translation (float, default = 0.0)')
    parser.add_argument('--tz', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the z-axis translation (float, default = 0.0)')

    parser.add_argument('--opacity', metavar='<opacity>', type=float, nargs='?',default=1.0,
                   help='opacity of the isosurface (float, default = 1.0)')

    parser.add_argument('--scale', metavar='<factor>', type=float, nargs='?',default=1.0,
                   help='the scaling factor (float, default = 1.0)')
    parser.add_argument('--no-montage', action="store_true",
                   help='call montage to combine images. (string, default = false)')
    parser.add_argument('--no-labels', action="store_true",
                   help='do not add labels to images. (string, default = false)')

    parser.add_argument('--imagew', metavar='<integer>', type=int, nargs='?', default=RESOLUTION,
                   help='the width of images (integer, default = 250)')
    parser.add_argument('--imageh', metavar='<integer>', type=int, nargs='?', default=RESOLUTION,
                   help='the height of images (integer, default = 250)')
    parser.add_argument('--fontsize', metavar='<integer>', type=int, nargs='?', default=20,
                   help='the font size (integer, default = 20)')

    parser.add_argument('--interactive', action="store_true",
                   help='run in interactive mode (default = false)')

    parser.add_argument('--gzip', action="store_true",
                   help='gzip cube files (default = false)')

    parser.add_argument('--national_scheme', action="store_true",
                   help='use a red/blue color scheme. (string, default = false)')
    parser.add_argument('--silver_scheme', action="store_true",
                   help='use a gray/white color scheme. (string, default = false)')
    parser.add_argument('--bright_scheme', action="store_true",
                   help='use a soft yellow/blue color scheme. (string, default = false)')
    parser.add_argument('--electron_scheme', action="store_true",
                   help='use a purple/green color scheme. (string, default = false)')

    args = parser.parse_args()

    options["CUBEDIR"][0] = str(args.data)
    options["ISOVALUE"][0] = args.isovalue
    options["ISOCOLOR"][0] = args.isocolor
    options["ISOCUT"][0] = str(args.isocut)
    options["RX"][0] = str(args.rx)
    options["RY"][0] = str(args.ry)
    options["RZ"][0] = str(args.rz)
    options["TX"][0] = str(args.tx)
    options["TY"][0] = str(args.ty)
    options["TZ"][0] = str(args.tz)
    options["OPACITY"][0] = str(args.opacity)
    options["SCALE"][0] = str(args.scale)
    options["LABEL_MOS"][0] = str(not args.no_labels)
    options["MONTAGE"][0] = str(not args.no_montage)
    options["FONTSIZE"][0] = str(args.fontsize)
    options["IMAGEW"][0] = str(args.imagew)
    options["IMAGEH"][0] = str(args.imageh)
    options["INTERACTIVE"][0] = str(args.interactive)
    options["GZIP"][0] = str(args.gzip)

    if args.national_scheme:
        options["ISOCOLOR"][0] = [23,30]

    if args.silver_scheme:
        options["ISOCOLOR"][0] = [24,32]  # silver = [2,8]

    if args.electron_scheme:
        options["ISOCOLOR"][0] = [13,12]

    if args.bright_scheme:
        options["ISOCOLOR"][0] = [32,22]

    print("\nParameters:")
    sorted_parameters = sorted(options.keys())
    for k in sorted_parameters:
        print("  %-20s %s" % (options[k][1],str(options[k][0])))

def find_cubes(options):
    # Find all the cube files in a given directory
    dir = options["CUBEDIR"][0]
    sorted_files = []
    zipped_files = []

    for f in listdir(options["CUBEDIR"][0]):
        if "\'" in f:
            nf = f.replace("\'", "p")
            os.rename(f,nf)
            f = nf
        if "\"" in f:
            nf = f.replace("\"", "pp")
            os.rename(f,nf)
            f = nf
        if f[-5:] == '.cube':
            sorted_files.append(f)
        elif f[-8:] == '.cube.gz':
            found_zipped = True
            # unzip file
            sorted_files.append(f[:-3])
            zipped_files.append(f)

    if len(zipped_files) > 0:
        print("\nDecompressing gzipped cube files")
        FNULL = open(os.devnull, 'w')
        subprocess.call(("gzip -d %s" % " ".join(zipped_files)),stdout=FNULL, shell=True)
        options["GZIP"][0] = 'True'

    return sorted(sorted_files)


def write_and_run_vmd_script(options,cube_files):
    vmd_script = open(vmd_script_name,"w+")

    # Define a map that contains all the values of the VMD parameters
    replacement_map = {}
    for (k, v) in iteritems(options):
        key = "PARAM_" + k.upper()
        replacement_map[key] = v[0]

    for n,f in enumerate(cube_files):
        replacement_map["PARAM_CUBENUM"] = '%03d' % n
        replacement_map["PARAM_CUBEFILE"] = options["CUBEDIR"][0] + '/' + f[:-5]

        # Default isocontour values or user-provided
        isovalue = options["ISOVALUE"][0][:]
        isocolor = options["ISOCOLOR"][0][:]
        
        # Read isocontour values from file, if available
        with open(f,'r') as file:
            l1 = file.readline()
            l2 = file.readline()
            m = re.search(r'density: \(([-+]?[0-9]*\.?[0-9]+)\,([-+]?[0-9]*\.?[0-9]+)\)',l2)
            if m:
                isovalue[0] = float(m.groups()[0])
                isovalue[1] = float(m.groups()[1])

        nisovalue = len(isovalue)
        nisocolor = len(isocolor)
        if nisovalue!= nisocolor:
            print("Quitting: Please specify the same number of isosurface values and colors.")
            quit()
        else:
            print("Plotting %s with isosurface values" % (f), str(isovalue))

        vmd_script_surface = ""
        surf = zip(isovalue,isocolor)
        for c in surf:
            if abs(c[0]) > float(options["ISOCUT"][0]):
                replacement_map["PARAM_ISOVALUE"] = str(c[0])
                replacement_map["PARAM_ISOCOLOR"] = str(c[1])
                vmd_script_surface += multigsub(replacement_map,vmd_template_surface)
            else:
                print(" * Skipping isosurface with isocontour value %f" % c[0])
        vmd_script_head = multigsub(replacement_map,vmd_template)
        
        if options["INTERACTIVE"][0] == 'True':
            vmd_script_render = multigsub(replacement_map, vmd_template_interactive)
        else:
            vmd_script_render = multigsub(replacement_map, vmd_template_render)

        vmd_script.write(vmd_script_head + "\n" + vmd_script_surface + "\n" + vmd_script_render)

    if options["INTERACTIVE"][0] == 'False':
        vmd_script.write("quit")
        vmd_script.close()
        # Call VMD in text mode
        FNULL = open(os.devnull, 'w')
        subprocess.call(("%s -dispdev text -e %s" % (options["VMDPATH"][0],vmd_script_name)),stdout=FNULL, shell=True)
    else:
        vmd_script.close()
        # Call VMD in graphic mode
        FNULL = open(os.devnull, 'w')
        subprocess.call(("%s -e %s" % (options["VMDPATH"][0],vmd_script_name)),stdout=FNULL, shell=True)


def call_montage(options,cube_files):
    if options["MONTAGE"][0] == 'True':
        # Optionally, combine all figures into one image using montage
        montage_exe = which("montage")
        if montage_exe:
            alpha_mos = []
            beta_mos = []
            densities = []
            basis_functions = []
            for f in cube_files:
                tga_file = f[:-5] + ".tga"
                if "Psi_a" in f:
                    alpha_mos.append(tga_file)
                if "Psi_b" in f:
                    beta_mos.append(tga_file)
                if "D" in f:
                    densities.append(tga_file)
                if "Phi" in f:
                    basis_functions.append(tga_file)

            # Sort the MOs
            sorted_mos = []
            for set in [alpha_mos,beta_mos]:
                sorted_set = []
                for s in set:
                    s_split = s.split('_')
                    sorted_set.append((int(s_split[2]),"Psi_a_%s_%s" % (s_split[2],s_split[3])))
                sorted_set = sorted(sorted_set)
                sorted_mos.append([s[1] for s in sorted_set])
           
            os.chdir(options["CUBEDIR"][0])
                    
            # Add labels
            if options["LABEL_MOS"][0] == 'True':
                for f in sorted_mos[0]:
                    f_split = f.split('_')
                    label = '%s\ \(%s\)' % (f_split[3][:-4],f_split[2])
                    subprocess.call(("montage -pointsize %s -label %s %s -geometry '%sx%s+0+0>' %s" %
                        (options["FONTSIZE"][0],label,f,options["IMAGEW"][0],options["IMAGEH"][0],f)), shell=True)

            # Combine together in one image
            if len(alpha_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 AlphaMOs.tga" % (montage_exe," ".join(sorted_mos[0]))), shell=True)
            if len(beta_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 BetaMOs.tga" % (montage_exe," ".join(sorted_mos[1]))), shell=True)
            if len(densities) > 0:
                subprocess.call(("%s %s -geometry +2+2 Densities.tga" % (montage_exe," ".join(densities))), shell=True)
            if len(basis_functions) > 0:
                subprocess.call(("%s %s -geometry +2+2 BasisFunctions.tga" % (montage_exe," ".join(basis_functions))), shell=True)


def zip_files(cube_files,options):
    """Gzip cube files if requested or necessary."""
    if options["GZIP"][0] == 'True':
        print("\nCompressing cube files")
        FNULL = open(os.devnull, 'w')
        subprocess.call(("gzip %s" % " ".join(cube_files)),stdout=FNULL, shell=True)


def get_cumulative_density_iso_value(file,sigma):
    """Find the isosurface values that capture a certain amount of the total density (sigma)."""
    cube_data = []
    norm = 0.0
    k = 0
    with open(file) as f:
        for line in f:
            if k > 6:
                for s in line.split():
                    value = float(s)
                    value_sqr = value * value
                    norm = norm + value_sqr
                    cube_data.append((value_sqr,value))
            k = k + 1

    cube_data.sort(reverse=True)

    sum = 0.0
    positive_iso = 0.0
    negative_iso = 0.0
    for (value_sqr,value) in cube_data:
        if sum < sigma:
            sum = sum + value_sqr / norm
            if value > 0:
                positive_iso = value
            else:
                negative_iso = value
        else:
            return (positive_iso, negative_iso)
    return (positive_iso, negative_iso)


def main(argv):
    find_vmd(options)
    read_options(options)
    save_setup_command(argv)
    cube_files = find_cubes(options)
    write_and_run_vmd_script(options,cube_files)
    call_montage(options,cube_files)
    zip_files(cube_files,options)

if __name__ == '__main__':
    main(sys.argv)

os.remove('.vmd_cube_command')
os.remove('.vmd_mo_script.vmd')

##############################################################################
""" This part is for stitching together the generated images. """

def stitch_images_together(filenames, spin):
    plt.figure(figsize=(2*len(filenames), len(filenames) + 1))
    for num, file in enumerate(filenames):
        img = PIL.Image.open(file)
        if file.startswith('filtered'):
            plt.subplot(2, len(filenames)//2, len(filenames)//2 + num//2 + 1)
            subfigure_name = ""
        else:
            plt.subplot(2, len(filenames)//2, num//2 + 1)
            subfigure_name = file.split('_')[3]
            subfigure_name = str(int(subfigure_name)) # Removes the zeros from the subfigure name.
            plt.title(f'State {subfigure_name}', fontsize=50)
        if num == 0:
            plt.ylabel('raw DFT', fontsize=40)
        elif num == 1:
            plt.ylabel('filtered DFT', fontsize=40)
        plt.xlim(0, RESOLUTION)
        plt.ylim(0, RESOLUTION)
        plt.xticks([])
        plt.yticks([])
        plt.box(False)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
#        plt.axis('square')
        plt.subplots_adjust(left=0.1)
        plt.imshow(img, interpolation='none', aspect='auto')
        plt.tight_layout()
        plt.savefig(f'orbitals_{spin}.svg')

filenames = {}
for spin in ['spin_1', 'spin_2']:
    filenames[spin] = []
    for file in raw_CUBE_files:
        if (file.find(spin) != -1):
            filenames[spin].append(file.replace('.cube', '.tga'))
            filenames[spin].append(f"filtered_{file.replace('.cube', '.tga')}")
    if bool(spin):
        stitch_images_together(filenames[spin], spin)

orbital.print_orbital_data()
