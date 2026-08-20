import sys
import numpy as np

class XYZ_file:
    # Represents a XYZ file with the coordinates of an atomic-molecular system.
    def __init__(self):
        self.elements = []
        self.x = []
        self.y = []
        self.z = []
        self.no_atoms = 0
        self.original_filename = ''

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
        self.original_filename = self.find_file()
        with open(self.original_filename) as file_in:
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

    def translate_molecules_centroid_to_origin(self):
        # Name is self-explanatory.
        centroid_x = sum(self.x)/self.no_atoms
        centroid_y = sum(self.y)/self.no_atoms
        centroid_z = sum(self.z)/self.no_atoms
        for index in range(self.no_atoms):
            self.x[index] -= centroid_x
            self.y[index] -= centroid_y
#            self.z[index] -= centroid_z

    def write_xyz_file(self):
        # Write a new XYZ file with the coordinates of the atoms translated to the centroid.
        file = open('centered_' + self.original_filename, "w") 
        file.write(str(self.no_atoms) + '\n')
        file.write('Centroid of the system is the origin.' + '\n')
        for index in range(self.no_atoms):
            if len(str(self.elements[index])) == 1:
                element = ' ' + str(self.elements[index])
            else:
                element = str(self.elements[index])
            line = f'{element}   {self.x[index]:12.8f}   {self.y[index]:12.8f}   {self.z[index]:12.8f}'
            file.write(line + '\n')
        file.write('\n')
        file.close()
        print(f'Done. File centered_{self.original_filename} created.\n')

file = XYZ_file()
file.read_xyz_file()
file.translate_molecules_centroid_to_origin()
file.write_xyz_file()
