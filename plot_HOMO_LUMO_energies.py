"""
Script to plot the HOMO and LUMO (plus any amount of neighboring orbitals) orbitals from an output file of FHI-aims.
Written by Ricardo Ruvalcaba at MONA-group in King Abdullah University of Science and Technology (KAUST).
Contact: ricardo.ruvalcababriones@kaust.edu.sa
Version 2 (19/03/23). No known bugs.
This is version 3 (04/05/23).
    - Added capability to read file name from command line input.

When executed, the script will ask for all the necessary information.
To execute, just run on your bash terminal:
python3 plot_HOMO_LUMO_energies.py FILENAME.out

Or the default:
python3 plot_HOMO_LUMO_energies.py  # ->  will seek for the file "aims.out",
                                         otherwise will ask for a correct file name.
"""

# 0. Import useful packages
import matplotlib.pyplot as plt
import sys

# 1. Define useful functions and variables
font_size = 16


def find_file():
    # Searches for a file with a given name in the current folder.
    arguments = sys.argv
    if len(arguments) > 1:
        try:
            filename = sys.argv[1]
            open(filename)
            print(f'\nFile {filename} found in the current folder.\nOpening now...')
            return filename
        except FileNotFoundError:
            print(f'\nERROR: File {filename} was not found in the current folder.')
            sys.argv[1] = input('Please enter the name of a valid FHI-aims output file: ')
            return find_file()
    else:
        try:
            open('aims.out')
            print(f'\nFile aims.out found in the current folder.\nOpening now...')
            return 'aims.out'
        except FileNotFoundError:
            print(f'\nFile aims.out was not found in the current folder.')
            sys.argv.append(input('Please enter the name of your FHI-aims output file: '))
            return find_file()


def store_files_in_list():
    filename = find_file()
    with open(filename) as file_in:
        lines = []
        for line in file_in:
            lines.append(line)
    if 'Invoking FHI-aims ...' not in lines[1]:
        print(f'\nERROR: File {filename} is not a FHI-aims output file.')
        sys.argv[1] = input('Please enter the name of a valid FHI-aims output file:\n')
        return store_files_in_list()
    return lines


def save_no_of_orbitals(message):
    try:
        no_orbitals = float(input(message))
        if no_orbitals < 0 or no_orbitals % 1 != 0:
            print(f'\nERROR: Number must be a non-negative integer.')
            return save_no_of_orbitals('Please enter a valid number: ')
        return int(no_orbitals)
    except ValueError:
        print(f'\nERROR: Number must be a non-negative integer.')
        return save_no_of_orbitals('Please enter a valid number: ')


def ask_type_of_energy():
    selection = input('\nSelect how you want to plot the energy values: (1) rounded or (2) precise: ')
    while selection != '1' and selection != '2':
        print('\nERROR: Only values 1 or 2 are allowed.')
        selection = input('Please enter a valid number: ')
    if selection == '1':
        is_rounded = True
    elif selection == '2':
        is_rounded = False
    return is_rounded


# 2. Determine the number of KS states and the energy tolerance in the system
file = store_files_in_list()

if 'Have a nice day.' not in file[-2]:
    print('\nWARNING: the calculation did not end properly.\nResults may be inaccurate.')

no_of_KS_states_found = False
is_energy_tolerance_found = False
for index, line in enumerate(file):
    if 'Number of Kohn-Sham states (occupied + empty):' in line:
        no_of_KS_states = file[index].split()
        no_of_KS_states_found = True
    if 'energy_tolerance' in line:
        energy_tolerance = file[index].split()
        is_energy_tolerance_found = True
    if no_of_KS_states_found and is_energy_tolerance_found:
        break

no_of_KS_states = int(no_of_KS_states[-1])
energy_tolerance = float(energy_tolerance[-1])  # in eV
energy_tolerance_decimal_string = ("%.10f" % energy_tolerance).rstrip('0').rstrip('.')
# rounding_criteria = energy_tolerance_decimal_string[::-1].find('.')

# 3. Determine the position fo the final KS states
file.reverse()
eigenvalues_index = len(file) - file.index('  Writing Kohn-Sham eigenvalues.\n') - 1
file.reverse()

#################### Spin-polarized calculation #####################
calculation_is_spin_polarized = 'Spin-up eigenvalues' in file[eigenvalues_index + 2]
if calculation_is_spin_polarized:
    # 4. Create a list containing the orbitals, their occupations, and energies; and find HOMO and LUMO indexes
    print('\nCalculation is spin polarized.\n')
    up_orbitals_data, dn_orbitals_data = [[0, 0, 0, 0]], [[0, 0, 0, 0]]
    up_LUMO_index, dn_LUMO_index = None, None
    up_HOMO_index, dn_HOMO_index = None, None
    for line in file[eigenvalues_index + 5: eigenvalues_index + 5 + no_of_KS_states]:
        up_orbital = list(map(float, line.split()))
        up_orbital[0] = int(up_orbital[0])
        up_orbitals_data.append(up_orbital)
        up_occupation = up_orbital[1]
        if (up_LUMO_index is None) and (up_occupation < 0.3):
            up_LUMO_index = up_orbital[0]
            up_HOMO_index = up_LUMO_index - 1
    for line in file[eigenvalues_index + 5 + no_of_KS_states + 4: eigenvalues_index + 5 + 2*no_of_KS_states + 4]:
        dn_orbital = list(map(float, line.split()))
        dn_orbital[0] = int(dn_orbital[0])
        dn_orbitals_data.append(dn_orbital)
        dn_occupation = dn_orbital[1]
        if (dn_LUMO_index is None) and (dn_occupation < 0.3):
            dn_LUMO_index = dn_orbital[0]
            dn_HOMO_index = dn_LUMO_index - 1

    # 5. Create a list containing the orbitals to be plotted, print to terminal, and create a dictionary containing the
    # degenerate orbitals according to energy_tolerance criterion
    up_no_OMOs = save_no_of_orbitals('Enter how many spin-up orbitals below HOMO you want to plot: ')
    if up_HOMO_index - up_no_OMOs <= 0:
        up_max_available_OMO_orbitals = up_HOMO_index - 1
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({up_max_available_OMO_orbitals}) will be used instead.\n')
        up_no_OMOs = up_max_available_OMO_orbitals

    up_no_UMOs = save_no_of_orbitals('Enter how many spin-up orbitals above LUMO you want to plot: ')
    if up_LUMO_index + up_no_UMOs > no_of_KS_states:
        up_max_available_UMO_orbitals = no_of_KS_states - up_LUMO_index
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({up_max_available_UMO_orbitals}) will be used instead.')
        up_no_UMOs = up_max_available_UMO_orbitals

    dn_no_OMOs = save_no_of_orbitals('Enter how many spin-down orbitals below HOMO you want to plot: ')
    if dn_HOMO_index - dn_no_OMOs <= 0:
        dn_max_available_OMO_orbitals = dn_HOMO_index - 1
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({dn_max_available_OMO_orbitals}) will be used instead.\n')
        dn_no_OMOs = dn_max_available_OMO_orbitals

    dn_no_UMOs = save_no_of_orbitals('Enter how many spin-down orbitals above LUMO you want to plot: ')
    if dn_LUMO_index + dn_no_UMOs > no_of_KS_states:
        dn_max_available_UMO_orbitals = no_of_KS_states - dn_LUMO_index
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({dn_max_available_UMO_orbitals}) will be used instead.')
        dn_no_UMOs = dn_max_available_UMO_orbitals
    
    print('\n  Selected eigenvalues:')
    print('  State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]')
    print('  Spin-up eigenvalues')
    line_lower_limit = eigenvalues_index + 5 + up_HOMO_index - up_no_OMOs - 1
    line_upper_limit = eigenvalues_index + 5 + up_LUMO_index + up_no_UMOs
    for line in file[line_lower_limit: line_upper_limit]:
        print(line[:-1])
    print('  Spin-down eigenvalues')
    line_lower_limit = eigenvalues_index + 5 + no_of_KS_states + 4 + dn_HOMO_index - dn_no_OMOs - 1
    line_upper_limit = eigenvalues_index + 5 + no_of_KS_states + 4 + dn_LUMO_index + dn_no_UMOs
    for line in file[line_lower_limit: line_upper_limit]:
        print(line[:-1])

    print(f'\nEnergy tolerance in your calculation is ({energy_tolerance_decimal_string} eV).')
    rounding_criteria = save_no_of_orbitals('Enter how many digits after the decimal point you want to'
                                            ' consider for the energy (in eV): ')
    up_plotted_orbitals = up_orbitals_data[up_HOMO_index - up_no_OMOs: up_LUMO_index + up_no_UMOs + 1]
    up_orbitals_degeneracy = dict()
    orbital_index = -1
    up_max_degeneracy = 0
    for orbital in up_plotted_orbitals:
        orbital_index += 1        
        energy_eV = orbital[3]
        rounded_energy = round(energy_eV, rounding_criteria)
        if rounded_energy not in up_orbitals_degeneracy:
            up_orbitals_degeneracy[rounded_energy] = [orbital_index]
        else:
            up_orbitals_degeneracy[rounded_energy].append(orbital_index)
        degeneracy = len(up_orbitals_degeneracy[rounded_energy])
        if degeneracy > up_max_degeneracy:
            up_max_degeneracy = degeneracy

    dn_plotted_orbitals = dn_orbitals_data[dn_HOMO_index - dn_no_OMOs: dn_LUMO_index + dn_no_UMOs + 1]
    dn_orbitals_degeneracy = dict()
    orbital_index = -1
    dn_max_degeneracy = 0
    for orbital in dn_plotted_orbitals:
        orbital_index += 1        
        energy_eV = orbital[3]
        rounded_energy = round(energy_eV, rounding_criteria)
        if rounded_energy not in dn_orbitals_degeneracy:
            dn_orbitals_degeneracy[rounded_energy] = [orbital_index]
        else:
            dn_orbitals_degeneracy[rounded_energy].append(orbital_index)
        degeneracy = len(dn_orbitals_degeneracy[rounded_energy])
        if degeneracy > dn_max_degeneracy:
            dn_max_degeneracy = degeneracy

    # 6. Plot orbital energies
    bar_width = 0.4
    HOMO_color = 'tab:blue'
    LUMO_color = 'tab:orange'
    line_width = 1
    with_rounded_energies = ask_type_of_energy()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    max_E_Ha = max(up_orbitals_data[up_LUMO_index + up_no_UMOs][2], dn_orbitals_data[dn_LUMO_index + dn_no_UMOs][2])
    min_E_Ha = min(up_orbitals_data[up_HOMO_index - up_no_OMOs][2], dn_orbitals_data[dn_HOMO_index - dn_no_OMOs][2])
    ax.plot([-up_max_degeneracy*101, -up_max_degeneracy*100], [max_E_Ha, max_E_Ha], color=LUMO_color, linewidth=2,
            label='Unoccupied orbitals')
    ax.plot([-up_max_degeneracy*101, -up_max_degeneracy*100], [min_E_Ha, min_E_Ha], color=HOMO_color, linewidth=2,
            label='Occupied orbitals')
    ax.legend(bbox_to_anchor=(1.25, -0.05))
    for rounded_energy in up_orbitals_degeneracy:
        indexes_in_same_energy = up_orbitals_degeneracy[rounded_energy]
        degeneracy = len(indexes_in_same_energy)
        step = -1
        for index in indexes_in_same_energy:
            step += 1
            if with_rounded_energies:
                plot_energy = round(up_plotted_orbitals[index][3], rounding_criteria)
            else:
                plot_energy = up_plotted_orbitals[index][3]
            if up_plotted_orbitals[index][0] >= up_LUMO_index:
                orbital_color = LUMO_color
            else:
                orbital_color = HOMO_color        
            x_position = -up_max_degeneracy/2 - degeneracy/2 + step
            ax.plot([x_position - bar_width, x_position + bar_width], [plot_energy, plot_energy], color=orbital_color,
                    linewidth=line_width)
    for rounded_energy in dn_orbitals_degeneracy:
        indexes_in_same_energy = dn_orbitals_degeneracy[rounded_energy]
        degeneracy = len(indexes_in_same_energy)
        step = -1
        for index in indexes_in_same_energy:
            step += 1
            if with_rounded_energies:
                plot_energy = round(dn_plotted_orbitals[index][3], rounding_criteria)
            else:
                plot_energy = dn_plotted_orbitals[index][3]
            if dn_plotted_orbitals[index][0] >= dn_LUMO_index:
                orbital_color = LUMO_color
            else:
                orbital_color = HOMO_color        
            x_position = dn_max_degeneracy/2 - degeneracy/2 + step + 1
            ax.plot([x_position - bar_width, x_position + bar_width], [plot_energy, plot_energy], color=orbital_color,
                    linewidth=line_width)
    ax2 = ax.twinx()
    ax2.plot([-up_max_degeneracy*101, -up_max_degeneracy*100], [max_E_Ha, max_E_Ha], color='k', linewidth=2)
    ax2.plot([-up_max_degeneracy*101, -up_max_degeneracy*100], [min_E_Ha, min_E_Ha], color='k', linewidth=2)

    ax.set(xlim=[-(up_max_degeneracy + 0.5) * 1.1, (dn_max_degeneracy + 0.5) * 1.1])
    ax.set_xlabel("Spin", fontsize=font_size, labelpad=0)
    ax.set_ylabel(r"Energy (eV)", fontsize=font_size, labelpad=2)
    ax2.set_ylabel(r"Energy (Ha)", fontsize=font_size, labelpad=5)
    ax.set_xticks([-up_max_degeneracy/2 - 0.5, dn_max_degeneracy/2 + 0.5])
    ax.set_xticklabels([r"$\uparrow$", r"$\downarrow$"], fontsize=font_size)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax2.tick_params(axis='both', which='major', labelsize=font_size)

    plt.tight_layout()
    plt.show()

################### Non-spin-polarized calculation ##################
else:
    # 4. Create a list containing the orbitals, their occupations, and energies; and find HOMO and LUMO indexes
    print('\nCalculation is non-spin polarized.\n')
    orbitals_data = [[0, 0, 0, 0]]
    LUMO_index = None
    HOMO_index = None
    for line in file[eigenvalues_index + 3: eigenvalues_index + 3 + no_of_KS_states]:
        orbital = list(map(float, line.split()))
        orbital[0] = int(orbital[0])
        orbitals_data.append(orbital)
        occupation = orbital[1]
        if (LUMO_index is None) and (occupation < 0.6):
            LUMO_index = orbital[0]
            HOMO_index = LUMO_index - 1

    # 5. Create a list containing the orbitals to be plotted, print to terminal, and create a dictionary containing the
    # degenerate orbitals according to energy_tolerance criterion
    no_OMOs = save_no_of_orbitals('Enter how many orbitals below HOMO you want to plot: ')
    if HOMO_index - no_OMOs <= 0:
        max_available_OMO_orbitals = HOMO_index - 1
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({max_available_OMO_orbitals}) will be used instead.\n')
        no_OMOs = max_available_OMO_orbitals

    no_UMOs = save_no_of_orbitals('Enter how many orbitals above LUMO you want to plot: ')
    if LUMO_index + no_UMOs > no_of_KS_states:
        max_available_UMO_orbitals = no_of_KS_states - LUMO_index
        print(f'This number exceeds the number of available orbitals.'
              f'\nThe maximum available ({max_available_UMO_orbitals}) will be used instead.')
        no_UMOs = max_available_UMO_orbitals
    
    print('\n  Selected eigenvalues:')
    print('  State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]')

    for line in file[eigenvalues_index + 3 + HOMO_index - no_OMOs - 1: eigenvalues_index + 3 + LUMO_index + no_UMOs]:
        print(line[:-1])

    print(f'\nEnergy tolerance in your calculation is ({energy_tolerance_decimal_string} eV).')
    rounding_criteria = save_no_of_orbitals('Enter how many digits after the decimal point you want to '
                                            'consider for the energy (in eV): ')
    plotted_orbitals = orbitals_data[HOMO_index - no_OMOs: LUMO_index + no_UMOs + 1]
    orbitals_degeneracy = dict()
    orbital_index = -1
    max_degeneracy = 0
    for orbital in plotted_orbitals:
        orbital_index += 1        
        energy_eV = orbital[3]
        rounded_energy = round(energy_eV, rounding_criteria)
        if rounded_energy not in orbitals_degeneracy:
            orbitals_degeneracy[rounded_energy] = [orbital_index]
        else:
            orbitals_degeneracy[rounded_energy].append(orbital_index)
        degeneracy = len(orbitals_degeneracy[rounded_energy])
        if degeneracy > max_degeneracy:
            max_degeneracy = degeneracy

    # 6. Plot orbital energies
    bar_width = 0.4
    HOMO_color = 'tab:blue'
    LUMO_color = 'tab:orange'
    line_width = 1
    with_rounded_energies = ask_type_of_energy()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    max_E_Ha = plotted_orbitals[-1][2]
    min_E_Ha = plotted_orbitals[-1][2]
    ax.plot([-max_degeneracy*1.1, -max_degeneracy], [max_E_Ha, max_E_Ha], color=LUMO_color, linewidth=2,
            label='Unoccupied orbitals')
    ax.plot([-max_degeneracy*1.1, -max_degeneracy], [min_E_Ha, min_E_Ha], color=HOMO_color, linewidth=2,
            label='Occupied orbitals')
    ax.legend(bbox_to_anchor=(1.25, 0))
    for rounded_energy in orbitals_degeneracy:
        indexes_in_same_energy = orbitals_degeneracy[rounded_energy]
        degeneracy = len(indexes_in_same_energy)
        step = -1
        for index in indexes_in_same_energy:
            step += 1
            if with_rounded_energies:
                plot_energy = round(plotted_orbitals[index][3], rounding_criteria)
            else:
                plot_energy = plotted_orbitals[index][3]
            if plotted_orbitals[index][0] >= LUMO_index:
                orbital_color = LUMO_color
            else:
                orbital_color = HOMO_color        
            x_position = -degeneracy/2 + step + 0.5
            ax.plot([x_position - bar_width, x_position + bar_width], [plot_energy, plot_energy], color=orbital_color,
                    linewidth=line_width)
    
    ax2 = ax.twinx()
    ax2.plot([-max_degeneracy*1.1, -max_degeneracy], [max_E_Ha, max_E_Ha], color='k', linewidth=2)
    ax2.plot([-max_degeneracy*1.1, -max_degeneracy], [min_E_Ha, min_E_Ha], color='k', linewidth=2)

    ax.set(xlim=[-max_degeneracy/2 * 1.1, max_degeneracy/2 * 1.1])
    ax.set_xlabel("Spin", fontsize=font_size, labelpad=0)
    ax.set_ylabel(r"Energy (eV)", fontsize=font_size, labelpad=2)
    ax2.set_ylabel(r"Energy (Ha)", fontsize=font_size, labelpad=5)
    ax.set_xticks([0])
    ax.set_xticklabels([r"$\uparrow \downarrow$"], fontsize=font_size)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax2.tick_params(axis='both', which='major', labelsize=font_size)

    plt.tight_layout()
    plt.show()
