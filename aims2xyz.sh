: 'Script that converts input and output FHI-aims files to .xyz format.

 Parameters (in positional order)
-----------------
 1: string
    Name/adress of the FHI-aims file to process. 
    It can be either aims.out for initial and final positions,
    or geometry.in for only initial positions,
    or geometry.in.next_step for only final positions.

 Returns
-----------------
 Files in .xyz format with the atomic coordinates of the system.
 Their name depends on the input file:
    geometry.in.xyz: initial coordinates. Generated from aims.out or geometry.in files.
    geometry.out.xyz: final coordinates. Generated from aims.out file.
    geometry.in.next_step.xyz: intermediate/final coordinales. Generated from geometry.in.next_step file.

 Warnings
-----------------
 Only works if the filename has geometry. as prefix or .out as suffix.

To execute, just run:
bash aims2xyz.sh FILENAME
'
# 0. Create temporal files containing format instructions for xyz files
cat > .geometry.in.awk << EOF
{
        x=\$5;
        y=\$6;
        z=\$7;
        lab=\$4;
        printf("%10s %20.8f %20.8f %20.8f\n", lab, x, y, z);
}
EOF

cat > .geometry.out.awk << EOF
{
        x=\$2;
        y=\$3;
        z=\$4;
        lab=\$5;
        printf("%10s %20.8f %20.8f %20.8f\n", lab, x, y, z);
}
EOF

# 1. Process files with geometry. as a prefix.
if [[ "$1" == geometry* ]]
then
    N=$(grep atom $1 | wc -l)
    echo $N > ${1}.xyz
    echo "Generated from FHI-aims input" >> ${1}.xyz
    grep atom $1 | awk -f .geometry.out.awk >> ${1}.xyz
    echo "File ${1}.xyz generated."
# 2. Process files with .out as a suffix.
elif [[ "$1" == *out ]]
then
    # Create .xyz file of the input coordinates.
    N=$(grep "Number of atoms        " $1 | awk '{print $6}')
    echo $N > geometry.in.xyz
    echo "Generated from FHI-aims input" >> geometry.in.xyz
    grep ": Species " $1 | awk -f .geometry.in.awk >> geometry.in.xyz
    echo "File geometry.in.xyz generated."

    # Set criteria to evaluate if calculation contains relaxed output coordinates.
    file_ending_line_1=$(tail -n 2 $1 | head -n 1 | tr -s " ") 
    file_ending_line_2=$(tail -n 3 $1 | head -n 1 | tr -s " ") # in case you printed date and time at the end of the aims.out file
    test_script=" Have a nice day."

    if [[ "$file_ending_line_1" == "$test_script" ]] || [[ "$file_ending_line_2" == "$test_script" ]]
    then
        # Create .xyz file of the output coordinates.
        echo $N > geometry.out.xyz
        echo "Generated from FHI-aims output" >> geometry.out.xyz
        grep "    atom       " $1 | tail -n ${N} | awk -f .geometry.out.awk >> geometry.out.xyz
        echo "File geometry.out.xyz generated."
    else
        # Warn user if there are no output coordinates in the file.
        printf "%$(bc -l <<< "$(tput cols)/2 - 7")s %9s %$(bc -l <<< "$(tput cols)/2 - 4")s" " " "_WARNING_" " " | tr " " "=" | tr "_" " "
        echo "File does not contain relaxed positions."
        echo "Please make sure the calculation ended properly."
        echo "If it did, check calculation parameters in control.in and make sure you relaxed the structure."
        echo "File geometry.out.xyz not generated."
        printf "%$(tput cols)s" " " | tr " " "="
    fi
fi

# 3. Remove temporal files
rm .geometry.in.awk .geometry.out.awk
