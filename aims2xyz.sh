: 'Script that converts input and output FHI-aims files to .xyz format.

Parameters (in positional order)
________________
1: string
    Name/adress of the FHI-aims file to process. 
    It can be either aims.out for initial and final positions,
    or geometry.in for only initial positions,
    or geometry.in.next_step for only final positions.

Warnings
________________
Only works if the filename has geometry. as prefix or .out as suffix.

To execute run:
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


# 1. Create xyz file of the initial positions
if [[ "$1" == *.out ]]
then
    N=$(grep "Number of atoms        " $1 | awk '{print $6}')
    echo $N > geometry.in.xyz
    echo "Generated from FHI-aims input" >> geometry.in.xyz
    grep ": Species " $1 | awk -f .geometry.in.awk >> geometry.in.xyz
elif [[ "$1" == geometry* ]]
then
    N=$(grep atom $1 | wc -l)
    echo $N > geometry.in.xyz
    echo "Generated from FHI-aims input" >> geometry.in.xyz
    grep atom $1 | awk -f .geometry.out.awk >> geometry.in.xyz
fi
echo "File geometry.in.xyz generated."

# 2. Create xyz file of the final positions or send message in case the calculation did not finish properly
file_ending_line_1=$(tail -n 2 $1 | head -n 1 | tr -s " ") 
file_ending_line_2=$(tail -n 3 $1 | head -n 1 | tr -s " ") # in case you printed date and time at the end of the aims.out file
test_script=" Have a nice day."
if [[ "$1" == *out ]]
then
    if [[ "$file_ending_line_1" == "$test_script" ]] || [[ "$file_ending_line_2" == "$test_script" ]]
    then
        echo $N > geometry.out.xyz
        echo "Generated from FHI-aims output" >> geometry.out.xyz
        grep "    atom       " $1 | tail -n ${N} | awk -f .geometry.out.awk >> geometry.out.xyz
        echo "File geometry.out.xyz generated."
    else
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
