: 'Script tht generates the xyz file of the initial and final positions from the aims calculation. It takes the output file from aims as only argument.

To execute run:
bash aims2xyz.sh OUTPUT_FILE.out
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
N=$(grep "Number of atoms        " $1 | awk '{print $6}')
echo $N > geometry.in.xyz
echo "Generated from FHI-aims input" >> geometry.in.xyz
grep ": Species " $1 | awk -f .geometry.in.awk >> geometry.in.xyz


# 2. Create xyz file of the final positions or send message in case the calculation did not finish properly
file_ending_line_1=$(tail -n 3 $1 | head -n 1 | tr -s " ")
file_ending_line_2=$(tail -n 2 $1 | head -n 1 | tr -s " ")
test_script=" Have a nice day."
if [[ "$file_ending_line_1" == "$test_script" ]] || [[ "$file_ending_line_2" == "$test_script" ]]
then
    echo $N > geometry.out.xyz
    echo "Generated from FHI-aims output" >> geometry.out.xyz
    grep "    atom       " $1 | tail -n ${N} | awk -f .geometry.out.awk >> geometry.out.xyz
else
    printf "%$(bc -l <<< "$(tput cols)/2 - 7")s %9s %$(bc -l <<< "$(tput cols)/2 - 4")s" " " "_WARNING_" " " | tr " " "=" | tr "_" " "
    printf "Calculation did not end properly.\ngeometry.out.xyz file not generated.\n"
    printf "%$(tput cols)s" " " | tr " " "="
fi


# 3. Remove temporal files
rm .geometry.in.awk .geometry.out.awk
