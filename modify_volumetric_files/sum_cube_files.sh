: 'Script that sums the volume data of two CUBE files.

 Parameters (in positional order)
-----------------
 1: string
    Name/adress of the first CUBE file to process.
 2: string
    Name/adress of the second CUBE file to process.

 Returns
-----------------
 Prints to terminal the same file in CUBE format with its volumetric data squared.
 
 Warnings
-----------------
 The program assumes both CUBE files have the same origin and grid density as FILE_1.

To execute, just run:
bash aims2xyz.sh FILE_1.cube FILE_2.cube

And to save it to a new file, run:
bash aims2xyz.sh FILE_1.cube FILE_2.cube > NEW_FILE.cube
'

cat > tmp.awk <<EOF
BEGIN {
    line = 1;
    atoms = 9999999;
}
{
    if (line == 3) {
	atoms = \$1;
    }

    if (line > atoms + 6) {
	    for (i = 1; i <= NF/2; i++) {
			 j=i+NF/2;
	         printf("%12.4e ", \$i+\$j);
            }
	    printf("\n");
    }   
    line++;
}
EOF

file_1=$1
file_2=$2

natoms1=$(head -n 3 $file_1 | tail -1 | awk '{print $1}')
natoms2=$(head -n 3 $file_2 | tail -1 | awk '{print $1}')

head -n 2 $file_1 > temp_file_1
head -n 3 $file_1 | tail -1 | awk '{print "   ", 0, $2, $3, $4}' >> temp_file_1
head -n 6 $file_1 | tail -3 >> temp_file_1
tail -n +$(echo "7 + $natoms1" | bc) $file_1 >> temp_file_1
head -n 2 $file_2 > temp_file_2
head -n 3 $file_2 | tail -1 | awk '{print "   ", 0, $2, $3, $4}' >> temp_file_2
head -n 6 $file_2 | tail -3 >> temp_file_2
tail -n +$(echo "7 + $natoms2" | bc) $file_2 >> temp_file_2

head -n 2 temp_file_1
natoms_tot=$(echo "$natoms1 + $natoms2" | bc)
head -n 3 temp_file_1 | tail -1 | awk -v var="$natoms_tot" '{print "   ", var, $2, $3, $4}'
head -n 6 temp_file_1 | tail -3
head -n $( echo "$natoms1 + 6" | bc) $file_1 | tail -$natoms1
head -n $( echo "$natoms2 + 6" | bc) $file_2 | tail -$natoms2

paste temp_file_1 temp_file_2 | awk -f tmp.awk

rm tmp.awk temp_file_1 temp_file_2

