#!/bin/bash

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
    } else {
		for (i = 1; i <= NF/2; i++) {
	         printf("%s ", \$i);
            }
			printf("\n");
    }   
    line++;
}
EOF

natoms1=$(head -n 3 $1 | tail -1 | awk '{print $1}')
natoms2=$(head -n 3 $2 | tail -1 | awk '{print $1}')

head -n 6 $1 > zeros_prefix
cp zeros_prefix prefix
sed -i -e "s/  $natoms1  /  0  /g" zeros_prefix
sed -i -e "s/  $natoms1  /  $(echo "$natoms1 + $natoms2" | bc)  /g" prefix

head -n $(echo "$natoms1 + 6" | bc) $1 | tail -n $natoms1 > file1_atoms
head -n $(echo "$natoms2 + 6" | bc) $2 | tail -n $natoms2 > file2_atoms

tail -n +$(echo "$natoms1 + 6 + 1" | bc) $1 > volume_data1
tail -n +$(echo "$natoms2 + 6 + 1" | bc) $2 > volume_data2

cat zeros_prefix volume_data1 > process1
cat zeros_prefix volume_data2 > process2

paste process1 process2 | awk -f tmp.awk > tmp_sum
tail -n +7 tmp_sum > volume_sum
cat prefix file{1,2}_atoms volume_sum > sum.cube

rm tmp.awk zeros_prefix prefix file{1,2}_atoms volume_data{1,2} process{1,2} {tmp,volume}_sum
