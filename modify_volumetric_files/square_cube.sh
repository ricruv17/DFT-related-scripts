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
	    for (i = 1; i <= NF; i++) {
	         printf("%12.4e ", \$i*\$i);
            }
	    printf("\n");
    } else {
	    print \$0;
    }   
    line++;
}
EOF

for FILE in $*
do

NFILE=`echo $FILE | sed 's/.cube/.sqr.cube/'`

awk -f tmp.awk ${FILE} > ${NFILE}

done

rm tmp.awk
