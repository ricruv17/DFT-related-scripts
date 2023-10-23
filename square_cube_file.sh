: 'Script that squares the volume data in a CUBE file.

 Parameters (in positional order)
-----------------
 1: string
    Name/adress of the CUBE file to process. 

 Returns
-----------------
 Prints to terminal the same file in CUBE format with its volumetric data squared.
 
To execute, just run:
bash aims2xyz.sh FILENAME.cube

And to save it to a new file, run:
bash aims2xyz.sh FILENAME.cube > NEW_FILE.cube
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

#for FILE in *{1,2}.cube # $*
#do
#NFILE=`echo $FILE | sed 's/.cube/.sqr.cube/'`
#awk -f tmp.awk ${FILE} > ${NFILE}
#done

FILE=$1
awk -f tmp.awk ${FILE}

rm tmp.awk
