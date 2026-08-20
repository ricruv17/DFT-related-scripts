: 'Script that converts .xyz files to geometry.in files in FHI-aims format.

 Parameters (in positional order)
-----------------
 1: string
    Name/adress of the .xyz file to process.
    It can have any name, as long as it has a .xyz extension.

 Returns
-----------------
 File in FHI-aims geometry.in format with the atomic coordinates of the system.
 The name is the original name with the .xyz extension changed for a geometry.in extension.

To execute, just run:
bash xyz2aims.sh FILENAME
'

filename=$1
cat $filename | tail -n +3 | awk '{print "atom      " $2 "      " $3 "      " $4 "      " $1}' > .tmp;
new_name=${filename::-4}.geometry.in
if [ $new_name == geometry.in.geometry.in ]
then
   new_name=geometry.in
fi
tr -d $'\r' < .tmp | rev | column -t | rev > $new_name
echo "File $new_name generated."
rm .tmp

