#!/bin/bash
#SBATCH --account=k1663
#SBATCH --partition=workq
#SBATCH --job-name="ppAFM"
#SBATCH --nodes=1
#SBATCH --sockets-per-node=32 ##standard value
#SBATCH --ntasks-per-node=16 ## do not change
#SBATCH --error=error.err
#SBATCH --output=output.out
#SBATCH --distribution=*:fcyclic:*
#SBATCH --time=2:00:00
#SBATCH --err=ppAFM.error
#SBATCH --output=x.o
#----------------------------------------------------------#
module load cray-python
export PYTHONPATH=/project/k1663/ruvalcrm/.local:$PYTHONPATH
location=/project/k1663/ruvalcrm/.local/bin
#----------------------------------------------------------#
function elapsed_time {
  myinput=$1
  start=$(head -n 1 $myinput)
  start=$(date --date "$start" +%s)
  end=$(tail -n 1 $myinput)
  end=$(date --date "$end" +%s)
  duration=$(echo "print($end - $start)" | python3)
  hours=$(echo "print($duration//3600)" | python3)
  minutes=$(echo "print($duration//60 - $hours*60)" | python3)
  seconds=$(echo "print($duration%60)" | python3)
  echo "$hours hours, $minutes minutes and $seconds seconds elapsed" >> $myinput
}

filename=geometry.out.xyz

date > ppAFM.out
python3 $location/ppafm-generate-ljff -i $filename >> ppAFM.out
python3 $location/ppafm-generate-elff-point-charges -i $filename >> ppAFM.out
python3 $location/ppafm-relaxed-scan >> ppAFM.out
python3 $location/ppafm-plot-results --df --cbar  >> ppAFM.out
date >> ppAFM.out
elapsed_time ppAFM.out
rm *.xsf */*.xsf

folder=$(pwd)
cd ..
mv $folder ~/ppAFM

