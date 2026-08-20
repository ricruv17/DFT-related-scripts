module unload xsorb

for a in $@
do
  a=$(echo $a | tr -d /)

  xmin=$(tail -n +3 $a/geometry.out.xyz | awk '{print $2}' | sed -r '/^\s*$/d' | sort -n | head -1)
  ymin=$(tail -n +3 $a/geometry.out.xyz | awk '{print $3}' | sed -r '/^\s*$/d' | sort -n | head -1)
  xmax=$(tail -n +3 $a/geometry.out.xyz | awk '{print $2}' | sed -r '/^\s*$/d' | sort -n | tail -1)
  ymax=$(tail -n +3 $a/geometry.out.xyz | awk '{print $3}' | sed -r '/^\s*$/d' | sort -n | tail -1)
  xgrid=$(echo "$xmax - $xmin" | bc -l)
  ygrid=$(echo "$ymax - $ymin" | bc -l)

  cp files_4_running/params.ini $a
  sed -i -e "s/__xgrid__/$xgrid/g" $a/params.ini
  sed -i -e "s/__ygrid__/$ygrid/g" $a/params.ini
  sed -i -e "s/__xmin__/$xmin/g" $a/params.ini
  sed -i -e "s/__xmax__/$xmax/g" $a/params.ini
  sed -i -e "s/__ymin__/$ymin/g" $a/params.ini
  sed -i -e "s/__ymax__/$ymax/g" $a/params.ini

  cp files_4_running/run_ibex.sh $a
  sed -i -e "s/__name__/$a/g" $a/run_ibex.sh

  mv $a /ibex/user/ruvalcrm
  cd /ibex/user/ruvalcrm/$a
  sbatch run_ibex.sh
  cd ~/ppAFM
done
