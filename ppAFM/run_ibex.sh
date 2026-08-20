#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="__name__"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32 ## Increase for large systems (require more memory)
#SBATCH --distribution=*:fcyclic:*
#SBATCH --time=2:00:00
#SBATCH --err=ppAFM.error
#SBATCH --output=x.o
#----------------------------------------------------------#
export PYTHONPATH=~/.local:$PYTHONPATH
location=~/.local/bin
module unload xsorb
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

function LJ_PointCharges {
python3 $location/ppafm-generate-ljff -i geometry.out.xyz >> ppAFM.out
python3 $location/ppafm-generate-elff-point-charges -i geometry.out.xyz >> ppAFM.out
python3 $location/ppafm-relaxed-scan >> ppAFM.out
python3 $location/ppafm-plot-results --df --WSxM >> ppAFM.out
#python3 $location/ppafm-plot-results --df --cbar --save_df >> ppAFM.out
#find . -name *xsf -delete
rm *.xsf */*.xsf */*/*.xsf
}

function LJ_Hartree {
python3 $location/ppafm-generate-ljff -i geometry.out.xyz >> ppAFM.out
python3 $location/ppafm-generate-elff -i cube_001_hartree_potential.cube >> ppAFM.out
python3 $location/ppafm-relaxed-scan >> ppAFM.out
python3 $location/ppafm-plot-results --df --cbar --save_df >> ppAFM.out
find . -name *xsf -delete
}

function FDBM {
python3 $location/ppafm-conv-rho -s $1 -t $2 -B 1.0 -E >> ppAFM.out
python3 $location/ppafm-generate-elff -i $3 --tip_dens $2 --Rcore 0.7 -E --doDensity >> ppAFM.out
python3 $location/ppafm-generate-dftd3 -i $3 --df_name PBE >> ppAFM.out
python3 $location/ppafm-relaxed-scan >> ppAFM.out
python3 $location/ppafm-plot-results --df --cbar --save_df >> ppAFM.out
find . -name *xsf -delete
}

function ppKPFM {
#argument after --Vref must be opposite sign to that of electric field in calculation, as per: https://github.com/Probe-Particle/ppafm/issues/334, https://github.com/Probe-Particle/ppafm/issues/250
python3 $location/ppafm-generate-elff -i cube_001_hartree_potential.cube -t dz2 --KPFM_sample cube_002_hartree_potential_0.1V.cube --KPFM_tip fit --Vref -0.1 --Rcore -1.0 --z0 0.0 >> ppAFM.out
python3 $location/ppafm-generate-ljff -i cube_001_hartree_potential.cube >> ppAFM.out

python3 $location/ppafm-relaxed-scan --Vrange -0.5 0.5 3 >> ppAFM.out
python3 $location/ppafm-plot-results --Vrange -0.5 0.5 3 --LCPD_maps --cbar --V0 0 --df --save_df >> ppAFM.out

find . -name *xsf -delete
#rm FFel_x.xsf FFel_y.xsf FFel_z.xsf FFkpfm_t0sV_x.xsf FFkpfm_t0sV_y.xsf FFkpfm_t0sV_z.xsf FFkpfm_tVs0_x.xsf FFkpfm_tVs0_y.xsf FFkpfm_tVs0_z.xsf FFLJ_x.xsf FFLJ_y.xsf FFLJ_z.xsf rhoTip.xsf Q*/OutFz.xsf Q*/A*/df.xsf

mv LCDP_HzperV.xsf LCPD_volts.xsf # this reflects correct units, as addressed by aureliojgc in ppAFM wiki (https://github.com/Probe-Particle/ppafm/pull/193)

mkdir PNGs
mv *png PNGs
}


date > ppAFM.out

#FDBM sample/cube_001_total_density.cube tip/CHGCAR.xsf sample/cube_004_total_potential.cube
#LJ_Hartree
LJ_PointCharges
#ppKPFM

date >> ppAFM.out
elapsed_time ppAFM.out

folder=$(pwd)
cd ..
mv $folder ~/ppAFM


