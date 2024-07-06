#!/bin/sh
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12
wl=1.550
icc -o make_lps make_lps.cpp
echo -n > Lc_${wl}.dat
for g in `seq 0.100 0.010 0.300`
do
    ./make_lps 0.420 0.450 0.322 $g 0.020
    sh shell_mesh2D_LTQN.sh wire
    sh exec_beo_c.sh

    cd wave_${wl}
    echo -e -n "${g}\t`tail -n1 _results.dat | head -n1 | cut -f2 `" >> ../Lc_${wl}.dat
    echo -e "\t`tail -n2 _results.dat | head -n1 | cut -f2 `" >> ../Lc_${wl}.dat
    cd ..

done