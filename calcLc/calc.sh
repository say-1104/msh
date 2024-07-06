#!/bin/sh
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12

icc -o make_lps make_lps.cpp

./make_lps 0.420 0.450 0.322 0.153 0.020
sh shell_mesh2D_LTQN.sh wire

sh exec_beo_a.sh
cp -r wave_1.530 ./beo_a
rm -r wave_1.530 
cp -r wave_1.540 ./beo_a
rm -r wave_1.540 
cp -r wave_1.560 ./beo_a
rm -r wave_1.560 
cp -r wave_1.570 ./beo_a
rm -r wave_1.570 

sh exec_beo_c.sh
cp -r wave_1.530 ./beo_c
rm -r wave_1.530 
cp -r wave_1.540 ./beo_c
rm -r wave_1.540 
cp -r wave_1.560 ./beo_c
rm -r wave_1.560 
cp -r wave_1.570 ./beo_c
rm -r wave_1.570 

sh exec_b1.sh
cp -r wave_1.530 ./b1
rm -r wave_1.530 
cp -r wave_1.540 ./b1
rm -r wave_1.540 
cp -r wave_1.560 ./b1
rm -r wave_1.560 
cp -r wave_1.570 ./b1
rm -r wave_1.570 

sh exec_b2_a.sh
cp -r wave_1.530 ./b2_a
rm -r wave_1.530
cp -r wave_1.540 ./b2_a
rm -r wave_1.540
cp -r wave_1.560 ./b2_a
rm -r wave_1.560
cp -r wave_1.570 ./b2_a
rm -r wave_1.570

sh exec_b2_c.sh
cp -r wave_1.530 ./b2_c
rm -r wave_1.530 
cp -r wave_1.540 ./b2_c
rm -r wave_1.540 
cp -r wave_1.560 ./b2_c
rm -r wave_1.560 
cp -r wave_1.570 ./b2_c
rm -r wave_1.570 
