#!/bin/sh
export OMP_NUM_THREADS=24
export MKL_NUM_THREADS=24

g++ LpsTools.cpp make_lps.cpp -I./  -std=c++11 

if [ ! -e ./a.out ]; then
	echo "コンパイルエラーが出たので終了します。"
	exit
fi

#for wide in `seq 0.300 0.100 1.200`
for wide in 0.400
do
    ./a.out ${wide}
    sh shell_mesh2D_LTQN.sh wire
    #sh exec2D_VFEM.sh
done
