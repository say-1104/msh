#!/bin/sh
dz=0.1
z0=0.0
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "コンパイルエラーが出たので終了します。"
	exit
fi
make all
for wl in `seq 1.550 0.010 1.550`
do
    echo -n > Leff0.dat
    for g in `seq 0.100 0.002 0.200`
    do
        echo -n > g${g}.dat
        z1=
        ./make_cfg ${dz} ${z0} ${g} ${z1} ${g}
        for Leff in `seq 0.0 0.2 26.6`
        do
            ./cmt -pcm1 -wl ${wl} -leff ${Leff} #2>output_err

            echo -e -n "${Leff}\t`tail -n1 output | head -n1 | cut -f3 `" >> ./g${g}.dat
            echo -e "\t`tail -n1 output | head -n1 | cut -f4 `" >> ./g${g}.dat
        done

    done
done