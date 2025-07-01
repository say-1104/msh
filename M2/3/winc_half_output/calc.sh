#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi

#DIR=$(cd $(dirname $0);pwd)
w=0.400
p=0.300
d=0.000

dz=0.01
g1=`echo "scale=3;0.200 - ${w} + 0.400" | bc`
g2=`echo "scale=3;1.200 - ${w} + 0.400" | bc`
Lsbend=10.0
filename="output_w${w}p${p}d${d}"
echo -n > ${filename}.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> ${filename}.dat
L1=14.4
L2=8.0
./make_cfg ${dz} ${g2} ${w} 4 sbend ${Lsbend} ${g1} 3 dc ${L1} 8 dc ${L2} 7 sbend ${Lsbend} ${g2} 3

for wl in `seq 1.530 0.010 1.570`
do
    ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
    cp output ./output_${wl}
done
i=2
for L in `seq 0.0 0.01 42.4`
do  
    echo -n "${L}" >> ${filename}.dat
    for wl in `seq 1.530 0.010 1.570`
    do
        echo -n "\t`head -n${i} output_${wl} | tail -n1 | cut -f3 `" >> ${filename}.dat
    done
    echo >> ${filename}.dat
    i=`expr $i + 1`
done