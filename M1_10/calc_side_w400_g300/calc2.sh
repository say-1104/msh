#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD50 calc_MD50.cpp -std=c++11
if [ ! -e ./calc_MD50 ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.300
L=50.0
echo -n > forMD.dat
for g in `seq 0.350 0.010 0.350`
do
	echo -n > forMD50.dat
	./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 2

	for wl in `seq 1.530 0.010 1.570`
	do
		~/msh/CMT/cmt -pcm2 -wl ${wl}
		echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
	done
	./calc_MD0 <forMD50.dat
	echo "${g}\t`tail -n1 MDlambda | head -n1 | cut -f1 `" >> forMD.dat
done