#!/bin/sh
# arguments
# 1: name of mesh file

if [ $# -ne 1 ]; then
    echo "指定された引数が$#個です。" 1>&2
    echo "引数は、1:ファイル名" 1>&2
    exit 1
fi

rm -f ./a.out
rm -f $1*msh
rm -f $1.lps

g++ LpsTools.cpp make_lps.cpp -I./  -std=c++11 

if [ ! -e ./a.out ]; then
	echo "コンパイルエラーが出たので終了します。"
	exit
fi

./a.out 0.400

/home/okazaki/Solver/for3dmesh/lps2bch-3d_aug/lps2bch-3d -kakihara $1
/home/share/Gid10.2.1/gid -n -b $1.bch
/home/okazaki/Solver/for3dmesh/gidmsh2msh/gidmsh2msh -kakihara $1