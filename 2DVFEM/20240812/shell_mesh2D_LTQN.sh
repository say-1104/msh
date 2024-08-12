#!/bin/sh

# 2D_VFEM用

FILE=`echo ${1} | sed -e "s/\.[^.]*$//"`
if [ $# -ne 1 ]; then
  echo "usage example : (for hogehoge.lps)"
  echo "  > mesh.sh hogehoge"
  exit 1
fi

if [ -f "${FILE}.lps" ]; then
  echo "[read ${FILE}.lps]"
else 
  echo "error : not exist \"${FILE}.lps\""
  exit 1
fi

saveMesh=`cat ${FILE}.lps | grep -n saveMesh | sed -e "s/:saveMesh.*//g"`
exportMESH=`cat ${FILE}.lps | grep -n exportMesh | sed -e "s/:exportMesh.*//g"`
#echo $exportMESH-$saveMesh
tmp=`echo ${exportMESH}-${saveMesh} | bc`
#echo tmp:$tmp
saveMesh=`cat ${FILE}.lps | sed -e "1,${saveMesh}d"`
saveMesh=`echo ${saveMesh} | head -n ${tmp} | sed -e 's/[^"]*"\([^"]*\)".*/\1/'`
exportMESH=`cat ${FILE}.lps | sed -e "1,${exportMESH}d" | sed -e 's/[^"]*"\([^"]*\)".*/\1/'`

echo saveMesh: $saveMesh.gid
echo exportMesh: $exportMESH.gid.msh
echo "========================="

#.lps to .bch
/home/okazaki/FEMrensyu/program/tools/lps2bch-3d/lps2bch-3d -kakihara -d 2 $FILE
echo "========================="

/home/share/Gid10.2.1/gid -n -b $FILE.bch
/home/okazaki/FEMrensyu/program/tools/gidmsh2msh_for_Gid10.2.1/ver.3.2/gidmsh2msh -detail -kakihara $exportMESH
/home/okazaki/FEMrensyu/program/tools/gid-de-mesh/VMesh/mod_kaki5/vfem-mesh -t 3 -detail  < $exportMESH.msh > $exportMESH.v2.msh


