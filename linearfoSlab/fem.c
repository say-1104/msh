/*
	Finite Element Method for 1-D Slab Optical Waveguides
	Produced by Dr. Takeshi FUJISAWA
	Dec 2007

��Ĺ�򿶤ä������������뤳�Ȥ��Ǥ��롥
��Ĺʬ����Ʊ���˷׻����롥

�ޥ���ե��󥿥��Ȥ��Ȥ��ޤ������ʤ���
����Ϥ�����ȤǤ��Ƥ��뤬propgation���̤���
���������ʤ롥����������

usage

1) sfem a.msh -a%X 
2) sfem a.msh -e 

test.msh
*/

#include "fem.h"

extern int main(int argc, char **argv)
{
  int	i, j, nf, count = 0, flag;
  double   xc_max, xc_temp, delta;
  DataTable	data;
  char string[256];
  FILE *fp;
    
  analysis(&data, argc, argv);

}
