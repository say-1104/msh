#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "iteration.h"

#define ALLOCATION(a,c,n) \
  if ((a = (c *)calloc((n),sizeof(c))) == NULL) {		\
    fprintf(stderr, "can't allocate memory\n");			\
    exit(0);							\
  }

static long long int sorting(long long int *aa, long long int nn) {
	
  long long int i, j, k, l, m, w;
  
  l = nn / 2;
  for (k = l; k > 0; k /= 2) {
    for (i = k; i < nn; i++) {
      j = i - k;
      while (1) {
	m = j + k;
	if (j < 0 || aa[j] <= aa[m])
	  break;
	w = aa[j];
	aa[j] = aa[m];
	aa[m] = w;
	j -= k;
      }
    }
  }
  
  return 0;
}

/* 各行の非零要素の個数を調べ，行列格納時に格納場所を指定するテーブルを作る */
void makeCGmatrixTable(ELEMENT *element, long long int ne, long long int np, 
		       CGtable *cg_table, long long int NK) {
  long long int i, j, k, n, ii, jj;
  long long int imax;
  long long int *aa, **bb;
  long long int tmp;

  /* 各行に含まれる非零要素の個数をカウントする */
  /* 重複を無視した最大個数をカウントする */
  ALLOCATION(aa, long long int, np);
  for (i = 0; i < np; i++)
    aa[i] = 0;
  for (k = 0; k < ne; k++) {
    for (i = 0; i < NK; i++) {
      if (element[k].kk[i] < np && element[k].kk[i] >= 0) {
	aa[element[k].kk[i]] += NK;
      }
    }
  }
  
  /* 最大個数を格納できる配列を用意する */
  for (imax = 0, i = 0; i < np; i++) {
    if (aa[i] > imax)
      imax = aa[i];
  }
  tmp = np * imax;
  fprintf(stderr, "np * imax = %lld  ", np * imax);
  fprintf(stderr, "Calculated Memory %lld GiB \n", 
	  np * imax * sizeof(long long int) / (1024*1024*1024) );

  ALLOCATION(bb, long long int*, np);
  for(i=0;i<np;i++){
    ALLOCATION(bb[i], long long int, imax);
  }
    
  /* 重複を考慮し各行の非零要素の列番号を bb に格納する */
  /* aa には各行の現在の非零要素数を格納している */
  fprintf(stderr, "aa bb initializing ... \n");
  for (i = 0; i < np; i++){
    aa[i] = 0;
    for (j = 0; j < imax; j++)
      bb[i][j] = 0;
  }
  fprintf(stderr, "aa bb initialized ... \n");
  for (k = 0; k < ne; k++) {
    for (i = 0; i < NK; i++) {
      ii = element[k].kk[i];
      if (ii < np && ii >= 0) {
	for (j = 0; j < NK; j++) {
	  jj = element[k].kk[j];
	  if (jj < np && jj >= 0) {
	    for (n = 0; n < aa[ii]; n++) {
	      if (bb[ii][n] == jj)
		break;
	    }
	    /* 未登録ならば登録する */
	    if (n == aa[ii]) {
	      bb[ii][n] = jj;
	      aa[ii]++;
	    }
	  }
	}
      }
    }
  }
    
  /* 各行ごとに非零要素のある列番号を小さい順番に並べ替えて
     要素行列から全体行列への組み込みに用いるテーブルを作る */
  for (i = 0; i < np; i++) {
    sorting(bb[i], aa[i]);
    cg_table[i].n = aa[i];
    ALLOCATION(cg_table[i].col, long long int, cg_table[i].n);
    for (j = 0; j < cg_table[i].n; j++) {
      cg_table[i].col[j] = bb[i][j];
    }
  }
  
  free(aa);
  for(i = 0; i < np; i++) free(bb[i]);
  free(bb);
  fprintf(stderr, "free aa bb \n");
  return;
}

/* 行列の列番号に対して圧縮形式での格納場所を返す */
long long int columnInCGmatrix(long long int ii, long long int jj, 
			       CGtable *cg_table, long long int np) {
  long long int i;
  
  /* ディリクレ境界等の場合は列番号 -1 を返し格納しない */
  if (ii >= np)
    return -1;
  
  /* 指定行を探索し非零要素の格納場所を調べる */
  for (i = 0; i < cg_table[ii].n; i++) {
    if (cg_table[ii].col[i] == jj)
      return i;
  }
  
  return -1;
}
