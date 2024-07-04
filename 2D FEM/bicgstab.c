/*
	A[n][?] : 
	cg_table[n] :
*/
#include "optfem.h"

static int sorting(int *aa, int nn)
{
	int i, j, k, l, m, w;

	l = nn/2;
	for (k = l; k > 0; k /= 2) {
		for (i = k; i < nn; i++) {
			j = i-k;
			while(1) {
				m = j+k;
				if (j < 0 || aa[j] <= aa[m]) break;
				w = aa[j]; aa[j] = aa[m]; aa[m] = w;
				j -= k;
			}
		}
	}

	return 0;
}

// 各行の非零要素の個数を調べ，行列格納時に格納場所を指定するテーブルを作る
extern int makeCGmatrixTable(Element3D *element, int ne, int np,
							CGtable *cg_table, int NK)
{
	int i, j, k, n, ii, jj;
	int imax;
	int *aa, *bb;

	/* 各行に含まれる非零要素の個数をカウントする */
	/* 重複を無視した最大個数をカウントする */
	ALLOCATION(aa, int, np)
	for (i = 0; i < np; i++) aa[i] = 0;
	for (k = 0; k < ne; k++) {
		for (i = 0; i < NK; i++) {
			if (element[k].kk[i] < np && element[k].kk[i] >= 0) {
				aa[element[k].kk[i]] += NK;
			}
		}
	}

	/* 最大個数を格納できる配列を用意する */
	for (imax = 0, i = 0; i < np; i++) {
		if (aa[i] > imax) imax = aa[i];
	}
	ALLOCATION(bb, int, np*imax)

	/* 重複を考慮し各行の非零要素の列番号を bb に格納する */
	/* aa には各行の現在の非零要素数を格納している */
	for (i = 0; i < np; i++) aa[i] = 0;
	for (i = 0;i < np*imax; i++) bb[i] = 0;
	for (k = 0; k < ne; k++) {
		for (i = 0; i < NK; i++) {
			ii = element[k].kk[i];
			if (ii < np && ii >= 0) {
				for (j = 0; j < NK; j++) {
					jj = element[k].kk[j];
					if (jj < np && jj >= 0) {
						for (n = 0; n < aa[ii]; n++) {
							if (bb[imax*ii+n] == jj) break;
						}
						/* 未登録ならば登録する */
						if (n == aa[ii]) {
							bb[imax*ii+n] = jj;
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
		sorting(bb+i*imax, aa[i]);
		cg_table[i].n = aa[i];
		ALLOCATION(cg_table[i].col, int, cg_table[i].n)
		for (j = 0; j < cg_table[i].n; j++) {
			cg_table[i].col[j] = bb[i*imax+j];
		}
	}

	free(aa); free(bb);

	return 0;
}

// 行列の列番号に対して圧縮形式での格納場所を返す
extern int columnInCGmatrix(int ii, int jj, CGtable *cg_table, int np)
{
	int i;

	/* ディリクレ境界等の場合は列番号 -1 を返し格納しない */
	if (ii >= np) return -1;

	/* 指定行を探索し非零要素の格納場所を調べる */
	for (i = 0; i < cg_table[ii].n; i++) {
		if (cg_table[ii].col[i] == jj) return i;
	}

	return -1;
}