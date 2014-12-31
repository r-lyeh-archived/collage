#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#define new new_

//#line 1 "bsdiff.c"
#include <stddef.h>
#include <stdint.h>

struct bsdiff_stream
{
	void* opaque;

	void* (*malloc)(size_t size);
	void (*free)(void* ptr);
	int (*write)(struct bsdiff_stream* stream, const void* buffer, int size);
};

int bsdiff(const uint8_t* old, int64_t oldsize, const uint8_t* new, int64_t newsize, struct bsdiff_stream* stream);

#if !defined(BSDIFF_HEADER_ONLY)

#include <limits.h>
#include <string.h>

#define MIN(x,y) (((x)<(y)) ? (x) : (y))

static void split(int64_t *I,int64_t *V,int64_t start,int64_t len,int64_t h)
{
	int64_t i,j,k,x,tmp,jj,kk;

	if(len<16) {
		for(k=start;k<start+len;k+=j) {
			j=1;x=V[I[k]+h];
			for(i=1;k+i<start+len;i++) {
				if(V[I[k+i]+h]<x) {
					x=V[I[k+i]+h];
					j=0;
				};
				if(V[I[k+i]+h]==x) {
					tmp=I[k+j];I[k+j]=I[k+i];I[k+i]=tmp;
					j++;
				};
			};
			for(i=0;i<j;i++) V[I[k+i]]=k+j-1;
			if(j==1) I[k]=-1;
		};
		return;
	};

	x=V[I[start+len/2]+h];
	jj=0;kk=0;
	for(i=start;i<start+len;i++) {
		if(V[I[i]+h]<x) jj++;
		if(V[I[i]+h]==x) kk++;
	};
	jj+=start;kk+=jj;

	i=start;j=0;k=0;
	while(i<jj) {
		if(V[I[i]+h]<x) {
			i++;
		} else if(V[I[i]+h]==x) {
			tmp=I[i];I[i]=I[jj+j];I[jj+j]=tmp;
			j++;
		} else {
			tmp=I[i];I[i]=I[kk+k];I[kk+k]=tmp;
			k++;
		};
	};

	while(jj+j<kk) {
		if(V[I[jj+j]+h]==x) {
			j++;
		} else {
			tmp=I[jj+j];I[jj+j]=I[kk+k];I[kk+k]=tmp;
			k++;
		};
	};

	if(jj>start) split(I,V,start,jj-start,h);

	for(i=0;i<kk-jj;i++) V[I[jj+i]]=kk-1;
	if(jj==kk-1) I[jj]=-1;

	if(start+len>kk) split(I,V,kk,start+len-kk,h);
}

/* [ref] Code commented because of http://stackoverflow.com/questions/12751775/why-does-bsdiff-exe-have-trouble-with-this-smaller-file */
/* [ref] See Graeme Johnson's answer */
#if 1//def BSDIFF_USE_SAIS

//#line 1 "sais.h"
#ifndef _SAIS_H
#define _SAIS_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* find the suffix array SA of T[0..n-1]
   use a working space (excluding T and SA) of at most 2n+O(lg n) */
int
sais(const unsigned char *T, int *SA, int n);

/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
   use a working space (excluding T and SA) of at most MAX(4k,2n) */
int
sais_int(const int *T, int *SA, int n, int k);

/* burrows-wheeler transform */
int
sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n);
int
sais_int_bwt(const int *T, int *U, int *A, int n, int k);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _SAIS_H */



//#line 1 "sais.c"
/*
 * sais.c for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <stdlib.h>

//#line 1 "sais.h"
#ifndef _SAIS_H
#define _SAIS_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* find the suffix array SA of T[0..n-1]
   use a working space (excluding T and SA) of at most 2n+O(lg n) */
int
sais(const unsigned char *T, int *SA, int n);

/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
   use a working space (excluding T and SA) of at most MAX(4k,2n) */
int
sais_int(const int *T, int *SA, int n, int k);

/* burrows-wheeler transform */
int
sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n);
int
sais_int_bwt(const int *T, int *U, int *A, int n, int k);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _SAIS_H */


#ifndef UCHAR_SIZE
# define UCHAR_SIZE 256
#endif
#ifndef MINBUCKETSIZE
# define MINBUCKETSIZE 256
#endif

#define sais_index_type int
#define sais_bool_type  int
#define SAIS_LMSSORT2_LIMIT 0x3fffffff

#define SAIS_MYMALLOC(_num, _type) ((_type *)malloc((_num) * sizeof(_type)))
#define SAIS_MYFREE(_ptr, _num, _type) free((_ptr))
#define chr(_a) (cs == sizeof(sais_index_type) ? ((sais_index_type *)T)[(_a)] : ((unsigned char *)T)[(_a)])

/* find the start or end of each bucket */
static
void
getCounts(const void *T, sais_index_type *C, sais_index_type n, sais_index_type k, int cs) {
  sais_index_type i;
  for(i = 0; i < k; ++i) { C[i] = 0; }
  for(i = 0; i < n; ++i) { ++C[chr(i)]; }
}
static
void
getBuckets(const sais_index_type *C, sais_index_type *B, sais_index_type k, sais_bool_type end) {
  sais_index_type i, sum = 0;
  if(end) { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum; } }
  else { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum - C[i]; } }
}

/* sort all type LMS suffixes */
static
void
LMSsort1(const void *T, sais_index_type *SA,
		 sais_index_type *C, sais_index_type *B,
		 sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j;
  sais_index_type c0, c1;

  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  --j;
  *b++ = (chr(j) < c1) ? ~j : j;
  for(i = 0; i < n; ++i) {
	if(0 < (j = SA[i])) {
	  assert(chr(j) >= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert(i < (b - SA));
	  --j;
	  *b++ = (chr(j) < c1) ? ~j : j;
	  SA[i] = 0;
	} else if(j < 0) {
	  SA[i] = ~j;
	}
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
	if(0 < (j = SA[i])) {
	  assert(chr(j) <= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert((b - SA) <= i);
	  --j;
	  *--b = (chr(j) > c1) ? ~(j + 1) : j;
	  SA[i] = 0;
	}
  }
}
static
sais_index_type
LMSpostproc1(const void *T, sais_index_type *SA,
			 sais_index_type n, sais_index_type m, int cs) {
  sais_index_type i, j, p, q, plen, qlen, name;
  sais_index_type c0, c1;
  sais_bool_type diff;

  /* compact all the sorted substrings into the first m items of SA
	  2*m must be not larger than n (proveable) */
  assert(0 < n);
  for(i = 0; (p = SA[i]) < 0; ++i) { SA[i] = ~p; assert((i + 1) < n); }
  if(i < m) {
	for(j = i, ++i;; ++i) {
	  assert(i < n);
	  if((p = SA[i]) < 0) {
		SA[j++] = ~p; SA[i] = 0;
		if(j == m) { break; }
	  }
	}
  }

  /* store the length of all substrings */
  i = n - 1; j = n - 1; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
	do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
	if(0 <= i) {
	  SA[m + ((i + 1) >> 1)] = j - i; j = i + 1;
	  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
	}
  }

  /* find the lexicographic names of all substrings */
  for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
	p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
	if((plen == qlen) && ((q + plen) < n)) {
	  for(j = 0; (j < plen) && (chr(p + j) == chr(q + j)); ++j) { }
	  if(j == plen) { diff = 0; }
	}
	if(diff != 0) { ++name, q = p, qlen = plen; }
	SA[m + (p >> 1)] = name;
  }

  return name;
}
static
void
LMSsort2(const void *T, sais_index_type *SA,
		 sais_index_type *C, sais_index_type *B, sais_index_type *D,
		 sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j, t, d;
  sais_index_type c0, c1;
  assert(C != B);

  /* compute SAl */
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  --j;
  t = (chr(j) < c1);
  j += n;
  *b++ = (t & 1) ? ~j : j;
  for(i = 0, d = 0; i < n; ++i) {
	if(0 < (j = SA[i])) {
	  if(n <= j) { d += 1; j -= n; }
	  assert(chr(j) >= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert(i < (b - SA));
	  --j;
	  t = c0; t = (t << 1) | (chr(j) < c1);
	  if(D[t] != d) { j += n; D[t] = d; }
	  *b++ = (t & 1) ? ~j : j;
	  SA[i] = 0;
	} else if(j < 0) {
	  SA[i] = ~j;
	}
  }
  for(i = n - 1; 0 <= i; --i) {
	if(0 < SA[i]) {
	  if(SA[i] < n) {
		SA[i] += n;
		for(j = i - 1; SA[j] < n; --j) { }
		SA[j] -= n;
		i = j;
	  }
	}
  }

  /* compute SAs */
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, d += 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
	if(0 < (j = SA[i])) {
	  if(n <= j) { d += 1; j -= n; }
	  assert(chr(j) <= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert((b - SA) <= i);
	  --j;
	  t = c0; t = (t << 1) | (chr(j) > c1);
	  if(D[t] != d) { j += n; D[t] = d; }
	  *--b = (t & 1) ? ~(j + 1) : j;
	  SA[i] = 0;
	}
  }
}
static
sais_index_type
LMSpostproc2(sais_index_type *SA, sais_index_type n, sais_index_type m) {
  sais_index_type i, j, d, name;

  /* compact all the sorted LMS substrings into the first m items of SA */
  assert(0 < n);
  for(i = 0, name = 0; (j = SA[i]) < 0; ++i) {
	j = ~j;
	if(n <= j) { name += 1; }
	SA[i] = j;
	assert((i + 1) < n);
  }
  if(i < m) {
	for(d = i, ++i;; ++i) {
	  assert(i < n);
	  if((j = SA[i]) < 0) {
		j = ~j;
		if(n <= j) { name += 1; }
		SA[d++] = j; SA[i] = 0;
		if(d == m) { break; }
	  }
	}
  }
  if(name < m) {
	/* store the lexicographic names */
	for(i = m - 1, d = name + 1; 0 <= i; --i) {
	  if(n <= (j = SA[i])) { j -= n; --d; }
	  SA[m + (j >> 1)] = d;
	}
  } else {
	/* unset flags */
	for(i = 0; i < m; ++i) {
	  if(n <= (j = SA[i])) { j -= n; SA[i] = j; }
	}
  }

  return name;
}

/* compute SA and BWT */
static
void
induceSA(const void *T, sais_index_type *SA,
		 sais_index_type *C, sais_index_type *B,
		 sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j;
  sais_index_type c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
	j = SA[i], SA[i] = ~j;
	if(0 < j) {
	  --j;
	  assert(chr(j) >= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert(i < (b - SA));
	  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
	}
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
	if(0 < (j = SA[i])) {
	  --j;
	  assert(chr(j) <= chr(j + 1));
	  if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert((b - SA) <= i);
	  *--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
	} else {
	  SA[i] = ~j;
	}
  }
}
static
sais_index_type
computeBWT(const void *T, sais_index_type *SA,
		   sais_index_type *C, sais_index_type *B,
		   sais_index_type n, sais_index_type k, int cs) {
  sais_index_type *b, i, j, pidx = -1;
  sais_index_type c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
	if(0 < (j = SA[i])) {
	  --j;
	  assert(chr(j) >= chr(j + 1));
	  SA[i] = ~((sais_index_type)(c0 = chr(j)));
	  if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert(i < (b - SA));
	  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
	} else if(j != 0) {
	  SA[i] = ~j;
	}
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
	if(0 < (j = SA[i])) {
	  --j;
	  assert(chr(j) <= chr(j + 1));
	  SA[i] = (c0 = chr(j));
	  if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
	  assert((b - SA) <= i);
	  *--b = ((0 < j) && (chr(j - 1) > c1)) ? ~((sais_index_type)chr(j - 1)) : j;
	} else if(j != 0) {
	  SA[i] = ~j;
	} else {
	  pidx = i;
	}
  }
  return pidx;
}

/* find the suffix array SA of T[0..n-1] in {0..255}^n */
static
sais_index_type
sais_main(const void *T, sais_index_type *SA,
		  sais_index_type fs, sais_index_type n, sais_index_type k, int cs,
		  sais_bool_type isbwt) {
  sais_index_type *C, *B, *D, *RA, *b;
  sais_index_type i, j, m, p, q, t, name, pidx = 0, newfs;
  sais_index_type c0, c1;
  unsigned int flags;

  assert((T != NULL) && (SA != NULL));
  assert((0 <= fs) && (0 < n) && (1 <= k));

  if(k <= MINBUCKETSIZE) {
	if((C = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
	if(k <= fs) {
	  B = SA + (n + fs - k);
	  flags = 1;
	} else {
	  if((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { SAIS_MYFREE(C, k, sais_index_type); return -2; }
	  flags = 3;
	}
  } else if(k <= fs) {
	C = SA + (n + fs - k);
	if(k <= (fs - k)) {
	  B = C - k;
	  flags = 0;
	} else if(k <= (MINBUCKETSIZE * 4)) {
	  if((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
	  flags = 2;
	} else {
	  B = C;
	  flags = 8;
	}
  } else {
	if((C = B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) { return -2; }
	flags = 4 | 8;
  }
  if((n <= SAIS_LMSSORT2_LIMIT) && (2 <= (n / k))) {
	if(flags & 1) { flags |= ((k * 2) <= (fs - k)) ? 32 : 16; }
	else if((flags == 0) && ((k * 2) <= (fs - k * 2))) { flags |= 32; }
  }

  /* stage 1: reduce the problem by at least 1/2
	 sort all the LMS-substrings */
  getCounts(T, C, n, k, cs); getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = 0; i < n; ++i) { SA[i] = 0; }
  b = &t; i = n - 1; j = n; m = 0; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
	do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
	if(0 <= i) {
	  *b = j; b = SA + --B[c1]; j = i; ++m;
	  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
	}
  }

  if(1 < m) {
	if(flags & (16 | 32)) {
	  if(flags & 16) {
		if((D = SAIS_MYMALLOC(k * 2, sais_index_type)) == NULL) {
		  if(flags & (1 | 4)) { SAIS_MYFREE(C, k, sais_index_type); }
		  if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }
		  return -2;
		}
	  } else {
		D = B - k * 2;
	  }
	  assert((j + 1) < n);
	  ++B[chr(j + 1)];
	  for(i = 0, j = 0; i < k; ++i) {
		j += C[i];
		if(B[i] != j) { assert(SA[B[i]] != 0); SA[B[i]] += n; }
		D[i] = D[i + k] = 0;
	  }
	  LMSsort2(T, SA, C, B, D, n, k, cs);
	  name = LMSpostproc2(SA, n, m);
	  if(flags & 16) { SAIS_MYFREE(D, k * 2, sais_index_type); }
	} else {
	  LMSsort1(T, SA, C, B, n, k, cs);
	  name = LMSpostproc1(T, SA, n, m, cs);
	}
  } else if(m == 1) {
	*b = j + 1;
	name = 1;
  } else {
	name = 0;
  }

  /* stage 2: solve the reduced problem
	 recurse if names are not yet unique */
  if(name < m) {
	if(flags & 4) { SAIS_MYFREE(C, k, sais_index_type); }
	if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }
	newfs = (n + fs) - (m * 2);
	if((flags & (1 | 4 | 8)) == 0) {
	  if((k + name) <= newfs) { newfs -= k; }
	  else { flags |= 8; }
	}
	assert((n >> 1) <= (newfs + m));
	RA = SA + m + newfs;
	for(i = m + (n >> 1) - 1, j = m - 1; m <= i; --i) {
	  if(SA[i] != 0) {
		RA[j--] = SA[i] - 1;
	  }
	}
	if(sais_main(RA, SA, newfs, m, name, sizeof(sais_index_type), 0) != 0) {
	  if(flags & 1) { SAIS_MYFREE(C, k, sais_index_type); }
	  return -2;
	}

	i = n - 1; j = m - 1; c0 = chr(n - 1);
	do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
	for(; 0 <= i;) {
	  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
	  if(0 <= i) {
		RA[j--] = i + 1;
		do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
	  }
	}
	for(i = 0; i < m; ++i) { SA[i] = RA[SA[i]]; }
	if(flags & 4) {
	  if((C = B = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
	}
	if(flags & 2) {
	  if((B = SAIS_MYMALLOC(k, int)) == NULL) {
		if(flags & 1) { SAIS_MYFREE(C, k, sais_index_type); }
		return -2;
	  }
	}
  }

  /* stage 3: induce the result for the original problem */
  if(flags & 8) { getCounts(T, C, n, k, cs); }
  /* put all left-most S characters into their buckets */
  if(1 < m) {
	getBuckets(C, B, k, 1); /* find ends of buckets */
	i = m - 1, j = n, p = SA[m - 1], c1 = chr(p);
	do {
	  q = B[c0 = c1];
	  while(q < j) { SA[--j] = 0; }
	  do {
		SA[--j] = p;
		if(--i < 0) { break; }
		p = SA[i];
	  } while((c1 = chr(p)) == c0);
	} while(0 <= i);
	while(0 < j) { SA[--j] = 0; }
  }
  if(isbwt == 0) { induceSA(T, SA, C, B, n, k, cs); }
  else { pidx = computeBWT(T, SA, C, B, n, k, cs); }
  if(flags & (1 | 4)) { SAIS_MYFREE(C, k, sais_index_type); }
  if(flags & 2) { SAIS_MYFREE(B, k, sais_index_type); }

  return pidx;
}

/*---------------------------------------------------------------------------*/

int
sais(const unsigned char *T, int *SA, int n) {
  if((T == NULL) || (SA == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; } return 0; }
  return sais_main(T, SA, 0, n, UCHAR_SIZE, sizeof(unsigned char), 0);
}

int
sais_int(const int *T, int *SA, int n, int k) {
  if((T == NULL) || (SA == NULL) || (n < 0) || (k <= 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; } return 0; }
  return sais_main(T, SA, 0, n, k, sizeof(int), 0);
}

int
sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main(T, A, 0, n, UCHAR_SIZE, sizeof(unsigned char), 1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = (unsigned char)A[i]; }
  for(i += 1; i < n; ++i) { U[i] = (unsigned char)A[i]; }
  pidx += 1;
  return pidx;
}

int
sais_int_bwt(const int *T, int *U, int *A, int n, int k) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (k <= 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main(T, A, 0, n, k, sizeof(int), 1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = A[i]; }
  for(i += 1; i < n; ++i) { U[i] = A[i]; }
  pidx += 1;
  return pidx;
}

#endif

static void qsufsort(int64_t *I,int64_t *V,const uint8_t *old,int64_t oldsize)
{
	int64_t buckets[256];
	int64_t i,h,len;

	for(i=0;i<256;i++) buckets[i]=0;
	for(i=0;i<oldsize;i++) buckets[old[i]]++;
	for(i=1;i<256;i++) buckets[i]+=buckets[i-1];
	for(i=255;i>0;i--) buckets[i]=buckets[i-1];
	buckets[0]=0;

	for(i=0;i<oldsize;i++) I[++buckets[old[i]]]=i;
#if 1//def BSDIFF_USE_SAIS
	/* Graeme Johnson's solution */
	I[0] = oldsize; sais(old, ((int *)I)+1, oldsize);	return;
#else
	I[0] = oldsize;
#endif
	for(i=0;i<oldsize;i++) V[i]=buckets[old[i]];
	V[oldsize]=0;
	for(i=1;i<256;i++) if(buckets[i]==buckets[i-1]+1) I[buckets[i]]=-1;
	I[0]=-1;

	for(h=1;I[0]!=-(oldsize+1);h+=h) {
		len=0;
		for(i=0;i<oldsize+1;) {
			if(I[i]<0) {
				len-=I[i];
				i-=I[i];
			} else {
				if(len) I[i-len]=-len;
				len=V[I[i]]+1-i;
				split(I,V,i,len,h);
				i+=len;
				len=0;
			};
		};
		if(len) I[i-len]=-len;
	};

	for(i=0;i<oldsize+1;i++) I[V[i]]=i;
}

static int64_t matchlen(const uint8_t *old,int64_t oldsize,const uint8_t *new,int64_t newsize)
{
	int64_t i;

	for(i=0;(i<oldsize)&&(i<newsize);i++)
		if(old[i]!=new[i]) break;

	return i;
}

static int64_t search(const int64_t *I,const uint8_t *old,int64_t oldsize,
		const uint8_t *new,int64_t newsize,int64_t st,int64_t en,int64_t *pos)
{
	int64_t x,y;

	if(en-st<2) {
		x=matchlen(old+I[st],oldsize-I[st],new,newsize);
		y=matchlen(old+I[en],oldsize-I[en],new,newsize);

		if(x>y) {
			*pos=I[st];
			return x;
		} else {
			*pos=I[en];
			return y;
		}
	};

	x=st+(en-st)/2;
	if(memcmp(old+I[x],new,MIN(oldsize-I[x],newsize))<0) {
		return search(I,old,oldsize,new,newsize,x,en,pos);
	} else {
		return search(I,old,oldsize,new,newsize,st,x,pos);
	};
}

static int offtout( int64_t ii, uint8_t *buf ) {
	/* taken from https://github.com/r-lyeh/vle */
	uint64_t i = (uint64_t)ii;
	i = i & (1ull << 63) ? ~(i << 1) : (i << 1);
	unsigned char *origin = buf;
	do {
		*buf++ = (unsigned char)( 0x80 | (i & 0x7f));
		i >>= 7;
	} while( i > 0 );
	*(buf-1) ^= 0x80;
	return buf - origin;
}

static int64_t writedata(struct bsdiff_stream* stream, const void* buffer, int64_t length)
{
	int64_t result = 0;

	while (length > 0)
	{
		const int smallsize = (int)MIN(length, INT_MAX);
		const int writeresult = stream->write(stream, buffer, smallsize);
		if (writeresult == -1)
		{
			return -1;
		}

		result += writeresult;
		length -= smallsize;
		buffer = (uint8_t*)buffer + smallsize;
	}

	return result;
}

struct bsdiff_request
{
	const uint8_t* old;
	int64_t oldsize;
	const uint8_t* new;
	int64_t newsize;
	struct bsdiff_stream* stream;
	int64_t *I;
	uint8_t *buffer;
};

static int bsdiff_internal(const struct bsdiff_request req)
{
	int64_t *I,*V;
	int64_t scan,pos,len;
	int64_t lastscan,lastpos,lastoffset;
	int64_t oldscore,scsc;
	int64_t s,Sf,lenf,Sb,lenb;
	int64_t overlap,Ss,lens;
	int64_t i;
	uint8_t *buffer;
	uint8_t buf[10 * 3], *ptr;

	if((V=(int64_t*)req.stream->malloc((req.oldsize+1)*sizeof(int64_t)))==NULL) return -1;
	I = req.I;

	qsufsort(I,V,req.old,req.oldsize);
	req.stream->free(V);

	buffer = req.buffer;

	/* Compute the differences, writing ctrl as we go */
	scan=0;len=0;pos=0;
	lastscan=0;lastpos=0;lastoffset=0;
	while(scan<req.newsize) {
		oldscore=0;

		for(scsc=scan+=len;scan<req.newsize;scan++) {
			len=search(I,req.old,req.oldsize,req.new+scan,req.newsize-scan,
					0,req.oldsize,&pos);

			for(;scsc<scan+len;scsc++)
			if((scsc+lastoffset<req.oldsize) &&
				(req.old[scsc+lastoffset] == req.new[scsc]))
				oldscore++;

			if(((len==oldscore) && (len!=0)) ||
				(len>oldscore+8)) break;

			if((scan+lastoffset<req.oldsize) &&
				(req.old[scan+lastoffset] == req.new[scan]))
				oldscore--;
		};

		if((len!=oldscore) || (scan==req.newsize)) {
			s=0;Sf=0;lenf=0;
			for(i=0;(lastscan+i<scan)&&(lastpos+i<req.oldsize);) {
				if(req.old[lastpos+i]==req.new[lastscan+i]) s++;
				i++;
				if(s*2-i>Sf*2-lenf) { Sf=s; lenf=i; };
			};

			lenb=0;
			if(scan<req.newsize) {
				s=0;Sb=0;
				for(i=1;(scan>=lastscan+i)&&(pos>=i);i++) {
					if(req.old[pos-i]==req.new[scan-i]) s++;
					if(s*2-i>Sb*2-lenb) { Sb=s; lenb=i; };
				};
			};

			if(lastscan+lenf>scan-lenb) {
				overlap=(lastscan+lenf)-(scan-lenb);
				s=0;Ss=0;lens=0;
				for(i=0;i<overlap;i++) {
					if(req.new[lastscan+lenf-overlap+i]==
					   req.old[lastpos+lenf-overlap+i]) s++;
					if(req.new[scan-lenb+i]==
					   req.old[pos-lenb+i]) s--;
					if(s>Ss) { Ss=s; lens=i+1; };
				};

				lenf+=lens-overlap;
				lenb-=lens;
			};

			ptr = buf;
			ptr += offtout(lenf,buf);
			ptr += offtout((scan-lenb)-(lastscan+lenf),ptr);
			ptr += offtout((pos-lenb)-(lastpos+lenf),ptr);

			/* Write control data */
			if (writedata(req.stream, buf, ptr - buf))
				return -1;

			/* Write diff data */
			for(i=0;i<lenf;i++)
				buffer[i]=req.new[lastscan+i]-req.old[lastpos+i];
			if (writedata(req.stream, buffer, lenf))
				return -1;

			/* Write extra data */
			for(i=0;i<(scan-lenb)-(lastscan+lenf);i++)
				buffer[i]=req.new[lastscan+lenf+i];
			if (writedata(req.stream, buffer, (scan-lenb)-(lastscan+lenf)))
				return -1;

			lastscan=scan-lenb;
			lastpos=pos-lenb;
			lastoffset=pos-scan;
		};
	};

	return 0;
}

int bsdiff(const uint8_t* old, int64_t oldsize, const uint8_t* new, int64_t newsize, struct bsdiff_stream* stream)
{
	int result;
	struct bsdiff_request req;

	if((req.I=(int64_t*)stream->malloc((oldsize+1)*sizeof(int64_t)))==NULL)
		return -1;

	if((req.buffer=(uint8_t*)stream->malloc(newsize+1))==NULL)
	{
		stream->free(req.I);
		return -1;
	}

	req.old = old;
	req.oldsize = oldsize;
	req.new = new;
	req.newsize = newsize;
	req.stream = stream;

	result = bsdiff_internal(req);

	stream->free(req.buffer);
	stream->free(req.I);

	return result;
}

#if defined(BSDIFF_EXECUTABLE)

#include <sys/types.h>

#include <bzlib.h>
#include <err.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static int bz2_write(struct bsdiff_stream* stream, const void* buffer, int size)
{
	int bz2err;
	BZFILE* bz2;

	bz2 = (BZFILE*)stream->opaque;
	BZ2_bzWrite(&bz2err, bz2, (void*)buffer, size);
	if (bz2err != BZ_STREAM_END && bz2err != BZ_OK)
		return -1;

	return 0;
}

int main(int argc,char *argv[])
{
	int fd;
	int bz2err;
	uint8_t *old,*new;
	off_t oldsize,newsize;
	uint8_t buf[8];
	FILE * pf;
	struct bsdiff_stream stream;
	BZFILE* bz2;

	memset(&bz2, 0, sizeof(bz2));
	stream.malloc = malloc;
	stream.free = free;
	stream.write = bz2_write;

	if(argc!=4) errx(1,"usage: %s oldfile newfile patchfile\n",argv[0]);

	/* Allocate oldsize+1 bytes instead of oldsize bytes to ensure
		that we never try to malloc(0) and get a NULL pointer */
	if(((fd=open(argv[1],O_RDONLY,0))<0) ||
		((oldsize=lseek(fd,0,SEEK_END))==-1) ||
		((old=malloc(oldsize+1))==NULL) ||
		(lseek(fd,0,SEEK_SET)!=0) ||
		(read(fd,old,oldsize)!=oldsize) ||
		(close(fd)==-1)) err(1,"%s",argv[1]);

	/* Allocate newsize+1 bytes instead of newsize bytes to ensure
		that we never try to malloc(0) and get a NULL pointer */
	if(((fd=open(argv[2],O_RDONLY,0))<0) ||
		((newsize=lseek(fd,0,SEEK_END))==-1) ||
		((new=malloc(newsize+1))==NULL) ||
		(lseek(fd,0,SEEK_SET)!=0) ||
		(read(fd,new,newsize)!=newsize) ||
		(close(fd)==-1)) err(1,"%s",argv[2]);

	/* Create the patch file */
	if ((pf = fopen(argv[3], "w")) == NULL)
		err(1, "%s", argv[3]);

	/* Write header (signature+newsize)*/
	offtout(newsize, buf);
	if (fwrite("ENDSLEY/BSDIFF43", 16, 1, pf) != 1 ||
		fwrite(buf, sizeof(buf), 1, pf) != 1)
		err(1, "Failed to write header");

	if (NULL == (bz2 = BZ2_bzWriteOpen(&bz2err, pf, 9, 0, 0)))
		errx(1, "BZ2_bzWriteOpen, bz2err=%d", bz2err);

	stream.opaque = bz2;
	if (bsdiff(old, oldsize, new, newsize, &stream))
		err(1, "bsdiff");

	BZ2_bzWriteClose(&bz2err, bz2, 0, NULL, NULL);
	if (bz2err != BZ_OK)
		err(1, "BZ2_bzWriteClose, bz2err=%d", bz2err);

	if (fclose(pf))
		err(1, "fclose");

	/* Free the memory we used */
	free(old);
	free(new);

	return 0;
}

#endif

#endif



//#line 1 "bspatch.c"
#include <stdint.h>

struct bspatch_stream
{
	void* opaque;
	int (*read)(const struct bspatch_stream* stream, void* buffer, int length);
};

int bspatch(const uint8_t* old, int64_t oldsize, uint8_t* new, int64_t newsize, struct bspatch_stream* stream);

#if !defined(BSPATCH_HEADER_ONLY)

static int64_t offtin( const uint8_t *buf ) {
	/* taken from https://github.com/r-lyeh/vle */
	uint64_t out = 0, j = -7;
	do {
		out |= (( ((uint64_t)(*buf)) & 0x7f) << (j += 7) );
	} while( ((uint64_t)(*buf++)) & 0x80 );
	return (int64_t)( out & (1) ? ~(out >> 1) : (out >> 1) );
}

int bspatch(const uint8_t* old, int64_t oldsize, uint8_t* new, int64_t newsize, struct bspatch_stream* stream)
{
	uint8_t buf[10], *ptr;
	int64_t oldpos,newpos;
	int64_t ctrl[3];
	int64_t i;

	oldpos=0;newpos=0;
	while(newpos<newsize) {
		/* Read control data */
		for(i=0;i<=2;i++) {
			ptr = &buf[0];
			for(;;) {
				if( stream->read(stream, ptr, 1) ) {
					return -1;
				}
				if( ((*ptr++) & 0x80) == 0x00 ) {
					break;
				}
			}
			ctrl[i]=offtin(buf);
		};

		/* Sanity-check */
		if(newpos+ctrl[0]>newsize)
			return -1;

		/* Read diff string */
		if (stream->read(stream, new + newpos, ctrl[0]))
			return -1;

		/* Add old data to diff string */
		for(i=0;i<ctrl[0];i++)
			if((oldpos+i>=0) && (oldpos+i<oldsize))
				new[newpos+i]+=old[oldpos+i];

		/* Adjust pointers */
		newpos+=ctrl[0];
		oldpos+=ctrl[0];

		/* Sanity-check */
		if(newpos+ctrl[1]>newsize)
			return -1;

		/* Read extra string */
		if (stream->read(stream, new + newpos, ctrl[1]))
			return -1;

		/* Adjust pointers */
		newpos+=ctrl[1];
		oldpos+=ctrl[2];
	};

	return 0;
}

#if defined(BSPATCH_EXECUTABLE)

#include <bzlib.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <err.h>
#include <unistd.h>
#include <fcntl.h>

static int bz2_read(const struct bspatch_stream* stream, void* buffer, int length)
{
	int n;
	int bz2err;
	BZFILE* bz2;

	bz2 = (BZFILE*)stream->opaque;
	n = BZ2_bzRead(&bz2err, bz2, buffer, length);
	if (n != length)
		return -1;

	return 0;
}

int main(int argc,char * argv[])
{
	FILE * f;
	int fd;
	int bz2err;
	uint8_t header[24];
	uint8_t *old, *new;
	int64_t oldsize, newsize;
	BZFILE* bz2;
	struct bspatch_stream stream;

	if(argc!=4) errx(1,"usage: %s oldfile newfile patchfile\n",argv[0]);

	/* Open patch file */
	if ((f = fopen(argv[3], "r")) == NULL)
		err(1, "fopen(%s)", argv[3]);

	/* Read header */
	if (fread(header, 1, 16, f) != 16) {
		if (feof(f))
			errx(1, "Corrupt patch\n");
		err(1, "fread(%s)", argv[3]);
	}

	/* Check for appropriate magic */
	if (memcmp(header, "ENDSLEY/BSDIFF43", 16) != 0)
		errx(1, "Corrupt patch\n");

	/* Read lengths from header */
	newsize=offtin(header+16);
	if(newsize<0)
		errx(1,"Corrupt patch\n");

	/* Close patch file and re-open it via libbzip2 at the right places */
	if(((fd=open(argv[1],O_RDONLY,0))<0) ||
		((oldsize=lseek(fd,0,SEEK_END))==-1) ||
		((old=malloc(oldsize+1))==NULL) ||
		(lseek(fd,0,SEEK_SET)!=0) ||
		(read(fd,old,oldsize)!=oldsize) ||
		(close(fd)==-1)) err(1,"%s",argv[1]);
	if((new=malloc(newsize+1))==NULL) err(1,NULL);

	if (NULL == (bz2 = BZ2_bzReadOpen(&bz2err, f, 0, 0, NULL, 0)))
		errx(1, "BZ2_bzReadOpen, bz2err=%d", bz2err);

	stream.read = bz2_read;
	stream.opaque = bz2;
	if (bspatch(old, oldsize, new, newsize, &stream))
		errx(1, "bspatch");

	/* Clean up the bzip2 reads */
	BZ2_bzReadClose(&bz2err, bz2);
	fclose(f);

	/* Write the new file */
	if(((fd=open(argv[2],O_CREAT|O_TRUNC|O_WRONLY,0666))<0) ||
		(write(fd,new,newsize)!=newsize) || (close(fd)==-1))
		err(1,"%s",argv[2]);

	free(new);
	free(old);

	return 0;
}

#endif
#endif

#undef new

#include <string>
#include <utility>

#include "collage.hpp"

namespace collage {

	namespace {
		/* variable length encoding */
		std::string vlebit( size_t i ) {
			std::string out;
			do {
				out += (unsigned char)( 0x80 | (i & 0x7f));
				i >>= 7;
			} while( i > 0 );
			*out.rbegin() ^= 0x80;
			return out;
		}
		size_t vlebit( const char *&i ) {
			size_t out = 0, j = -7;
			do {
				out |= ((size_t(*i) & 0x7f) << (j += 7) );
			} while( size_t(*i++) & 0x80 );
			return out;
		}

		/* wrappers */
		static int bs_write(struct bsdiff_stream* stream, const void* buffer, int size) {
			return *((std::string*)stream->opaque) += std::string((char *)buffer,size), 0;
		}

		static int bs_read(const struct bspatch_stream * stream, void* buffer, int size) {
			std::pair<const char *,const char *> &opaque = *( (std::pair<const char *,const char *> *)stream->opaque );
			if( opaque.first + size <= opaque.second ) {
				memcpy( buffer, opaque.first, size );
				opaque.first += size;
				return 0;
			} else {
				return -1;
			}
		}
	}

	bool diff( std::string &result, const char *from0, const char *from1, const char *to0, const char *to1, unsigned Q ) {
		result = std::string();
		switch( Q ) {
			default:
			case BSDIFF: {
				struct bsdiff_stream diff_stream;
				diff_stream.opaque = (void *)&result;
				diff_stream.malloc = malloc;
				diff_stream.free = free;
				diff_stream.write = bs_write;
				if( 0 == bsdiff( (const uint8_t *)from0, from1 - from0, (const uint8_t *)to0, to1 - to0, &diff_stream ) ) {
					return true;
				}
			}
		}
		return false;
	}

	/* result must be resized in advance */
	bool patch( std::string &result, const char *from0, const char *from1, const char *diff0, const char *diff1, unsigned Q ) {
		switch( Q ) {
			default:
			case BSDIFF: {
				std::pair<const char *,const char *> pair( diff0, diff1 );
				struct bspatch_stream patch_stream;
				patch_stream.opaque = (void *)&pair;
				patch_stream.read = bs_read;
				if( 0 == bspatch( (const uint8_t *)from0, from1 - from0, (uint8_t *)(&result[0]), result.size(), &patch_stream ) ) {
					return true;
				}
			}
		}
		return false;
	}

	std::string diff( const std::string &from, const std::string &to, unsigned Q ) {
		std::string result;
		if( diff( result, from.c_str(), from.c_str() + from.size(), to.c_str(), to.c_str() + to.size(), Q ) ) {
			return std::string() + char(Q & 0x7F) + vlebit(to.size()) + result;
		} else {
			return std::string();
		}
	}

	std::string patch( const std::string &from, const std::string &diff ) {
		const char *diff8 = diff.c_str();
		unsigned Q = *diff8++;
		std::string result( vlebit(diff8), '\0' );
		if( patch( result, from.c_str(), from.c_str() + from.size(), diff8, diff8 + (diff.size() - (diff8 - diff.c_str()) ), Q ) ) {
			return result;
		} else {
			return std::string();
		}
	}
}

