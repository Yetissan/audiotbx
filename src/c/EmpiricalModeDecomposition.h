#pragma once

/*
 This code implements empirical mode decomposition in C.
 Required paramters include:
 - order: the number of IMFs to return
 - iterations: the number of iterations per IMF
 - locality: in samples, the nearest two extrema may be
 If it is not specified, there is no limit (locality = 0).

 Typical use consists of calling emdCreate(), followed by
 emdDecompose(), and then using the struct's "imfs" field
 to retrieve the data. Call emdClear() to deallocate memory
 inside the struct.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// MSVC needs "__inline" instead of "inline"
#if defined( _MSC_VER ) && !defined( __cplusplus )
# define inline __inline
#endif

#define cnew(type, size) ((type*) malloc((size) * sizeof(type)))
#define cdelete(ptr) free(ptr)

typedef struct {
	int iterations, order, locality;
	int *minPoints, *maxPoints;
	double *min, *max, **imfs, *residue;
	int size, minSize, maxSize;
} emdData;

void emdSetup(emdData* emd, int order, int iterations, int locality);
void emdResize(emdData* emd, int size);
void emdCreate(emdData* emd, int size, int order, int iterations, int locality);
void emdClear(emdData* emd);
void emdDecompose(emdData* emd, const double* signal);
void emdMakeExtrema(emdData* emd, const double* curImf);
void emdInterpolate(emdData* emd, const double* in, double* out, int* points, int pointsSize);
void emdUpdateImf(emdData* emd, double* imf);
void emdMakeResidue(emdData* emd, const double* cur);
inline int mirrorIndex(int i, int size);
